#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <pthread.h>
#include <unistd.h>
#include <float.h>
#include <math.h>

#include <assert.h>

#include "kmacros.h"
#include "rfmix.h"
#include "mm.h"

extern rfmix_opts_t rfmix_opts;

#define SAMPLES_PER_BLOCK (8)
typedef struct {
  int next_sample;
  int samples_completed;
  input_t *input;

  pthread_mutex_t lock;
} thread_args_t;

#define IDX(window,subpop) ((window)*n_subpops + (subpop))
static double viterbi(sample_t *sample, int haplotype, crf_window_t *crf_windows,
		    int n_windows, int n_subpops, mm *ma) {
  int i, j, k;

  /* Decode the estimated probabilities from random forest and cache in an array
     as doubles. If we do not use the IDX macro to do the indexing, we will typically
     use as much or more memory setting up the pointers as we will storing data */
  double *p = (double *) ma->allocate(sizeof(double)*n_windows*n_subpops, WHEREFROM);
  for(i=0; i < n_windows; i++) {
    for(k=0; k < n_subpops; k++) {
      double tmp = DF8(sample->est_p[haplotype][i][k]);
      if (tmp <= 0.) tmp = 0.0001;
      p[IDX(i,k)] = log(tmp);
    }
  }

  /* We may want to have a vector with the sample that stores their current global
     ancestry proportions and use that for the initial probability vector. For now
     set to equal probability. */
  double initial_p[n_subpops];
  for(k=0; k < n_subpops; k++) initial_p[k] = 1./(double) n_subpops;

  int *phi = (int *) ma->allocate(sizeof(int)*n_subpops*n_windows, WHEREFROM);
  double *d = (double *) ma->allocate(sizeof(double)*n_subpops, WHEREFROM);
  double *nd = (double *) ma->allocate(sizeof(double)*n_subpops, WHEREFROM);
  double *swap_d;
  double log_change, log_stay, max_d;
  int max_state;
  
  for(k=0; k < n_subpops; k++)
    d[k] = log(initial_p[k]) + p[ IDX(0,k) ];

  for(i=1; i < n_windows; i++) {
    double genetic_distance = crf_windows[i].genetic_pos - crf_windows[i-1].genetic_pos;
    if (genetic_distance > 0.5) genetic_distance = 0.5;
    log_change = log(genetic_distance);
    log_stay = log(1.0 - genetic_distance);
      
    for(j=0; j < n_subpops; j++) {
      double p_obs = p[ IDX(i,j) ];

      max_state = -1; max_d = -DBL_MAX;
      for(k=0; k < n_subpops; k++) {
	double tmp_d = d[k] + p_obs + ( (j==k) ? log_stay : log_change );
	if (tmp_d > max_d) {
	  max_state = k;
	  max_d = tmp_d;
	}
      }

      nd[j] = max_d;
      phi[ IDX(i,j) ] = max_state;
    }

    swap_d = d;
    d = nd;
    nd = swap_d;
  }

  int8_t *msp = sample->msp[haplotype];
  max_state = 0;
  for(k=1; k < n_subpops; k++)
    if (d[k] > d[max_state]) max_state = k;
  double logl = d[max_state];

  i = n_windows - 1;
  msp[i] = max_state;
  while(i > 0) {
    max_state = phi[ IDX(i, max_state) ];
    msp[i-1] = max_state;
    i--;
  }

  return logl;
}

static void forward_backward(sample_t *sample, int haplotype, crf_window_t *crf_windows, int n_windows, int n_subpops, mm *ma) {
  int i, j, k;
  double change, stay;

  if (haplotype > 1) return;
  
}


static void *crf_thread(void *targ) {
  thread_args_t *args;
  input_t *input;
  int start_sample, end_sample;
  mm *ma = new mm(16, WHEREFROM);
  
  args = (thread_args_t *) targ;
  input = args->input;
  
  pthread_mutex_lock(&args->lock);
  for(;;) {
    start_sample = args->next_sample;
    if (args->next_sample >= input->n_samples) break;
    
    end_sample = start_sample + SAMPLES_PER_BLOCK;
    if (end_sample > input->n_samples) end_sample = input->n_samples;
    args->next_sample = end_sample;
    pthread_mutex_unlock(&args->lock);

    for(int i=start_sample; i < end_sample; i++) {
      if (input->samples[i].apriori_subpop != 0) continue;
      for(int h=0; h < 4; h++) {
	viterbi(input->samples + i, h, input->crf_windows, input->n_windows,
		input->n_subpops-1, ma);
	forward_backward(input->samples + i, h, input->crf_windows, input->n_windows,
			 input->n_subpops-1, ma);
      }
    }

    pthread_mutex_lock(&args->lock);
    args->samples_completed += end_sample - start_sample;
    if (isatty(2)) {
      fprintf(stderr,"\rConditional random field ...       %6d/%6d (%1.1f%%)    ",
	      args->samples_completed, input->n_samples,
	      args->samples_completed / (double) input->n_samples * 100.);
    }
  }
  pthread_mutex_unlock(&args->lock);

  return NULL;
}

void crf(input_t *input) {
  thread_args_t args;
  pthread_t threads[rfmix_opts.n_threads];
  
  args.next_sample = 0;
  args.samples_completed = 0;
  args.input = input;
  pthread_mutex_init(&args.lock, NULL);

  for(int i=0; i < rfmix_opts.n_threads; i++)
    pthread_create(threads + i, NULL, crf_thread, (void *) &args);

  for(int i=0; i < rfmix_opts.n_threads; i++)
    pthread_join(threads[i], NULL);
  fprintf(stderr,"\n");

#if 0
  for(int i=0; i < input->n_samples; i++) {
    if (input->samples[i].apriori_subpop !=0) continue;
    for(int j=0; j < input->n_windows; j++) {
      fprintf(stderr,"sample %20.20s window %6d -", input->samples[i].sample_id, j);
      for(int h=0; h < 4; h++)
	fprintf(stderr,"\t%d", input->samples[i].msp[h][j]);
      fprintf(stderr,"\n");
    }
  }
#endif
#if 0
  for(int i=0; i < input->n_samples; i++) {
    if (input->samples[i].apriori_subpop !=0) continue;
    for(int h=0; h < 2; h++) {
      for(int j=0; j < input->n_windows; j++) {
	fprintf(stderr,"sample %20.20s haplotype %d window %6d -", input->samples[i].sample_id, h, j);
	for(int k=0; k < input->n_subpops; k++) 
	  fprintf(stderr,"\t%1.3f", DF8(input->samples[i].current_p[h][i][k]));
	fprintf(stderr,"\n");
      }
    }
  }
#endif
}

  
