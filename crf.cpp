/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

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
extern int em_iteration;

static void normalize_vector(double *p, int n) {
  double sum_p = 0.;

  for(int k=0; k < n; k++) {
    if (p[k] < 0.) p[k] = 0.;
    sum_p += p[k];
  }
  for(int k=0; k < n; k++)
    p[k] /= sum_p;
}

static void __attribute__((unused))lnormalize_vector(double *p, int n) {
  double d = p[0];
  for(int i=1; i < n; i++)
    if (p[i] > d) d = p[i];

  double pd_sum = 0.;
  for(int i=0; i < n; i++) {
    p[i] = exp(p[i] - d);
    pd_sum += p[i];
  }
  
  for(int i=0; i < n; i++)
    p[i] = p[i]/pd_sum;  
}

typedef struct {
  int next_sample;
  int samples_completed;
  input_t *input;
  double viterbi_logl;
  double crf_weight;
  
  pthread_mutex_t lock;
} thread_args_t;

static void compute_state_change(double *r_stay, double *r_change, double d, int g, double w) {
  if (d < MINIMUM_GENETIC_DISTANCE) d = MINIMUM_GENETIC_DISTANCE;
  double r = 1.0 - exp(-d*(g-1));
  double nr = 1.0 - r;
  r = pow(r, w);

  *r_stay = nr;
  *r_change = r;
}

static double viterbi(sample_t *sample, int haplotype, crf_window_t *crf_windows,
		      int n_windows, int n_subpops, snp_t *snps, double w, mm *ma) {
  int i, j, k;

  /* Decode the estimated probabilities from random forest and cache in an array
     as doubles. If we do not use the IDX macro to do the indexing, we will typically
     use as much or more memory setting up the pointers as we will storing data */
  double *p = (double *) ma->allocate(sizeof(double)*n_windows*n_subpops, WHEREFROM);
  for(i=0; i < n_windows; i++) {
    for(k=0; k < n_subpops; k++) {
      double tmp = DF16(sample->est_p[haplotype][ IDX(i,k) ]);
      if (tmp <= 0.) tmp = 0.000001;
      p[IDX(i,k)] = log(tmp);
    }
  }

  /* We may want to have a vector with the sample that stores their current global
     ancestry proportions and use that for the initial probability vector. For now
     set to equal probability. */
  double initial_p[n_subpops];
  if (em_iteration > 1) {
    for(k=0; k < n_subpops; k++)
      initial_p[k] = 1.;
    for(i=0; i < n_windows; i++)
      initial_p[sample->msp[haplotype][i]]++;
    normalize_vector(initial_p, n_subpops);
    for(k=0; k < n_subpops; k++)
      initial_p[k] = log(initial_p[k]);
  } else {
    for(k=0; k < n_subpops; k++) {
      initial_p[k] = log(1./(double) n_subpops);
    }
  }
  
  int *phi = (int *) ma->allocate(sizeof(int)*n_subpops*n_windows, WHEREFROM);
  double *d = (double *) ma->allocate(sizeof(double)*n_subpops, WHEREFROM);
  double *nd = (double *) ma->allocate(sizeof(double)*n_subpops, WHEREFROM);
  double *swap_d;
  double max_d;
  int max_state;
  
  for(k=0; k < n_subpops; k++)
    d[k] = initial_p[k] + p[ IDX(0,k) ];

  for(i=1; i < n_windows; i++) {
    double stay, change;
    
    compute_state_change(&stay, &change,
			 crf_windows[i].genetic_pos - crf_windows[i-1].genetic_pos,
			 rfmix_opts.n_generations,w);
    double log_stay = (double) log(stay);
    double log_change = (double) log(change);
    
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

  /* Do not accept the viterbi results if the log likelihood is worse than last
     time. The maximum state path is left as it already is and the previous log
     likelihood returned */
  if (em_iteration > 0 && logl < sample->logl[haplotype]) return sample->logl[haplotype];
  
  i = n_windows - 1;
  msp[i] = max_state;
  while(i > 0) {
    max_state = phi[ IDX(i, max_state) ];
    msp[i-1] = max_state;
    i--;
  }

  sample->logl[haplotype] = logl;
  return logl;
}

//#define DEBUG
static void forward_backward(sample_t *sample, int haplotype, crf_window_t *crf_windows,
			     int n_windows, int n_subpops, double w, mm *ma) {
  int i, j, k;
  double change, stay;
  
  if (haplotype > 1) return;

  double *alpha = (double *) ma->allocate(sizeof(double)*n_subpops*n_windows, WHEREFROM);
  for(k=0; k < n_subpops; k++)
    alpha[IDX(0,k)] = DF16(sample->est_p[haplotype][ IDX(0,k) ]);
  normalize_vector(alpha, n_subpops);

  for(i=1; i < n_windows; i++) {
    double gd = crf_windows[i].genetic_pos - crf_windows[i-1].genetic_pos;
    compute_state_change(&stay, &change, gd, rfmix_opts.n_generations, w);
    for(j=0; j < n_subpops; j++) {
      alpha[ IDX(i,j) ] = 0.;
      for(k=0; k < n_subpops; k++)
	alpha[ IDX(i,j) ] += alpha[ IDX(i-1,k) ]*( (j==k) ? stay : change );
      alpha[ IDX(i,j) ] = alpha[ IDX(i,j) ]*DF16(sample->est_p[haplotype][ IDX(i,j) ]);
    }
    normalize_vector(alpha + i*n_subpops, n_subpops);
#ifdef DEBUG
    fprintf(stderr,"sample %s  haplotype %d   window %5d", sample->sample_id, haplotype, i);
    for(k=0; k < n_subpops; k++)
      fprintf(stderr,"\t%1.3f %1.5f",alpha[ IDX(i,k) ], DF16(sample->est_p[haplotype][ IDX(i,k) ]));
    fprintf(stderr,"\n");
#endif      
  }

  double *beta = (double *) ma->allocate(sizeof(double)*n_subpops*(n_windows), WHEREFROM);
  for(k=0; k < n_subpops; k++)
    beta[ IDX(n_windows-1,k) ] = 1.;
  normalize_vector(beta + (n_windows-1)*n_subpops, n_subpops);
  
  for(i=n_windows-2; i >=0 ; i--) {
    double gd = crf_windows[i+1].genetic_pos - crf_windows[i].genetic_pos;
    compute_state_change(&stay, &change, gd, rfmix_opts.n_generations, w);

    for(j=0; j < n_subpops; j++) {
      beta[ IDX(i,j) ] = 0.;
      for(k=0; k < n_subpops; k++) {
	beta[ IDX(i,j) ] +=  beta[ IDX(i+1,k) ] * DF16(sample->est_p[haplotype][ IDX(i+1,k) ]) *
	  ( (j==k) ? stay : change );
      }
      beta[ IDX(i,j) ] = beta[ IDX(i,j) ];
    }
    normalize_vector(beta + i*n_subpops, n_subpops);
#ifdef DEBUG
    fprintf(stderr,"sample %s  haplotype %d   window %5d pos %8.5f", sample->sample_id, haplotype, i,
	    crf_windows[i].genetic_pos*100.);
    for(k=0; k < n_subpops; k++)
      fprintf(stderr,"\t%1.3f",beta[ IDX(i,k) ]);
    fprintf(stderr,"\n");
#endif
  }

  for(i=0; i < n_windows; i++) {
    double p[n_subpops];
    for(k=0; k < n_subpops; k++) {
      p[k] = alpha[ IDX(i,k) ]*beta[ IDX(i,k) ];
    }
    normalize_vector(p, n_subpops);

    for(k=0; k < n_subpops; k++) {
      //      p[k] /= sum_p;
      sample->current_p[haplotype][ IDX(i,k) ] = ef16(p[k]);
    }
  }

  /* Calculate and set Suyash stay-in-state forward-backward probabilities */
  for(i=0; i < n_windows; i++) {
    double s = 0.; // stay-in-state
    double x = 0.; // change state
    for(int k=0; k < n_subpops; k++) {
      for(int l=0; l < n_subpops; l++) {
	if (k == l) {
	  s +=  alpha[ IDX(i,k) ] * beta[ IDX(i,l) ];
	} else {
	  x += alpha[ IDX(i,k) ] * beta[ IDX(i,l) ];
	}
      }
    }
    /* s + x should equal 1 (or floating point rounding error close enough), 
       but I'll normalize here anyway in case normalization of alpha and
       beta changes above */
    sample->sis_p[haplotype][i] = s/(s+x);
  }
}

static void *crf_thread(void *targ) {
  thread_args_t *args;
  input_t *input;
  int start_sample, end_sample;
  mm *ma = new mm(64, WHEREFROM);
  double total_logl = 0.;
  double logl;
  
  args = (thread_args_t *) targ;
  input = args->input;
  
  pthread_mutex_lock(&args->lock);
  for(;;) {
    start_sample = args->next_sample;
    if (args->next_sample >= input->n_samples) break;
    
    end_sample = start_sample + CRF_SAMPLES_PER_BLOCK;
    if (end_sample > input->n_samples) end_sample = input->n_samples;
    args->next_sample = end_sample;
    pthread_mutex_unlock(&args->lock);

    for(int i=start_sample; i < end_sample; i++) {
      /* During the internal simulation for finding optimum CRF weight, we are only
	 concerned with doing the CRF on the internal simulation samples */
      if (em_iteration == -1 && input->samples[i].s_sample != 1) continue;

      /* Skip the CRF for reference samples if --reanalyze-reference is not on */
      if (input->samples[i].apriori_subpop != -1 && rfmix_opts.reanalyze_reference == 0) continue;
      
      for(int h=0; h < 4; h++) {
	logl = viterbi(input->samples + i, h, input->crf_windows, input->n_windows,
		       input->n_subpops, input->snps, args->crf_weight, ma);
	if (h < 2) total_logl += logl;
	if (em_iteration != -1)
	  forward_backward(input->samples + i, h, input->crf_windows, input->n_windows,
			   input->n_subpops, args->crf_weight, ma);
      }
    }
    ma->recycle();

    pthread_mutex_lock(&args->lock);
    args->samples_completed += end_sample - start_sample;
    if (isatty(2)) {
      fprintf(stderr,"\rConditional random field ...       %6d/%6d (%1.1f%%)    ",
	      args->samples_completed, input->n_samples,
	      args->samples_completed / (double) input->n_samples * 100.);
    }
  }
  args->viterbi_logl += total_logl;
  pthread_mutex_unlock(&args->lock);

  delete ma;
  return NULL;
}

/* Note, does not show viterbi msp for haplotypes 2 and 3, the phase-flip windows */
static void __attribute__((unused))dump_results(input_t *input) {
  int n_subpops = input->n_subpops;

  for(int i=0; i < input->n_samples; i++) {
    if (input->samples[i].apriori_subpop != -1) continue;
    for(int h=0; h < 2; h++) {
      for(int j=0; j < input->n_windows; j++) {
	fprintf(stderr,"sample %20.20s haplotype %d window %6d - %d", input->samples[i].sample_id, h, j,
		input->samples[i].msp[h][j]);
	for(int k=0; k < input->n_subpops; k++) 
	  fprintf(stderr,"\t%1.5f", DF16(input->samples[i].current_p[h][ IDX(j,k) ]));
	fprintf(stderr,"\n");
      }
    }
  }
}

double crf(input_t *input, double w) {
  thread_args_t args;
  pthread_t threads[rfmix_opts.n_threads];

  args.next_sample = 0;
  args.samples_completed = 0;
  args.input = input;
  args.viterbi_logl = 0.;
  args.crf_weight = w;
 
  pthread_mutex_init(&args.lock, NULL);
  for(int i=0; i < rfmix_opts.n_threads; i++)
    pthread_create(threads + i, NULL, crf_thread, (void *) &args);

  for(int i=0; i < rfmix_opts.n_threads; i++)
    pthread_join(threads[i], NULL);
  
#if 0
  dump_results(input);
#endif

  return args.viterbi_logl;
}

  
