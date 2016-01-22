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

static void *crf_thread(void *targ) {
  thread_args_t *args;
  input_t *input;
  int start_sample, end_sample;
  
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
}

  
