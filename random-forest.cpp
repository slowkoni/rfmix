#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <pthread.h>
#include <unistd.h>

#include "kmacros.h"
#include "rfmix.h"

extern rfmix_opts_t rfmix_opts;

#define THREAD_WINDOW_CHUNK_SIZE (16)
typedef struct {
  input_t *input;
  int next_window;
  int windows_complete;

  pthread_mutex_t lock;
} thread_args_t;

static void *random_forest_thread(void *targ) {
  thread_args_t *args = (thread_args_t *) targ;
  input_t *input = args->input;

  /* args object is always locked at the loop start point or when loop exits */
  pthread_mutex_lock(&args->lock);
  for(;;) {
    /* Get the next chunk of windows to process and unlock the shared args object */
    int start_window = args->next_window;
    int end_window = start_window + THREAD_WINDOW_CHUNK_SIZE;
    if (end_window > input->n_windows) end_window = input->n_windows;
    args->next_window = end_window;
    pthread_mutex_unlock(&args->lock);

    for(int window=start_window; window < end_window; window++) {
      
    /* set up window_t object and decoded/unpacked information from the input_t object */

    /* Call window function to process the window */

    /* Repack window_t object results into input_t and reset/reinitialize window_t 
       for next window */

    }
    
    pthread_mutex_lock(&args->lock);
    args->windows_complete += end_window - start_window;
    if (isatty(2))
      fprintf(stderr, "\rGrowing Random Forest Trees -- (%d/%d) %5.1f%%   ", args->windows_complete,
	      input->n_windows, args->windows_complete / (double) input->n_windows * 100.);

    if (args->next_window >= input->n_windows) break;
  }
  pthread_mutex_unlock(&args->lock);

  return NULL;
}

void random_forest(input_t *input) {
  thread_args_t *args;
  MA(args, sizeof(thread_args_t), thread_args_t);
  args->input = input;
  args->next_window = 0;
  args->windows_complete = 0;
  pthread_mutex_init(&args->lock, NULL);
  
  pthread_t *threads;
  MA(threads, sizeof(pthread_t)*rfmix_opts.n_threads, pthread_t);

  for(int i=0; i < rfmix_opts.n_threads; i++)
    pthread_create(threads + i, NULL, random_forest_thread, (void *) args);

  for(int i=0; i < rfmix_opts.n_threads; i++)
    pthread_join(threads[i], NULL);
  fprintf(stderr,"\n");
  
  free(threads);
  free(args);
}
