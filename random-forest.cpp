#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <pthread.h>
#include <unistd.h>

#include "kmacros.h"
#include "rfmix.h"
#include "mm.h"

extern rfmix_opts_t rfmix_opts;

#define THREAD_WINDOW_CHUNK_SIZE (16)
typedef struct {
  input_t *input;
  int next_window;
  int windows_complete;

  pthread_mutex_t lock;
} thread_args_t;

typedef struct {
  int sample_idx;
  int *haplotype[4];
  double *est_p[4];
} wsample_t;

typedef struct {
  int n_snps;
  snp_t *snps;

  int n_query_samples;
  wsample_t *query_samples;

  int n_ref_haplotypes;
  int **ref_haplotypes;
  double **current_p;
} window_t;

/* In version 2 of RFMIX, an explicit tree structure is built classifying the
   reference haplotypes only, and then the query haplotypes are evaluated on
   this tree. This is in contrast to RFMIX version 1 never actually building
   an explicit tree but having an ephemeral one exist on the recursive call
   stack which is destroyed as the stack unwinds. This is to allow, ultimately,
   pre-training of the algorithm and storing the trees in serialized form in
   a file. In version 2 however I do not have plans to actually implement this,
   only to allow adaptation to support it in the future. */
typedef struct {
  /* snp_id is the index of the SNP that divides decendents into left or right.
     0 goes left, 1 goes right. */
  int snp_id;
  struct node *left;
  struct node *right;

  /* This is the sum of the current_p for each remaining reference haplotype.
     p is never normalized to sum to one, not even at node termination. Thus,
     if on this tree 10 reference haplotypes support assignment of a query
     haplotype falling at the terminal node for this tree, it counts 10 times
     as much as another tree for which only one reference haplotype is 
     supporting subpop assignment at the terminal node. Note that each query
     haplotype evaluated on the tree can only fall at exactly one terminal
     node of the tree */
  double *p;
} node_t;

typedef struct {
  int **haplotypes;
  int n_haplotypes;

  double *p;
  node_t *root;

  /* These fields are needed during building of the tree, but are not required
     once the tree is built. */
  int n_remain;
  int *remain_q;
  int n_snps;
  int *snplist;
} tree_t;
  
/* information that does not change as the tree progresses is in tree_t, except for the
   est_p field which is the ultimate result of the tree. est_p is merged into est_p in 
   wsample_t when grow_tree returns from the top most call. The information that changes
   on every decent call is stored in node_t structs. It is not implemented as a class
   for the reasons given above, and as a result there is some manual work to create new
   nodes.*/
static void grow_tree(tree_t *tree, node_t *node) {

}

/* Unpacks alleles from sample_t haplotypes into a copy here in ints (for speed) and
   also creates the two additional flipped phase haplotypes to analyze */
static void setup_query_sample(wsample_t *wsample, sample_t *sample, int n_subpops,
			       int snp_start, int snp_end, int crf_snp,
			       mm *ma) {
  int k, t, s;

  for(k=0; k < 4; k++)
    wsample->haplotype[k] = (int *) ma->allocate(sizeof(int)*(snp_end - snp_start + 1), WHEREFROM);
  
  /* Haplotypes 0 and 1 in wsample are phased as given in input->sample[i] */
  for(k=0; k < 2; k++) {
    for(s=snp_start,t=0; s <= snp_end; s++,t++)
      wsample->haplotype[k][t] = sample->haplotype[k][s];
  }

  /* Haplotypes 2 and 3 in wsample are flipped phase at the crf_window snp_idx */
  for(k=0; k < 2; k++) {
    for(s=snp_start,t=0; s < crf_snp; s++,t++)
      wsample->haplotype[k+2][t] = sample->haplotype[k][s];

    /* second half of window is reverse phase */
    int l = k == 0 ? 1 : 0;
    for(s=crf_snp; s <= snp_end; s++,t++) 
      wsample->haplotype[k+2][t] = sample->haplotype[l][s];
  }
  
  for(int j=0; j < 4; j++) {
    for(k=0; k < n_subpops-1; k++)
      wsample->est_p[j][k] = 0.;
  }
}

static void setup_ref_haplotype(int **ref_haplotypes, double **current_p, int n_subpops,
				int w, sample_t *sample, int snp_start, int snp_end,
				mm *ma) {
  int s, t;

  ref_haplotypes[0] = (int *) ma->allocate(sizeof(int)*(snp_end - snp_start + 1), WHEREFROM);
  ref_haplotypes[1] = (int *) ma->allocate(sizeof(int)*(snp_end - snp_start + 1), WHEREFROM);
  for(s=snp_start,t=0; s <= snp_end; s++,t++) {
    ref_haplotypes[0][t] = sample->haplotype[0][s];
    ref_haplotypes[1][t] = sample->haplotype[1][s];
  }
  
  for(int k=0; k < n_subpops-1; k++) {
    current_p[0][k] = DF8(sample->current_p[0][w][k]);
    current_p[1][k] = DF8(sample->current_p[1][w][k]);
  }
}

				
static void *random_forest_thread(void *targ) {
  thread_args_t *args = (thread_args_t *) targ;
  input_t *input = args->input;
  window_t window;
  int i;
  
  mm *ma = new mm(2, WHEREFROM);
  
  window.n_query_samples = 0;
  window.n_ref_haplotypes = 0;
  for(i=0; i < input->n_samples; i++) {
    if (input->samples[i].apriori_subpop == 0)
      window.n_query_samples++;
    else
      window.n_ref_haplotypes += 2;
  }

  /* This memory allocation will not need to be done and redone with every window.
     These fields will not change in length and can be reused. The haplotypes 
     themselves though may change with each window if the rf window size is 
     variable. Those are allocated in the loop using mm->allocate() */
  MA(window.query_samples, sizeof(wsample_t)*window.n_query_samples, wsample_t);
  for(i=0; i < window.n_query_samples; i++) {
    MA(window.query_samples[i].est_p[0], sizeof(double)*4*(input->n_subpops-1), double);
    for(int j=1; j < 4; j++)
      window.query_samples[i].est_p[j] = window.query_samples[i].est_p[j-1] + (input->n_subpops - 1);
  }
  MA(window.ref_haplotypes, sizeof(int *)*window.n_ref_haplotypes, int *);
  MA(window.current_p, sizeof(double *)*window.n_ref_haplotypes, double *);
  MA(window.current_p[0], sizeof(double)*window.n_ref_haplotypes*(input->n_subpops-1), double);
  for(i=1; i < window.n_ref_haplotypes; i++)
    window.current_p[i] = window.current_p[i-1] + (input->n_subpops - 1);
  
  /* args object is always locked at the loop start point or when loop exits */
  pthread_mutex_lock(&args->lock);
  for(;;) {
    
    /* Get the next chunk of windows to process and unlock the shared args object */
    int start_window = args->next_window;
    int end_window = start_window + THREAD_WINDOW_CHUNK_SIZE;
    if (end_window > input->n_windows) end_window = input->n_windows;
    args->next_window = end_window;
    pthread_mutex_unlock(&args->lock);

    for(int w=start_window; w < end_window; w++) {
      /* At the beginning of each window we can recycle all the memory that was allocated
	 for the previous window's needs. This is done effectively in one step by the
	 mm class. */
      ma->recycle();

      crf_window_t *crf = input->crf_windows + w;
      
    /* set up window_t object and decoded/unpacked information from the input_t object */
      window.n_snps = crf->rf_end_idx - crf->rf_start_idx + 1;
      window.snps = input->snps + crf->rf_start_idx;

      int q = 0;
      int r = 0;
      for(i=0; i < input->n_samples; i++) {
	if (input->samples[i].apriori_subpop == 0) {
	  window.query_samples[q].sample_idx = i;
	  setup_query_sample(window.query_samples + q, input->samples + i, input->n_subpops,
			     crf->rf_start_idx, crf->rf_end_idx, crf->snp_idx, ma);
	  q++;
	} else {
	  setup_ref_haplotype(window.ref_haplotypes + r, window.current_p + r, input->n_subpops,
			      w, input->samples + i, crf->rf_start_idx, crf->rf_end_idx, ma);
	  r += 2;
	}
      }

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

  delete ma;

  free(window.current_p[0]);
  free(window.current_p);
  free(window.ref_haplotypes);
  for(i=0; i < window.n_query_samples; i++) 
    free(window.query_samples[i].est_p[0]);
  free(window.query_samples);
  
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
