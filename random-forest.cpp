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
#include "md5rng.h"
#include "rfmix.h"
#include "mm.h"

extern rfmix_opts_t rfmix_opts;
extern int em_iteration;

typedef struct {
  input_t *input;
  int next_window;
  int windows_complete;
  md5rng *rng;
  
  pthread_mutex_t lock;
} thread_args_t;

typedef struct {
  int sample_idx;
  int *haplotype[4];
  double *est_p[4];
  //double *current_p[4];
} wsample_t;

typedef struct {
  int *haplotype;
  double *current_p;
  int max_idx;
} ref_haplotype_t;

typedef struct {
  int idx;
  int rng_idx;
  int n_snps;
  snp_t *snps;;

  int n_query_samples;
  wsample_t *query_samples;

  int n_ref_haplotypes;
  ref_haplotype_t *ref_haplotypes;

  int n_subpops;
  int *n_ref_haplotypes_by_subpop;
  int **ref_haplotype_list;
} window_t;

/* In version 2 of RFMIX, an explicit tree structure is built classifying the
   reference haplotypes only, and then the query haplotypes are evaluated on
   this tree. This is in contrast to RFMIX version 1 never actually building
   an explicit tree but having an ephemeral one exist on the recursive call
   stack which is destroyed as the stack unwinds. This is to allow, ultimately,
   pre-training of the algorithm and storing the trees in serialized form in
   a file. In version 2 however I do not have plans to actually implement this,
   only to allow adaptation to support it in the future. */
typedef struct node {
  /* snp_id is the index of the SNP that divides decendents into left or right.
     0 goes left, 1 goes right. */
  int snp_id;
  int level;
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
  int n_subpops;
  md5rng *rng;
  int rng_idx;
  int window_idx;

  int n_snps;
  int **haplotypes;
  int n_haplotypes;
  double **current_p;

  node_t *root;
} tree_t;

static void __attribute__((unused))output_vector(FILE *f, double *p, int n, char delim) {
  fprintf(f,"%1.3f",p[0]);
  for(int i=1; i < n; i++)
    fprintf(f,"%c%1.3f", delim, p[i]);
}

static void normalize_vector(double *p, int n) {
  double p_sum = 0.;
  for(int i=0; i < n; i++)
    p_sum += p[i];
  for(int i=0; i < n; i++)
    p[i] /= p_sum;
}

static void output_node(FILE *f, node_t *node, int n) {

  if (node->snp_id == -1) {
    fprintf(f,"(-1,[");
    fprintf(f,"%1.1f", node->p[0]);
    for(int k=1; k < n-1; k++)
      fprintf(f,",%1.1f", node->p[k]);
    fprintf(f,"])");
  }
  else {
    fprintf(f,"(%d,", node->snp_id);
    output_node(f, node->left, n);
    fprintf(f,",");
    output_node(f, node->right, n);
    fprintf(f,")");
  }
}

static void __attribute__((unused))output_tree(FILE *f, tree_t *tree) {
  output_node(f, tree->root, tree->n_subpops);
  fprintf(f,"\n\n");
}

static inline void add_current_p(double *p, double *current_p, int n) {
  for(int k=0; k < n; k++) p[k] += current_p[k];
}

/* Normalizes a probability vector and returns the shannon information - any
   debug code outputting the value should divide by M_LN2 to report in bits. I
   am not calling normalize_vector() above because code calling this function
   does not want p[] modified. */
static double shannon_information(double *p, int np, int n) {
 double sum_p = 0.;

  for(int k=0; k < n; k++) 
    sum_p += p[k];

  double si = 0.;
  for(int k=0; k < n; k++) {
    double tmp = p[k]/sum_p;
    if (tmp > 0.) si += tmp*log(tmp);
  }

  return -si*n;    
}

/* Return the "gini index" of a probabily vector. The vector typically needs to
   be normalized first but calling code does not want p[] modified so not calling
   normalize_vector() above. The gini index is a measure of how different p[] is
   from equal probability among componenets. From what I understand, minimizing
   the gini index in choosing decisions for the tree is less susceptible to bias
   toward selecting attributes which uniquely identify single or small number of
   items among a large number of categories than shannon information */
static double gini_index(double *p, int n) {
  double sum_p = 0.;

  for(int k=0; k < n; k++) 
    sum_p += p[k];
  
  double gi = 0.;
  for(int i=0; i < n; i++) {
    gi += p[i]/sum_p * p[i]/sum_p;
  }

  return 1.0 - gi;
}
  
/* calculates the sum of the shannon information content of the resulting child 
   nodes if we split the reference haplotypes remaining <ref_q> on snp <snp>. This
   function and its caller can also use the gini index for this purpose */
static double evaluate_snp(double *si, double *n_child, tree_t *tree, int snp, int *ref_q, int n_ref, mm *ma) {
  int j,k;
  double p[2][tree->n_subpops];
  
  for(k=0; k < tree->n_subpops; k++) {
    p[0][k] = 0.;
    p[1][k] = 0.;
  }
  n_child[0] = 0.;
  n_child[1] = 0.;
  
  for(j=0; j < n_ref; j++) {
    if (tree->haplotypes[ref_q[j]][snp] == 0) {
      add_current_p(p[0], tree->current_p[ref_q[j]], tree->n_subpops);
      n_child[0] += 1.0;
    }
    else if (tree->haplotypes[ref_q[j]][snp] == 1) {
      add_current_p(p[1], tree->current_p[ref_q[j]], tree->n_subpops);
      n_child[1] += 1.0;
    }
    else {
      /* Missing data code - perhaps we should test explicitly for it? */
      /* We handle missing data by sending the haplotypes having missing data down
	 to both child nodes */
      add_current_p(p[0], tree->current_p[ref_q[j]], tree->n_subpops);
      add_current_p(p[1], tree->current_p[ref_q[j]], tree->n_subpops);
      n_child[0] += 0.5;
      n_child[1] += 0.5;
    }
  }

  si[0] = gini_index(p[0], tree->n_subpops)*n_child[0];
  si[1] = gini_index(p[1], tree->n_subpops)*n_child[1];
  //si[0] = shannon_information(p[0], n_child[0], tree->n_subpops);
  //si[1] = shannon_information(p[1], n_child[1], tree->n_subpops);

  //return sq_difference(p[0], p[1], tree->n_subpops);
  return (si[0] + si[1])/(n_child[0] + n_child[1]);
}

static node_t* add_node(tree_t *tree, int *snp_q, int n_snps, int *ref_q, int n_ref,
			double si, mm *ma, int level) {

  node_t *node = (node_t *) ma->allocate(sizeof(node_t), WHEREFROM);
  node->p = NULL;
  node->level = level;
  
  /* Terminate recursion if there is only 1 reference haplotype left, no snps left
     to divide on, or the information content is less than half a bit per haplotype */
  if (n_snps == 0 || n_ref <= 1 || level >= 20) {
#ifndef DEBUG_L1
    assert(n_ref != 0);
#endif
    node->p = (double *) ma->allocate(sizeof(double)*tree->n_subpops, WHEREFROM);
    for(int k=0; k < tree->n_subpops; k++)
      node->p[k] = 0.;
    for(int i=0; i < n_ref; i++) {
      for(int k=0; k < tree->n_subpops; k++)
	node->p[k] += tree->current_p[ref_q[i]][k];
    }
#ifdef DEBUG_L1
    /* This should not happen as the termination condition below should be triggered
       because of lack of a SNP to choose which still results in at least one haplotype
       in each branch. If we do get here, something is incorrect in the code below */
    if (n_ref == 0) fprintf(stderr,"N HAPLOTYPES IS ZERO!\n");
    fprintf(stderr,"Force terminate tree at level %d - %6.1f  %3d snps  %3d haplotypes\n",
	    level, si, n_snps, n_ref);
    fprintf(stderr,"[ %4.1f", node->p[0]);
    for(int k=1; k < tree->n_subpops; k++)
      fprintf(stderr,", %4.1f",node->p[k]);
    fprintf(stderr," ]\n");
#endif
    
    node->left = NULL;
    node->right = NULL;
    node->snp_id = -1;
    return node;
  }

  /* In original version 1 RFMIX, the number of SNPs which are randomly evaluated from
     the total number remaining is some fraction, called the mtry factor, which defaults
     to 1/sqrt(2) of those SNPs remaining. This was a parameter of the original program
     but is currently fixed here. All SNPs should not be evaluated, or the only randomness
     in the tree is due to bootstrap selection of haplotypes */
  int n_try = (int) n_snps/2.;//(n_snps * M_SQRT1_2);
  //if (n_try < 10) n_try = n_snps;

  /* Randomly permute the array of SNPs, then evaluate the first n_try of them. The random
     number generator is key'd to the window index for repeatability of runs on the same
     input with multithreading enabled */
  for(int i=0; i < n_snps; i++) {
    int j = tree->rng->uniform_int(RFOREST_RNG_KEY, tree->window_idx, tree->rng_idx++, 0, n_snps);
    int tmp = snp_q[i];
    snp_q[i] = snp_q[j];
    snp_q[j] = tmp;
  }

  int best_snp = -1;
  double best_split = DBL_MAX;;
  double child_si[2];
  double best_si[2];
  double child_n[2];
  for(int i=0; i < n_try && i < n_snps; i++) {
    int snp = snp_q[i];

#define NODE_SIZE (2.0)
    double split_information = evaluate_snp(child_si, child_n, tree, snp, ref_q, n_ref, ma);
    if (split_information < best_split &&
	child_n[0] > rfmix_opts.node_size && child_n[1] > rfmix_opts.node_size) {
      //#define DEBUG_L2
#ifdef DEBUG_L2
      fprintf(stderr,"level %d - try %3d/%3d/%3d  split %6.3f in %6.3f,%6.3f  parent %6.3f  %6.3f\n", level,
	      i, n_try, n_snps, split_information, child_si[0], child_si[1],si, split_information - si);
#endif
      best_snp = i;
      best_si[0] = child_si[0];
      best_si[1] = child_si[1];
      best_split = split_information;      
    } else if (child_n[0] <= rfmix_opts.node_size || child_n[1] <= rfmix_opts.node_size ||
	       (fabs(child_n[0] - n_ref) < 0.1 && fabs(child_n[1] - n_ref) < 0.1)) {
      snp_q[i] = snp_q[n_snps-1];
      snp_q[n_snps-1] = snp;
      n_snps--;
      i--;
    }
  }
  
  /* If no SNP was selected, then no division is possible with the n_try SNPs evaluated such 
     that both branches have at least one haplotype, with the n_try SNPs evaluated. This
     terminates the tree at this node. */
  if (best_snp == -1) {
    node->p = (double *) ma->allocate(sizeof(double)*tree->n_subpops, WHEREFROM);
    for(int k=0; k < tree->n_subpops; k++)
      node->p[k] = 0.;
    for(int i=0; i < n_ref; i++) {
      for(int k=0; k < tree->n_subpops; k++)
	node->p[k] += tree->current_p[ref_q[i]][k];
    }

    //#define DEBUG_L2
#ifdef DEBUG_L2
    fprintf(stderr,"Terminate tree at level %d (%d snps remain)\n", level, n_snps);
    fprintf(stderr,"[ %4.1f", node->p[0]);
    for(int k=1; k < tree->n_subpops; k++)
      fprintf(stderr,", %4.1f",node->p[k]);
    fprintf(stderr," ]\n");
#endif
    //    normalize_vector(node->p, tree->n_subpops);
    
    node->left = NULL;
    node->right = NULL;
    node->snp_id = -1;
    return node;
  }

  /* Normal condition - we divide the haplotypes based on SNP snp_q[best_snp] */
  node->snp_id = snp_q[best_snp];
  
  int *child_ref_q[2];
  int child_n_ref[2];
  
  for(int j=0; j < 2; j++) child_n_ref[j] = 0;
  for(int i=0; i < n_ref; i++) {
    if (tree->haplotypes[ref_q[i]][snp_q[best_snp]] == 0) {
      child_n_ref[0]++;
    }
    else if (tree->haplotypes[ref_q[i]][snp_q[best_snp]] == 1) {
      child_n_ref[1]++;
    }
    else {
      child_n_ref[0]++;
      child_n_ref[1]++;
    }
  }

  for(int j=0; j < 2; j++) {
    child_ref_q[j] = (int *) ma->allocate(sizeof(int)*(child_n_ref[j]+1), WHEREFROM);
    child_n_ref[j] = 0;
  }
  for(int i=0; i < n_ref; i++) {
    if (tree->haplotypes[ref_q[i]][snp_q[best_snp]] == 0) {
      child_ref_q[0][child_n_ref[0]++] = ref_q[i];
    }
    else if (tree->haplotypes[ref_q[i]][snp_q[best_snp]] == 1) {
      child_ref_q[1][child_n_ref[1]++] = ref_q[i];
    }
    else {
      /* Missing data results in the haplotype going down both branches */
      child_ref_q[0][child_n_ref[0]++] = ref_q[i];
      child_ref_q[1][child_n_ref[1]++] = ref_q[i];
    }
  }

  /* Move the selected SNP to the end of the queue. We'll tell the child calls the 
     queue is one shorter. Note that this automatically makes the SNP selected at 
     this level available again for selection when this function returns, meaning
     higher levels in the tree that start going down the right branch can use it
     to divide the haplotypes that went that way. */
  int tmp = snp_q[best_snp];
  snp_q[best_snp] = snp_q[n_snps-1];
  snp_q[n_snps-1] = tmp;

  node->left = add_node(tree, snp_q, n_snps-1, child_ref_q[0], child_n_ref[0], best_si[0], ma, level+1);
  node->right = add_node(tree, snp_q, n_snps-1, child_ref_q[1], child_n_ref[1], best_si[1], ma, level+1);
  
  return node;
}

static void flat_bootstrap(tree_t *tree, window_t *window) {
  int i;
  for(i=0; i < window->n_ref_haplotypes; i++) {
    int j = tree->rng->uniform_int(RFOREST_RNG_KEY, window->idx, window->rng_idx++, 0, window->n_ref_haplotypes);

    tree->haplotypes[i] = window->ref_haplotypes[j].haplotype;
    tree->current_p[i] = window->ref_haplotypes[j].current_p;
  }
  tree->n_haplotypes = i;
}

static void hierarchical_bootstrap(tree_t *tree, window_t *window) {
  int i;
  
  tree->n_haplotypes = 0;
  for(i=0; i < window->n_ref_haplotypes; i++) {
    int k = tree->rng->uniform_int(RFOREST_RNG_KEY, window->idx, window->rng_idx++, 0, window->n_subpops);
    int n = window->n_ref_haplotypes_by_subpop[k];
    if (n == 0) { i--; continue; }
    
    int j = tree->rng->uniform_int(RFOREST_RNG_KEY, window->idx, window->rng_idx++, 0, n);
    int h = window->ref_haplotype_list[k][j];

    tree->haplotypes[i] = window->ref_haplotypes[h].haplotype;
    tree->current_p[i] = window->ref_haplotypes[h].current_p;
  }
  tree->n_haplotypes = i;
}

static void stratified_bootstrap(tree_t *tree, window_t *window) {

  int h = 0;
  for(int k=0; k < window->n_subpops; k++) {
    int n = window->n_ref_haplotypes_by_subpop[k];
    for(int i=0; i < n; i++) {
      int j = tree->rng->uniform_int(RFOREST_RNG_KEY, window->idx, window->rng_idx++, 0, n);
      int t = window->ref_haplotype_list[k][j];
      
      tree->haplotypes[h] = window->ref_haplotypes[t].haplotype;
      tree->current_p[h] = window->ref_haplotypes[t].current_p;
      h++;
    }
  }
  tree->n_haplotypes = h;
  
}

/* An essential part of the random forest method is that each tree has a "bootstrapped" random 
   selection, with replacement, of haplotypes. Some haplotypes may be represented 2 or more 
   times, others none, but which ones are different for each tree. This will produce different
   decisions in each tree, as well as that of the random evaluation of only some of the 
   available SNPs to make a division on. This is to reduce over-fitting to characteristics of 
   the reference data which are not real characteristics of the population they are sampled 
   from but rather artifacts of the random sampling process. This is like how the mean of a
   randomly drawn sample estimates, but typically does not equal, the mean of the entire
   population. */
static void bootstrap_haplotypes(input_t *input, tree_t *tree, window_t *window, mm *ma) {
  tree->haplotypes = (int **) ma->allocate(sizeof(int *)*window->n_ref_haplotypes, WHEREFROM);
  tree->current_p = (double **) ma->allocate(sizeof(double *)*window->n_ref_haplotypes, WHEREFROM);
  switch(rfmix_opts.bootstrap_mode) {
  case RF_BOOTSTRAP_FLAT:
    flat_bootstrap(tree, window);
    break;
  case RF_BOOTSTRAP_HIERARCHICAL:
    hierarchical_bootstrap(tree, window);
    break;
  case RF_BOOTSTRAP_STRATIFIED:
    stratified_bootstrap(tree, window);
    break;
  default:
    fprintf(stderr,"Runtime error: Unknown bootstrap mode %d in random forest\n", rfmix_opts.bootstrap_mode);
    break;
  }
}


static tree_t *build_tree(input_t *input, window_t *window, md5rng *rng, mm *ma) {
  int i;
  tree_t *tree = (tree_t *) ma->allocate(sizeof(tree_t), WHEREFROM);

  /* Because we are not going to pass the window to the functions that build the
     tree, these variables are going to be copied to the tree structure. All but
     n_subpops are not essential properties of the tree built */
  tree->window_idx = window->idx,
  tree->n_subpops = window->n_subpops;
  tree->rng = rng;
  tree->rng_idx = window->rng_idx;

  /* Randomly select, with replacement, reference haplotypes for this tree. See above. */
  bootstrap_haplotypes(input, tree, window, ma);
  int *ref_q = (int *) ma->allocate(sizeof(int)*tree->n_haplotypes, WHEREFROM);
  for(i=0; i < tree->n_haplotypes; i++)
    ref_q[i] = i;

  /* All SNPs in the RF window are candidates, only a random subsample are evaluated
     at each node. See add_node() above. */
  int *snp_q = (int *) ma->allocate(sizeof(int)*window->n_snps, WHEREFROM);
  for(i=0; i < window->n_snps; i++)
    snp_q[i] = i;

  /* Initialize the shannon information of all bootstrap-selected reference haplotypes
     present at the start (root node) of the tree, and build the tree */
  double p[window->n_subpops];
  for(int k=0; k < window->n_subpops; k++)
    p[k] = 0.;
  for(i=0; i < tree->n_haplotypes; i++) {
    for(int k=0; k < window->n_subpops; k++)
      p[k] += tree->current_p[i][k];
  }
  double si = shannon_information(p, tree->n_haplotypes, window->n_subpops);
  tree->root = add_node(tree, snp_q, window->n_snps, ref_q, tree->n_haplotypes,
			si, ma, 0);

  /* The random number generator sequence index is increased during tree building, copy
     the incremented value back to the window */
  window->rng_idx = tree->rng_idx;
  return tree;
}

static void evaluate_tree(double *p, int *d, int *haplotype, node_t *node, int n_subpops) {
  
  if (node->snp_id == -1) {
    for(int k=0; k < n_subpops; k++)
      p[k] += node->p[k];
    (*d)++;
    return;
  }
  
  int allele = haplotype[node->snp_id];
  if (allele == 0)
    evaluate_tree(p, d, haplotype, node->left, n_subpops);
  else if (allele == 1)
    evaluate_tree(p, d, haplotype, node->right, n_subpops);
  else if (allele == 2) {
    evaluate_tree(p, d, haplotype, node->left, n_subpops);
    evaluate_tree(p, d, haplotype, node->right, n_subpops);
  }
}

static void evaluate_sample(wsample_t *wsample, window_t *window, tree_t **trees, int n_trees) {
  
  for(int h=0; h < 4; h++) {
    for(int t=0; t < n_trees; t++) {
      double p[window->n_subpops];
      for(int k=0; k < window->n_subpops; k++) p[k] = 0.;
      
      /* This is tricky - because of missing data, some trees might add the terminal node p vector
	 of several nodes, not just one. d is incremented by evaluate_tree() at every terminal 
	 node that it adds to p[]. After the return of the call below, we then merge the p vector
	 into est_p[h][] and divide it by the number of terminal nodes that were used by the tree.
	 This is necessary otherwise trees where missing data results in multiple terminal nodes
	 being used would in effect have a higher weight where fewer (one ideally) terminal 
         nodes were reached. */
      int d = 0;
      evaluate_tree(p, &d, wsample->haplotype[h], trees[t]->root, window->n_subpops);
      
      for(int k=0; k < window->n_subpops; k++) {
	  wsample->est_p[h][k] += p[k]/(double) d;
      }
    }

    /* Normalize to probabilities that sum to one across all subpops. */
    normalize_vector(wsample->est_p[h], window->n_subpops);
#ifdef DEBUG_L2
    fprintf(stderr,"window %d sample %d haplotype %d - %4.2f", window->idx, wsample->sample_idx,
	    h, wsample->est_p[h][0]);
    for(int k=1; k < window->n_subpops; k++)
      fprintf(stderr, ", %4.2f", wsample->est_p[h][k]);
    fprintf(stderr,"\n");
#endif
  }
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
    //    wsample->current_p[j] = (double *) ma->allocate(sizeof(double)*n_subpops, WHEREFROM);
    for(k=0; k < n_subpops; k++) {
      wsample->est_p[j][k] = 0.0;
      //      wsample->current_p[j][k] = DF16(sample->current_p[j & 0x1][k]);
    }
  }
}

static void setup_ref_haplotypes(window_t *w, input_t *input, int start_snp, int end_snp,
				 mm *ma) {
  int i, k, h;

  int window_idx = w->idx;
  int n_samples = input->n_samples;
  int n_subpops = input->n_subpops;
  int n_snps = end_snp - start_snp + 1;
  sample_t *samples = input->samples;
  
  double *p_tmp = new double[input->n_subpops];
  
  ref_haplotype_t *rh = (ref_haplotype_t *) ma->allocate(sizeof(ref_haplotype_t)*input->n_samples*2, WHEREFROM);
  int *n_by_subpop = (int *) ma->allocate(sizeof(int)*n_subpops, WHEREFROM);
  for(k=0; k < n_subpops; k++) n_by_subpop[k] = 0;

  int nrh = 0;
  for(i=0; i < n_samples; i++) {
    /* Do not include the parents for the internal control simulation during
       learning for an optimum CRF weight */
    if (em_iteration == -1 && samples[i].s_parent == 1) continue;

    /* Never use an internally simulated sample as a reference sample */
    if (samples[i].s_sample == 1) continue;

    /* In the initial, non simulation control iteration, we have no estimates for
       ancestry on the query individuals, so skip them */
    if (em_iteration <= 0 && samples[i].apriori_subpop == -1) continue;

    /* At this point, we the sample may be included if its haplotypes have strong
       enough estimates - we check that for the individual haplotype */
    for(h=0; h < 2; h++) {

      // unpack the haplotype's current_p and find the subpop with the max probability
      int max = 0;
      for(k=0; k < n_subpops; k++) {
	p_tmp[k] = DF16(samples[i].current_p[h][ IDX(window_idx, k) ]);
	if (p_tmp[k] > p_tmp[max]) max = k;
      }

      // If this haplotype is suitable for use as reference, add it
      if (p_tmp[max] > P_MINIMUM_FOR_REF) {
	rh[nrh].haplotype = (int *) ma->allocate(sizeof(int)*n_snps, WHEREFROM);
	rh[nrh].current_p = (double *) ma->allocate(sizeof(double)*n_subpops, WHEREFROM);
	
	/* copy out the alleles */
	for(int s = start_snp, t=0; s <= end_snp; s++, t++)
	  rh[nrh].haplotype[t] = (int) samples[i].haplotype[h][s];

	if (em_iteration > 1) {
	/* copy over the current_p that we already unpacked */
	  for(int k=0; k < n_subpops; k++)
	    rh[nrh].current_p[k] = p_tmp[k];
	} else {
	  double d = 0.1/(2. + em_iteration);
	  for(int k=0; k < n_subpops; k++)
	    rh[nrh].current_p[k] = d/(n_subpops-1);
	  rh[nrh].current_p[ samples[i].msp[h][window_idx] ] = 1. - d;
	}
	rh[nrh].max_idx = max;
	
	n_by_subpop[max]++;
	nrh++;
      }
    }
  }

  delete[] p_tmp;
  
  int **subpop_rh_list = (int **) ma->allocate(sizeof(int *)*n_subpops, WHEREFROM);
  int c[n_subpops];
  for(k=0; k < n_subpops; k++) {
    subpop_rh_list[k] = (int *) ma->allocate(sizeof(int)*n_by_subpop[k], WHEREFROM);
    c[k] = 0;
  }
  for(i=0; i < nrh; i++) {
    k = rh[i].max_idx;
    subpop_rh_list[k][ c[k]++ ] = i;
  }
  
  w->ref_haplotypes = rh;
  w->n_ref_haplotypes = nrh;
  w->n_ref_haplotypes_by_subpop = n_by_subpop;
  w->ref_haplotype_list = subpop_rh_list;
}


static void *random_forest_thread(void *targ) {
  thread_args_t *args = (thread_args_t *) targ;
  input_t *input = args->input;
  int n_subpops = input->n_subpops;
  window_t window;
  int i;
  
  mm *ma = new mm(16, WHEREFROM);
  
  window.n_query_samples = 0;
  for(i=0; i < input->n_samples; i++) {
    if (rfmix_opts.reanalyze_reference == 0 && input->samples[i].apriori_subpop >= 0) continue;
    if (em_iteration == -1 && input->samples[i].s_sample != 1) continue;
    if (em_iteration != -1 && input->samples[i].s_sample == 1) continue;
    
    window.n_query_samples++;
  }

  /* This memory allocation will not need to be done and redone with every window.
     These fields will not change in length and can be reused. The haplotypes 
     themselves though may change with each window if the rf window size is 
     variable. Those are allocated in the loop using mm->allocate() */
  MA(window.query_samples, sizeof(wsample_t)*window.n_query_samples, wsample_t);
  for(i=0; i < window.n_query_samples; i++) {
    MA(window.query_samples[i].est_p[0], sizeof(double)*4*(n_subpops), double);
    for(int j=1; j < 4; j++)
      window.query_samples[i].est_p[j] = window.query_samples[i].est_p[j-1] + n_subpops;
  }
  
  /* args object is always locked at the loop start point or when loop exits */
  pthread_mutex_lock(&args->lock);
  for(;;) {
    
    /* Get the next chunk of windows to process and unlock the shared args object */
    int start_window = args->next_window;
    int end_window = start_window + RF_THREAD_WINDOW_CHUNK_SIZE;
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
      window.n_subpops = n_subpops;
      window.idx = w;
      window.rng_idx = 1;
      window.n_snps = crf->rf_end_idx - crf->rf_start_idx + 1;
      window.snps = input->snps + crf->rf_start_idx;
#ifdef DEBUG_L2
      fprintf(stderr,"window %d  %d snps   %d to %d\n", window.idx, window.n_snps, crf->rf_start_idx,
	      crf->rf_end_idx);
#endif
      int q = 0;
      for(i=0; i < input->n_samples; i++) {
	if (rfmix_opts.reanalyze_reference == 0 && input->samples[i].apriori_subpop >= 0) continue;
	if (em_iteration == -1 && input->samples[i].s_sample != 1) continue;
	if (em_iteration != -1 && input->samples[i].s_sample == 1) continue;
	
	window.query_samples[q].sample_idx = i;
	setup_query_sample(window.query_samples + q, input->samples + i, n_subpops,
			   crf->rf_start_idx, crf->rf_end_idx, crf->snp_idx, ma);
	q++;
      }
      setup_ref_haplotypes(&window, input, crf->rf_start_idx, crf->rf_end_idx, ma);

      /* Build trees */
      tree_t **trees = (tree_t **) ma->allocate(sizeof(tree_t *)*rfmix_opts.n_trees, WHEREFROM);
      for(i=0; i < rfmix_opts.n_trees; i++)
	trees[i] = build_tree(input, &window, args->rng, ma);

#if 0
      pthread_mutex_lock(&args->lock);
      for(i=0; i < rfmix_opts.n_trees; i++) 
	output_tree(stdout, trees[i]);
      pthread_mutex_unlock(&args->lock);
#endif
      
      /* evaluate trees */
      for(i=0; i < window.n_query_samples; i++)
	evaluate_sample(window.query_samples + i, &window, trees, rfmix_opts.n_trees);

	/* Repack window_t object results into input_t and reset/reinitialize window_t 
       for next window */
      for(i=0; i < window.n_query_samples; i++) {
	wsample_t *wsample = window.query_samples + i;

	for(int j=0; j < 4; j++) {

	  for(int k=0; k < n_subpops; k++) {
	    double p = wsample->est_p[j][k];
	    input->samples[wsample->sample_idx].est_p[j][ IDX(w,k) ] = ef16(p);
	  }
	}
      }
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

  for(i=0; i < window.n_query_samples; i++) 
    free(window.query_samples[i].est_p[0]);
  free(window.query_samples);
  
  return NULL;
}

static void __attribute__((unused))dump_results(input_t *input) {
  int n_subpops = input->n_subpops;
  
  for(int i=0; i < input->n_samples; i++) {
    if (input->samples[i].apriori_subpop == -1) {
      for(int h=0; h < 4; h++) {
	for(int j=0; j < input->n_windows; j++) {
	  fprintf(stderr,"sample %d/%d\twindow %d", i, h, j);
	  for(int k=0; k < n_subpops; k++) {
	    fprintf(stderr,"\t%1.3f", DF16(input->samples[i].est_p[h][ IDX(j,k) ]));
	  }
	  fprintf(stderr,"\n");
	}
      }
    }
  }
}

void random_forest(input_t *input) {
  /* These are file scope variables used to ensure we do not underflow or overflow int8_t 
     floating point encoding when setting values into input->samples[].est_p[][ IDX(,) ]. These
     are the end result of this entire file's computations */
  
  thread_args_t *args;
  MA(args, sizeof(thread_args_t), thread_args_t);
  args->input = input;
  args->next_window = 0;
  args->windows_complete = 0;
  args->rng = new md5rng(rfmix_opts.random_seed);
  
  pthread_mutex_init(&args->lock, NULL);
  
  pthread_t *threads;
  MA(threads, sizeof(pthread_t)*rfmix_opts.n_threads, pthread_t);

  for(int i=0; i < rfmix_opts.n_threads; i++)
    pthread_create(threads + i, NULL, random_forest_thread, (void *) args);

  for(int i=0; i < rfmix_opts.n_threads; i++)
    pthread_join(threads[i], NULL);
  fprintf(stderr,"\n");

#if 0
  dump_results(input);
#endif
  
  delete args->rng;
  free(threads);
  free(args);
}
