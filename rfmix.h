#ifndef RFMIX_H
#define RFMIX_H

#include "genetic-map.h"

/* Program command line and configuration options - see rfmix.c for option definitions
   and default values set in init_options(). The global variable rfmix_opts, declared 
   and set in rfmix.c is referenced all over the program for these values where needed */
typedef struct {
  char *qvcf_fname;
  char *rvcf_fname;
  char *genetic_fname;
  char *class_fname;
  char *output_basename;

  int rf_window_size;
  int crf_window_size;
  int generations;
  int n_trees;
  int reanalyze_reference;
  
  int n_threads;
  char *chromosome;
  char *random_seed_str;
  int random_seed;  /* set by parsing random_seed_str which might be "clock" or a hex number */
} rfmix_opts_t;

/* I am using AF_TYPE to mean either float or double, depending on how set here, so
   that it is simple to change the program to operate either entirely in floats for 
   the large arrays of floating point numbers, or in doubles. float (32 bit) are 
   preferred for memory conservation, but doubles (64 bit) are usually faster in
   execution on modern systems */
#define AF_TYPE float

/* The chromosome is broken up into discrete segments, possibly one at each input SNP, on
   which the conditional random field is defined. For training the random forests estimating
   the probabilities at each CRF point/window, the SNPs used may come from a larger region 
   overlapping more than one CRF window. rf_start_idx and rf_end_idx indicate the first and
   last SNP (inclusive) to be included in training the random forest */
typedef struct {
  int snp_idx;
  int rf_start_idx;
  int rf_end_idx;
} crf_window_t;

typedef struct {
  int pos;
  AF_TYPE genetic_pos;
  int crf_index; // CRF window index for this snp
} snp_t;

typedef struct {
  char *sample_id;
  int apriori_subpop; // 0 means query/admixed/unknown sample. 1 through K, reference sample
  int8_t *haplotype[2];
  AF_TYPE **current_p[2]; // current estimate of probability of subpop [hap][crf_window][subpop]
  AF_TYPE **est_p; // new estimate of probability of subpop estimate [crf_window][crf state]
} sample_t;

typedef struct {
  int n_subpops;
  char **reference_subpops; // string names of the reference subpops

  int n_samples;
  sample_t *samples;
  
  int n_snps;
  snp_t *snps;
  
  int n_windows;
  crf_window_t *crf_windows;
  
  GeneticMap *genetic_map;
} input_t;

#endif
