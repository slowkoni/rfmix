#ifndef RFMIX_H
#define RFMIX_H

#include "genetic-map.h"
#include "hash-table.h"

#include <math.h>
/* As a critical feature to trim memory usage, we are going to use a log odds
   encoding of a floating point number restricted to range 0.0 - 1.0 by 
   8-bit integers in range -127 to 127. This effectively creates 8 bit floats.
   The error for representation is at maximum about 0.5% and peaks in the
   center of the 0.0 to 1.0 range. The minimum possible and maximum possible
   floating point numbers are 0.005008156 and 0.994991844 respectively. The
   macro DF8(x) decodes the 8-bit integer to a double, and EF8(p) encodes
   a number in range min to max above to an integer. It is responsibility
   of the user of EF8(p) to set range to -127 to 127 and cast to int8_t */
#define DF8(x) ( 1.0/(1.0+exp(((double) (x))/-24.0)) )
#define EF8(p) ( (int) ( -24.0*log( (1.0-(p))/(p) ) ) )


/* Alternative with uniform rounding error over the range, and full range 0.0
   to 1.0. Also less expensive to encode, decode. The log-odds formulation
   above can be tweaked to give better accuracy at the tails at the expense
   of rounding error in the middle of the range, or vice-versa. */
//#define DF8(x) ((x)/255.0)
//#define EF8(p) ((p)*255.0)

/* Program command line and configuration options - see rfmix.c for option definitions
   and default values set in init_options(). The global variable rfmix_opts, declared 
   and set in rfmix.c is referenced all over the program for these values where needed */
typedef struct {
  char *qvcf_fname;
  char *rvcf_fname;
  char *genetic_fname;
  char *class_fname;
  char *output_basename;

  double maximum_missing_data_freq;
  int rf_window_size;
  int crf_spacing;
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
  double genetic_pos;
} crf_window_t;

typedef struct {
  int pos;
  AF_TYPE genetic_pos;
  int crf_index; // CRF window index for this snp
} snp_t;

/* IMPORTANT: the current_p and est_p arrays are using the 8 bit float scheme discussed above */
typedef struct {
  char *sample_id;
  int apriori_subpop; // 0 means query/admixed/unknown sample. 1 through K, reference sample
  int8_t *haplotype[2];
  int8_t *msp[4];
  int8_t **current_p[2]; // current estimate of probability of subpop [hap][crf_window][subpop]
  int8_t **est_p[4]; // new estimate of probability of subpop estimate [hap][crf_window][subpop]
} sample_t;

typedef struct {
  int n_subpops;
  char **reference_subpops; // string names of the reference subpops

  int n_samples;
  sample_t *samples;
  HashTable *sample_hash;
  
  int n_snps;
  snp_t *snps;
  
  int n_windows;
  crf_window_t *crf_windows;

  GeneticMap *genetic_map;
} input_t;

/* This can be anything. The value I put here I pulled out of my backside. */
#define RFOREST_RNG_KEY 0x949FC1AD

void crf(input_t *input);

#endif
