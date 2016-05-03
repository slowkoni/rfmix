/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

#ifndef RFMIX_H
#define RFMIX_H

/* From GNU autoconf/automake */
#include "config.h"

#include <vector>
#include "genetic-map.h"
#include "hash-table.h"

#include <math.h>
/* As a critical feature to trim memory usage, we are going to use a log odds
   encoding of a floating point number restricted to range 0.0 - 1.0 by 
   8-bit integers in range -127 to 127. This effectively creates 8 bit floats.
   The error for representation is at maximum about 2.0% and peaks in the
   center of the 0.0 to 1.0 range. The minimum possible and maximum possible
   floating point numbers are 0.000025334 and 0.9999747 respectively. The
   macro DF8(x) decodes the 8-bit integer to a double, and the inline function
   ef8 encodes a number in range min to max above to an int8_t. */
#define DF8(x) ( 1.0/(1.0+exp(((double) (x))/-25.0)) )
#define EF8(p) ( (int) ( -25.0*log( (1.0-(p))/(p) ) ) )

static inline int8_t ef8(double p) {
  int tmp = (int) ( -25.0*log( (1.0-p)/p ) );
  if (tmp < -127) tmp = -127;
  if (tmp >  127) tmp =  127;
  return tmp;
}

/* Likewise with int8_t float encodings above, but with int16_t giving much
   higher precision and range closer to 0.0 and 1.0. The maximum error is
   0.024% peaking at 0.5 and the range is 1.26765e-14 to what rounds to 1.0
   in R. This level of precision is needed for current_p[2][] because these
   arrays store the results of the forward-backward algorithm both for 
   feedback to EM and for output of results. The forward-backward calculations
   can produce a much finer degree of precision than random forest does with
   a typical 100 trees or so */
#define DF16(x) ( 1.0/( 1.0 + exp((double) (x)/-1024.0) ))
#define EF16(p) ( (int) ( -1024*log( (1.0-(p))/(p) ) ) )

static inline int16_t ef16(double p) {
  if (p <= 0.) return -32767;
  if (p >= 1.) return +32767;
  
  int tmp = (int) ( -1024.0 * log( (1.0 - p)/p ) );
  if (tmp < -32767) tmp = -32767;
  if (tmp >  32767) tmp =  32767;
  return tmp;
}

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
  double n_generations;
  double rf_window_size;
  double crf_spacing;
  int n_trees;
  int node_size;
  int reanalyze_reference;
  int em_iterations;
  int bootstrap_mode;
  int minimum_snps;
  int analyze_range[2];
  char *analyze_str;
  double crf_weight;

  int debug;
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
  char *ref;
  char *alt;
  char *snp_id;
} snp_t;

/* IMPORTANT: Because pointers take 8 bytes and the number of subpops is often small but 
              number of windows large, the array of pointers for windows to arrays of 
	      values for subpops per each window takes up as much or more memory than
	      the values themselves. So the marco IDX(w,s) performs the translation from
	      two dimensional indexing to one dimensional so we do not need the array 
	      of pointers. The variable n_subpops must be defined and set appropriately
	      in the function that uses IDX().

	      An alternative is to have the outer dimension be subpops and the inner one
	      be windows, but since most or all loops will access all subpops in a row
	      for each window, looping subpops within the window loop, doing it this way
	      produces better L2 cache performance by keeping the subpop values together */
#define IDX(w,s) ( (w)*n_subpops + (s) )
/* IMPORTANT: the current_p and est_p arrays are using the 8 bit float scheme discussed above */
typedef struct {
  char *sample_id;
  int apriori_subpop; // 0 means query/admixed/unknown sample. 1 through K, reference sample
  int8_t *haplotype[2];
  int8_t *msp[4];
  int8_t *ksp[2]; /* known state path, allocated and set only for internal simulated samples */
  double logl[4];
  int16_t *current_p[2]; // current estimate of probability of subpop [hap][ IDX(crf_window,subpop) ]
  int16_t *est_p[4]; // new estimate of probability of subpop estimate [hap][ IDX(crf_window,subpop) ]
  float *sis_p[2]; // Suyash stay-in-state forward-backward probability [hap][ crf_window ]

  int column_idx;
  int sample_idx;
  int s_parent;
  int s_sample;
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
enum { RF_BOOTSTRAP_FLAT=0, RF_BOOTSTRAP_HIERARCHICAL, RF_BOOTSTRAP_STRATIFIED, N_RF_BOOTSTRAP };

#define MINIMUM_GENETIC_DISTANCE (0.00001)
#define P_MINIMUM_FOR_REF (0.0)
#define RF_THREAD_WINDOW_CHUNK_SIZE (3)
#define CRF_SAMPLES_PER_BLOCK (32)
#define SIM_PARENT_PROPORTION (0.10)
#define SIM_GROWTH_RATE (1.20)
#define SIM_SAMPLES_PER_SUBPOP (200)

double crf(input_t *input, double w);

void msp_output(input_t *input);
void fb_output(input_t *input);
void fb_stay_in_state_output(input_t *input);
void output_Q(input_t *input);
  
#endif
