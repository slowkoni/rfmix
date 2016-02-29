#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <time.h>

/* Local includes */
#include "cmdline-utils.h"
#include "kmacros.h"

#include "rfmix.h"
#include "load-input.h"
#include "random-forest.h"

rfmix_opts_t rfmix_opts;
int em_iteration;

static option_t options[] = {
  /* Input and output specification options (all are required) */
  { 'f', "query-file", &rfmix_opts.qvcf_fname, OPT_STR, 1, 1,
    "VCF file with samples to analyze                      (required)" },
  { 'r', "reference-file", &rfmix_opts.rvcf_fname, OPT_STR, 1, 1,
    "VCF file with reference individuals                   (required)" },
  { 'm', "sample-map", &rfmix_opts.class_fname, OPT_STR, 1, 1,
    "Reference panel sample population classification map  (required)" },
  { 'g', "genetic-map", &rfmix_opts.genetic_fname, OPT_STR, 1, 1,
    "Genetic map file                                      (required)" },
  { 'o', "output-basename", &rfmix_opts.output_basename, OPT_STR, 1, 1,
    "Basename (prefix) for output files                    (required)" },
  { 0, "chromosome", &rfmix_opts.chromosome, OPT_STR, 1, 1,
    "Execute only on specified chromosome                  (required)\n" },

  /* Tunable algorithm parameters (none are required - defaults are reasonable)*/
  { 'c', "crf-spacing", &rfmix_opts.crf_spacing, OPT_DBL, 0, 1,
    "Conditional Random Field spacing (# of SNPs)" },
  { 's', "rf-window-size", &rfmix_opts.rf_window_size, OPT_DBL, 0, 1,
    "Random forest window size (class estimation window size)" },
  { 'w', "crf-weight", &rfmix_opts.crf_weight, OPT_DBL, 0, 1,
    "Weight of observation term relative to transition term in conditional random field" },
  { 'G', "generations", &rfmix_opts.n_generations, OPT_DBL, 0, 1,
    "Average number of generations since expected admixture" },
  { 'e', "em-iterations", &rfmix_opts.em_iterations, OPT_INT, 0, 1,
    "Maximum number of EM iterations" },
  {  0, "reanalyze-reference", &rfmix_opts.reanalyze_reference, OPT_FLAG, 0, 0,
     "In EM, analyze local ancestry of the reference panel and reclassify it\n" },

  { 'n', "node-size", &rfmix_opts.node_size, OPT_INT, 0, 1,
    "Terminal node size for random forest trees" },
  { 't', "trees", &rfmix_opts.n_trees, OPT_INT, 0, 1,
    "Number of tree in random forest to estimate population class probability" },
  {  0, "max-missing", &rfmix_opts.maximum_missing_data_freq, OPT_DBL, 0, 1,
      "Maximum proportion of missing data allowed to include a SNP" },
  { 'b', "bootstrap-mode", &rfmix_opts.bootstrap_mode, OPT_INT, 0, 1,
    "Specify random forest bootstrap mode as integer code (see manual)" },
  { 0, "rf-minimum-snps", &rfmix_opts.minimum_snps, OPT_INT, 0, 1,
    "With genetic sized rf windows, include at least this many SNPs regardless of span" },
  { 0, "analyze-range", &rfmix_opts.analyze_str, OPT_STR, 0, 1,
    "Physical position range, specified as <start pos>-<end pos>, in Mbp (decimal allowed)\n" },
  
  /* Runtime execution control options (only specifies how the program runs)*/
  { 0, "n-threads", &rfmix_opts.n_threads, OPT_INT, 0, 1,
    "Force number of simultaneous thread for parallel execution" },
  { 0, "random-seed", &rfmix_opts.random_seed_str, OPT_STR, 0, 1,
    "Seed value for random number generation (integer)\n"
    "\t(maybe specified in hexadecimal by preceeding with 0x), or the string \"clock\" \n"
    "\tto seed with the current system time." },
  { 0, NULL, NULL, 0, 0, 0, NULL }
};

static void init_options(void) {
  rfmix_opts.qvcf_fname = (char *) "";
  rfmix_opts.rvcf_fname = (char *) "";
  rfmix_opts.genetic_fname = (char *) "";
  rfmix_opts.class_fname = (char *) "";
  rfmix_opts.output_basename = (char *) "";

  rfmix_opts.maximum_missing_data_freq = 0.05;
  rfmix_opts.rf_window_size = 50;
  rfmix_opts.crf_spacing = 5;
  rfmix_opts.n_generations = 8;
  rfmix_opts.n_trees = 100;
  rfmix_opts.node_size = 2;
  rfmix_opts.bootstrap_mode = 1;
  rfmix_opts.em_iterations = 0;
  rfmix_opts.minimum_snps = 10;
  rfmix_opts.analyze_str = (char *) "";
  rfmix_opts.analyze_range[0] = INT_MIN;
  rfmix_opts.analyze_range[1] = INT_MAX;
  rfmix_opts.crf_weight = 25.0;
  rfmix_opts.reanalyze_reference = 0;
  
  rfmix_opts.n_threads = sysconf(_SC_NPROCESSORS_CONF);
  rfmix_opts.chromosome = (char *) "";
  rfmix_opts.random_seed_str = (char *) "0xDEADBEEF";
}

static void print_banner(void) {
  fprintf(stderr,
"\n"
"RFMIX %s - Local Ancestry and Admixture Inference\n"
"(c) 2016 Mark Koni Hamilton Wright\n"
"Bustamante Lab - Stanford University School of Medicine\n"
"Based on concepts developed in RFMIX v1 by Brian Keith Maples, et al.\n"
"\n"
"This version is licensed for non-commercial academic research use only\n"
"For commercial licensing, please contact cdbadmin@stanford.edu\n"
"\n"
"--- For use in scientific publications please cite original publication ---\n"
"Brian Maples, Simon Gravel, Eimear E. Kenny, and Carlos D. Bustamante (2013).\n"
"RFMix: A Discriminative Modeling Approach for Rapid and Robust Local-Ancestry\n"
"Inference. Am. J. Hum. Genet. 93, 278-288\n"
"\n", VERSION);
}

static void verify_options(void) {
  int stop = 0;
  
  if (strcmp(rfmix_opts.qvcf_fname,"") == 0) {
    fprintf(stderr,"\nSpecify query/admixed VCF input file with -f option");
    stop = 1;
  }
  if (strcmp(rfmix_opts.rvcf_fname,"") == 0) {
    fprintf(stderr,"\nSpecify reference VCF input file with -r option");
    stop = 1;
  }
  if (strcmp(rfmix_opts.genetic_fname,"") == 0) {
    fprintf(stderr,"\nSpecify genetic map file with -g option");
    stop = 1;
  }
  if (strcmp(rfmix_opts.class_fname,"") == 0) {
    fprintf(stderr,"\nSpecify reference sample subpopulation mapping with -m option");
    stop = 1;
  }
  if (strcmp(rfmix_opts.output_basename,"") == 0) {
    fprintf(stderr,"\nSpecify output files basename (prefix) with -o option");
    stop = 1;
  }

  if (rfmix_opts.maximum_missing_data_freq < 0.0 || rfmix_opts.maximum_missing_data_freq > 1.0) {
    fprintf(stderr,"\nRange for --max-missing option is 0.0 to 1.0");
    stop = 1;
  }
  if (rfmix_opts.rf_window_size <= 0.) {
    fprintf(stderr,"\nRandom Forest window size must be greater than 0");
    stop = 1;
  }
  if (rfmix_opts.crf_spacing <= 0) {
    fprintf(stderr,"\nConditional random field size must be larger than 0");
    stop = 1;
  }
  if (rfmix_opts.n_generations <= 0.) {
    // and it really only makes sense 2 or larger, but smaller values useful for testing
    // penalizing recombination
    fprintf(stderr,"\nNumber of generations since putative admixture must be larger than 0.");
    stop = 1;
  }
  if (rfmix_opts.n_trees < 10) {
    fprintf(stderr,"\nNumber of random forest trees must be at least 10");
    stop = 1;
  }
  if (rfmix_opts.node_size < 2) {
    fprintf(stderr,"\nRandom forest node size must be at least 2");
    stop = 1;
  }
  if (rfmix_opts.bootstrap_mode < 0 || rfmix_opts.bootstrap_mode >= N_RF_BOOTSTRAP) {
    fprintf(stderr,"\nBootstrap mode (-b) out of valid range - see manual");
    stop = 1;
  }
  if (strcmp(rfmix_opts.analyze_str, "") != 0) {
    char *p, *start, *end;
    end = p = strdup(rfmix_opts.analyze_str);
    start = strsep(&end,"-");
    if (start == NULL || end == NULL) {
      fprintf(stderr,"Invalid physical range to analyze (--analyze-range)\n");
      stop = 1;
    } else {
      rfmix_opts.analyze_range[0] = atof(start)*1e6;
      rfmix_opts.analyze_range[1] = atof(end)*1e6;
      fprintf(stderr,"NOTICE: Analysis restricted to positions in range %d to %d\n", rfmix_opts.analyze_range[0],
	      rfmix_opts.analyze_range[1]);
    }
    free(p);
  }
  
  if (rfmix_opts.n_threads < 1) rfmix_opts.n_threads = 1;
  if (strcmp(rfmix_opts.chromosome,"") == 0) {
    fprintf(stderr,"\nSpecify VCF chromosome to analyze with -c option");
    stop = 1;
  }
  
  if (strcmp(rfmix_opts.random_seed_str, "clock") == 0) {
    rfmix_opts.random_seed = time(NULL);
  } else {
    rfmix_opts.random_seed = strtod(rfmix_opts.random_seed_str,0);
  }
  /* set random seed of the system random number generator for any cases it is
     used, usually temporary hacks/tests/debugging. Otherwise md5rng is used 
     for repeatability even when multiple threads are used */
  srand(rfmix_opts.random_seed);
  
  if (stop != 0) {
    fprintf(stderr,"\n\nCorrect command line errors to run rfmix. Run program with no options for help\n");
    exit(-1);
  }
}

int main(int argc, char *argv[]) {

  print_banner();
  init_options();
  cmdline_getoptions(options, argc, argv);
  verify_options();

  fprintf(stderr,"\n");
  input_t *rfmix_input = load_input();

  fprintf(stderr,"\n");
  em_iteration = 0;
  random_forest(rfmix_input);
  double logl = crf(rfmix_input);
  msp_output(rfmix_input);
  fb_output(rfmix_input);
  fb_stay_in_state_output(rfmix_input);
  output_Q(rfmix_input);
  fprintf(stderr,"Initial analysis - logl %1.1f\n\n", logl);

  for(int i=0; i < rfmix_opts.em_iterations; i++) {
    em_iteration = i + 1;
    double last_logl = logl;
    random_forest(rfmix_input);
    logl = crf(rfmix_input);
    double d = logl - last_logl;

    msp_output(rfmix_input);
    fb_output(rfmix_input);
    fb_stay_in_state_output(rfmix_input);
    output_Q(rfmix_input);
    
    fprintf(stderr,"EM iteration %d - logl %1.1f (%+1.1f)\n\n", i+1, logl, d);
    if (i > 0 && d < 0.1) break;
  }
  
  
  free_input(rfmix_input);
  return 0;
}
