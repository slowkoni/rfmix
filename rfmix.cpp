#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>
#include <time.h>

#include "cmdline-utils.h"
#include "kmacros.h"

#include "rfmix.h"
#include "load-input.h"
#include "random-forest.h"

rfmix_opts_t rfmix_opts;

static option_t options[] = {
  /* Input and output specification options (all are required) */
  { 'f', "query-file", &rfmix_opts.qvcf_fname, OPT_STR, 1, 1,
    "VCF file with samples to analyze (required)" },
  { 'r', "reference-file", &rfmix_opts.rvcf_fname, OPT_STR, 1, 1,
    "VCF file with reference individuals (required)" },
  { 'm', "sample-map", &rfmix_opts.class_fname, OPT_STR, 1, 1,
    "Reference panel sample population classification map (required)" },
  { 'g', "genetic-map", &rfmix_opts.genetic_fname, OPT_STR, 1, 1,
    "Genetic map file (required)" },
  { 'o', "output-basename", &rfmix_opts.output_basename, OPT_STR, 1, 1,
    "Basename (prefix) for output files (required)" },

  /* Tunable algorithm parameters (none are required - defaults are reasonable)*/
  {   0, "max-missing", &rfmix_opts.maximum_missing_data_freq, OPT_DBL, 0, 1,
      "Maximum proportion of missing data allowed to include a SNP" },
  { 'w', "rf-window-size", &rfmix_opts.rf_window_size, OPT_DBL, 0, 1,
    "Random forest window size (class estimation window size)" },
  { 'c', "crf-spacing", &rfmix_opts.crf_spacing, OPT_DBL, 0, 1,
    "Conditional Random Field spacing (# of SNPs)" },
  { 'G', "generations", &rfmix_opts.n_generations, OPT_DBL, 0, 1,
    "Average number of generations since expected admixture" },
  { 't', "trees", &rfmix_opts.n_trees, OPT_INT, 0, 1,
    "Number of tree in random forest to estimate population class probability" },
  {   0, "reanalyze-reference", &rfmix_opts.reanalyze_reference, OPT_FLAG, 0, 0,
      "After first iteration, include reference panel in analysis and reclassify" },

  /* Runtime execution control options (only specifies how the program runs)*/
  { 0, "n-threads", &rfmix_opts.n_threads, OPT_INT, 0, 1,
    "Force number of simultaneous thread for parallel execution" },
  { 0, "chromosome", &rfmix_opts.chromosome, OPT_STR, 1, 1,
    "Execute only on specified chromosome (currently required)" },
  { 0, "random-seed", &rfmix_opts.random_seed_str, OPT_STR, 0, 1,
    "Seed value for random number generation - integer value (maybe specified in"
    "hexadecimal by preceeding with 0x), or the string \"clock\" to seed with "
    "the current system time." },
  { 0, NULL, NULL, 0, 0, 0, NULL }
};
 
static void init_options(void) {
  rfmix_opts.qvcf_fname = (char *) "";
  rfmix_opts.rvcf_fname = (char *) "";
  rfmix_opts.genetic_fname = (char *) "";
  rfmix_opts.class_fname = (char *) "";
  rfmix_opts.output_basename = (char *) "";

  rfmix_opts.maximum_missing_data_freq = 0.05;
  rfmix_opts.rf_window_size = 0.2;
  rfmix_opts.crf_spacing = 0.1;
  rfmix_opts.n_generations = 8;
  rfmix_opts.n_trees = 100;

  rfmix_opts.n_threads = sysconf(_SC_NPROCESSORS_CONF);
  rfmix_opts.chromosome = (char *) "";
  rfmix_opts.random_seed_str = (char *) "0xDEADBEEF";
}

static void verify_options(void) {
  int stop = 0;
  
  if (rfmix_opts.qvcf_fname == NULL) {
    fprintf(stderr,"\nSpecify query/admixed VCF input file with -f option");
    stop = 1;
  }
  if (rfmix_opts.rvcf_fname == NULL) {
    fprintf(stderr,"\nSpecify reference VCF input file with -r option");
    stop = 1;
  }
  if (rfmix_opts.genetic_fname == NULL) {
    fprintf(stderr,"\nSpecify genetic map file with -g option");
    stop = 1;
  }
  if (rfmix_opts.class_fname == NULL) {
    fprintf(stderr,"\nSpecify reference sample subpopulation mapping with -m option");
    stop = 1;
  }
  if (rfmix_opts.output_basename == NULL) {
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
  if (rfmix_opts.n_generations < 2.) {
    fprintf(stderr,"\nNumber of generations since putative admixture must be 2 or larger");
    stop = 1;
  }
  if (rfmix_opts.n_trees < 10) {
    fprintf(stderr,"\nNumber of random forest trees must be at least 10");
    stop = 1;
  }
  if (rfmix_opts.n_threads < 1) rfmix_opts.n_threads = 1;
  if (rfmix_opts.chromosome == NULL) {
    fprintf(stderr,"\nSpecify VCF chromosome to analyze with -c option");
    stop = 1;
  }

  if (strcmp(rfmix_opts.random_seed_str, "clock") == 0) {
    rfmix_opts.random_seed = time(NULL);
  } else {
    rfmix_opts.random_seed = strtod(rfmix_opts.random_seed_str,0);
  }

  if (stop != 0) {
    fprintf(stderr,"\n\nCorrect command line errors to run rfmix. Run program with no options for help\n");
    exit(-1);
  }
}

int main(int argc, char *argv[]) {

  init_options();
  cmdline_getoptions(options, argc, argv);
  verify_options();

  input_t *rfmix_input = load_input();

  random_forest(rfmix_input);
  crf(rfmix_input);
  msp_output(rfmix_input);
  fb_output(rfmix_input);
  
  free_input(rfmix_input);
  return 0;
}
