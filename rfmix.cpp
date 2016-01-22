#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>

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
  { 'w', "rf-window-size", &rfmix_opts.rf_window_size, OPT_INT, 0, 1,
    "Random forest window size (class estimation window size)" },
  { 'c', "crf-spacing", &rfmix_opts.crf_spacing, OPT_INT, 0, 1,
    "Conditional Random Field spacing (# of SNPs)" },
  { 'g', "generations", &rfmix_opts.generations, OPT_INT, 0, 1,
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
  rfmix_opts.rf_window_size = 20;
  rfmix_opts.crf_spacing = 5;
  rfmix_opts.generations = 8;
  rfmix_opts.n_trees = 100;

  rfmix_opts.n_threads = sysconf(_SC_NPROCESSORS_CONF);
  rfmix_opts.chromosome = (char *) "";
  rfmix_opts.random_seed_str = (char *) "0xDEADBEEF";
}

static void verify_options(void) {
}

int main(int argc, char *argv[]) {

  init_options();
  cmdline_getoptions(options, argc, argv);
  verify_options();

  input_t *rfmix_input = load_input();

  random_forest(rfmix_input);
  
  free_input(rfmix_input);
  return 0;
}
