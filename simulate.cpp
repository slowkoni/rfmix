#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <stdint.h>

#include "kmacros.h"
#include "cmdline-utils.h"
#include "vcf.h"
#include "s-sample.h"

typedef struct {
  char *vcf_fname;
  char *class_fname;
  char *genetic_fname;
  char *output_basename;
  char *chromosome;
  int n_generations;
  char *random_seed_str;
  int32_t random_seed;
} opts_t;

opts_t opts;

static option_t options[] = {
  { 'f', "vcf", &opts.vcf_fname, OPT_STR, 1, 1,
    "Name of input VCF file" },
  { 'm', "sample-map", &opts.class_fname, OPT_STR, 1, 1,
    "Sample subpop mapping file - also selects which samples will be used for simulation" },
  { 'g', "genetic-map", &opts.genetic_fname, OPT_STR, 1, 1,
    "Genetic map file (required)" },
  { 'o', "output-basename", &opts.output_basename, OPT_STR, 1, 1,
    "Basename (prefix) for output files (required)" },

  { 'c', "chromosome", &opts.chromosome, OPT_STR, 1, 1,
    "Chromosome to select from the VCF file" },
  { 'G', "generations", &opts.n_generations, OPT_INT, 0, 1,
    "Number of generations to simulate random mating admixture" },

  { 0, "random_seed", &opts.random_seed_str, OPT_STR, 0, 1,
    "Seed value for random number generation - integer value (maybe specified in"
    "hexadecimal by preceeding with 0x), or the string \"clock\" to seed with "
    "the current system time." },
  { 0, NULL, NULL, 0, 0, 0, NULL }
};

static void init_options(void) {
  opts.vcf_fname = NULL;
  opts.class_fname = NULL;
  opts.genetic_fname = NULL;
  opts.output_basename = NULL;

  opts.chromosome = NULL;
  opts.n_generations = 8;
  opts.random_seed_str = (char *) "0xDEADBEEF";
}

static void verify_options(void) {

  if (opts.vcf_fname == NULL) {
    fprintf(stderr,"\nSpecify VCF input file for source data with -f option\n\n");
    exit(-1);
  }
  if (opts.output_basename == NULL) {
    fprintf(stderr,"\nSpecify the output basename (prefix) with -o option\n\n");
    exit(-1);
  }
  if (opts.genetic_fname == NULL) {
    fprintf(stderr,"\nA genetic map is required, specify with -g option\n\n");
    exit(-1);
  }
  if (opts.class_fname == NULL) {
    fprintf(stderr,"\nSpecify sample to subpopulation mapping file with -m option\n\n");
    exit(-1);
  }

  if (opts.chromosome == NULL) {
    fprintf(stderr,"\nSpecify chromosome to select from VCF input with -c option\n\n");
    exit(-1);
  }

  if (strcmp(opts.random_seed_str, "clock") == 0) {
    opts.random_seed = time(NULL);
  } else {
    opts.random_seed = strtod(opts.random_seed_str,0);
  }
}


int main(int argc, char *argv[]) {

  init_options();
  cmdline_getoptions(options, argc, argv);
  verify_options();
  
  VCF *vcf = new VCF(opts.vcf_fname, opts.chromosome);

  vcf->load_snps(opts.chromosome);
  vcf->load_haplotypes(opts.chromosome);

  Sample **parents = new Sample*[vcf->n_samples];
  for(int i=0; i < vcf->n_samples; i++) {
    parents[i] = new Sample(vcf->samples[i].sample_id, 0, vcf->snps, vcf->n_snps,
			    vcf->samples[i].haplotypes[0], vcf->samples[i].haplotypes[1]);
  }

  for(int g=0; g < opts.n_generations; g++) {
    Sample **children = new Sample*[vcf->n_samples];
    for(int i=0; i < vcf->n_samples; i++) {
      int p1_idx = rand()/(RAND_MAX + 1.0) * vcf->n_samples;
      int p2_idx = rand()/(RAND_MAX + 1.0) * vcf->n_samples;
      
      children[i] = new Sample(parents[p1_idx], parents[p2_idx]);
    }

    for(int i=0; i < vcf->n_samples; i++)
      delete parents[i];
    delete[] parents;

    parents = children;
  }
  
  return 0;
}
  
