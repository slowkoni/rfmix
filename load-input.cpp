#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "kmacros.h"
#include "rfmix.h"
#include "genetic-map.h"
#include "load-input.h"
#include "inputline.h"

extern rfmix_opts_t rfmix_opts;

#define SUBPOP_ALLOC_STEP (8)
#define SAMPLE_ALLOC_STEP (256)
static void load_samples(input_t *input) {
  sample_t *samples;
  int n_samples, i, ref_idx;
  char *sample_id, *reference_pop, *p;
  
  Inputline *f = new Inputline(rfmix_opts.class_fname);

  /* Define the 0th "reference" subpop to mean a non-reference sample to be analyzed */
  RA(input->reference_subpops, sizeof(char *)*(input->n_subpops + SUBPOP_ALLOC_STEP), char *);
  input->reference_subpops[0] = strdup("unknown/admixed");
  input->n_subpops++;

  /* Now scan the sample map file and determine the reference subpops and sample mapping to them */
  while((p = f->nextline(INPUTLINE_NOCOPY)) != NULL) {
    CHOMP(p);
    if (p[0] == 0 || p[0] == '#') continue;

    sample_id = strsep(&p, "\t");
    reference_pop = strsep(&p, "\t");

    if (reference_pop[0] != 0 && reference_pop[0] != ' ' &&
	strcasecmp(reference_pop, "admixed")  != 0 &&
	strcasecmp(reference_pop, "query")    != 0 &&
	strcasecmp(reference_pop, "unknown")  != 0) {

      /* Search for this reference subpop in the already known list */
      for(i=0; i < input->n_subpops; i++)
	if (strcasecmp(input->reference_subpops[i], reference_pop) == 0) break;
      
      /* Add the subpop name to the list of reference subpops if it is not found */
      if (i == input->n_subpops) {
	if (input->n_subpops % SUBPOP_ALLOC_STEP == 0)
	  RA(input->reference_subpops, sizeof(char *)*(input->n_subpops + SUBPOP_ALLOC_STEP), char *);
	input->reference_subpops[input->n_subpops] = strdup(reference_pop);
	input->n_subpops++;
      }
      /* Whether just added or found, the reference subpop index is i */
      ref_idx = i;
    } else {
      /* empty fields, tokens starting with whitespace, and the 3 keywords above explicitly
	 indicate query/admixed/unknown samples. */
      ref_idx = 0;
    }

    if (n_samples % SAMPLE_ALLOC_STEP == 0)
      RA(samples, sizeof(sample_t)*(SAMPLE_ALLOC_STEP + n_samples), sample_t);
	
    samples[n_samples].sample_id = strdup(sample_id);
    samples[n_samples].apriori_subpop = ref_idx;
    samples[n_samples].haplotype[0] = NULL;
    samples[n_samples].haplotype[1] = NULL;
    samples[n_samples].current_p[0] = NULL;
    samples[n_samples].current_p[1] = NULL;
    samples[n_samples].est_p = NULL;
    n_samples++;
  }
  input->samples = samples;
  input->n_samples = n_samples;
  
  delete f;
}


input_t *load_input(void) {
  input_t *input;
  MA(input, sizeof(input_t), input_t);

  input->genetic_map = new GeneticMap();
  input->genetic_map->load_map(rfmix_opts.genetic_fname, rfmix_opts.chromosome);

  return input;
}


void free_input(input_t *input) {

  delete input->genetic_map;
  free(input);
}
