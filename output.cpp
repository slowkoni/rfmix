#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>
#include <math.h>

#include "kmacros.h"
#include "rfmix.h"

extern rfmix_opts_t rfmix_opts;

static void msp_output_leader(FILE *f, snp_t *snps, crf_window_t *crfw, int n_windows,
			      int start, int end) {
  int start_snp, end_snp, n_snps;
  
  start_snp = start == 0 ? crfw[start].rf_start_idx : crfw[start].snp_idx;
  if (end == n_windows) {
    end--;
    end_snp = crfw[end].rf_end_idx;
    n_snps = end_snp - start_snp + 1;
  } else {
    end_snp = crfw[end].snp_idx;
    n_snps = end_snp - start_snp;
  }
    
  fprintf(f,"%s\t%d\t%d", rfmix_opts.chromosome, snps[start_snp].pos, snps[end_snp].pos);
  fprintf(f,"\t%1.2f\t%1.2f\t%d", snps[start_snp].genetic_pos, snps[end_snp].genetic_pos,
	    n_snps);
}

static void msp_output_result(FILE *f, sample_t *samples, int n_samples, int w) {
  for(int j=0; j < n_samples; j++) {
    if (samples[j].apriori_subpop != -1) continue;
    fprintf(f,"\t%d\t%d", samples[j].msp[0][w], samples[j].msp[1][w]);
  }
}

static int msp_compare(sample_t *samples, int n_samples, int n_windows, int a, int b) {
  for(int j=0; j < n_windows; j++) {
    if (samples[j].apriori_subpop != -1) continue;
    if (samples[j].msp[0][a] != samples[j].msp[0][b] ||
	samples[j].msp[1][a] != samples[j].msp[1][b]) return 1;
  }

  return 0;
}

void msp_output(input_t *input) {
  int fname_length = strlen(rfmix_opts.output_basename) + strlen(".msp.tsv") + 1;
  char fname[fname_length];

  sprintf(fname,"%s.msp.tsv", rfmix_opts.output_basename);
  FILE *f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open output file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  }
  
  fprintf(f,"#chm\tspos\tepos\tsgpos\tegpos\tn snps");
  for(int j=0; j < input->n_samples; j++) {
    sample_t *sample = input->samples + j;
    if (sample->apriori_subpop != -1) continue;

    fprintf(f,"\t%s.0\t%s.1", sample->sample_id, sample->sample_id);
  }
  fprintf(f,"\n");


  int i = 0;
  while(i < input->n_windows) {
    int j = i + 1;
    while(j < input->n_windows &&
	  msp_compare(input->samples, input->n_samples,	input->n_windows, i, j) == 0) j++;
    msp_output_leader(f, input->snps, input->crf_windows, input->n_windows, i, j);
    msp_output_result(f, input->samples, input->n_samples, i);
    fprintf(f,"\n");
    i = j;
  }
 
  fclose(f);
}



void output(input_t *input) {

  msp_output(input);
  //  fb_output(input);
}
