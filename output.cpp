#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>
#include <math.h>

#include "kmacros.h"
#include "rfmix.h"

extern rfmix_opts_t rfmix_opts;

static void msp_output_leader(FILE *f, snp_t *snps, int n_snps, crf_window_t *crfw, int n_windows,
			      int start, int end) {
  int start_snp, end_snp, n;

#if 1
  if (start == 0) {
    start_snp = 0;
  } else {
    start_snp = (crfw[start-1].snp_idx + crfw[start].snp_idx + 0.5)/2.;
  }
  if (end == n_windows) {
    end_snp = n_snps - 1;
    n = end_snp - start_snp + 1;
  } else {
    end_snp = (crfw[end].snp_idx + crfw[end-1].snp_idx)/2.;
    n = end_snp - start_snp;
  }
#else
  start_snp = start == 0 ? crfw[start].rf_start_idx : crfw[start].snp_idx;
  if (end == n_windows) {
    end--;
    end_snp = crfw[end].rf_end_idx;
    n = end_snp - start_snp + 1;
  } else {
    end_snp = crfw[end].snp_idx;
    n = end_snp - start_snp;
  }
#endif
  fprintf(f,"%s\t%d\t%d", rfmix_opts.chromosome, snps[start_snp].pos, snps[end_snp].pos);
  fprintf(f,"\t%1.2f\t%1.2f\t%d", snps[start_snp].genetic_pos, snps[end_snp].genetic_pos,
	    n);
}

static void msp_output_result(FILE *f, sample_t *samples, int n_samples, int w) {
  for(int j=0; j < n_samples; j++) {
    if (samples[j].apriori_subpop != -1) continue;
    fprintf(f,"\t%d\t%d", samples[j].msp[0][w], samples[j].msp[1][w]);
  }
}

static int msp_compare(sample_t *samples, int n_samples, int n_windows, int a, int b) {
  for(int j=0; j < n_samples; j++) {
    if (samples[j].apriori_subpop != -1) continue;
    if (samples[j].msp[0][a] != samples[j].msp[0][b] ||
	samples[j].msp[1][a] != samples[j].msp[1][b]) return 1;
  }

  return 0;
}

#define MSP_EXTENSION ".msp.tsv"
void msp_output(input_t *input) {
  int fname_length = strlen(rfmix_opts.output_basename) + strlen(MSP_EXTENSION) + 1;
  char fname[fname_length];

  sprintf(fname,"%s%s", rfmix_opts.output_basename, MSP_EXTENSION);
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
    msp_output_leader(f, input->snps, input->n_snps, input->crf_windows, input->n_windows, i, j);
    msp_output_result(f, input->samples, input->n_samples, i);
    fprintf(f,"\n");
    i = j;
  }
 
  fclose(f);
}

static void fb_output_haplotype(FILE *f, int16_t *p, int n) {
  fprintf(f,"%1.5f",DF16(p[0]));
  for(int k=1; k < n; k++)
    fprintf(f," %1.5f",DF16(p[k]));
}

#define FB_EXTENSION ".fb.tsv"
void fb_output(input_t *input) {
  fprintf(stderr,"Outputting forward-backward results.... \n");
  int fname_length = strlen(rfmix_opts.output_basename) + strlen(FB_EXTENSION) + 1;
  char fname[fname_length];

  sprintf(fname,"%s%s", rfmix_opts.output_basename, FB_EXTENSION);
  FILE *f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open output file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  }
  
  fprintf(f,"#chm\tpos\tgpos");
  for(int j=0; j < input->n_samples; j++) {
    sample_t *sample = input->samples + j;
    if (sample->apriori_subpop != -1) continue;

    fprintf(f,"\t%s.0\t%s.1", sample->sample_id, sample->sample_id);
  }
  fprintf(f,"\n");

  for(int i=0; i < input->n_windows; i++) {
    fprintf(f,"%s\t%d\t%1.2f", rfmix_opts.chromosome, input->snps[input->crf_windows[i].snp_idx].pos,
	    input->crf_windows[i].genetic_pos*100.);
    for(int j=0; j < input->n_samples; j++) {
      if (input->samples[j].apriori_subpop != -1) continue;
      for(int h=0; h < 2; h++) {
	fprintf(f,"\t");
	fb_output_haplotype(f, input->samples[j].current_p[h] + i*input->n_subpops, input->n_subpops);
      }
    }
    fprintf(f,"\n");
  }
 
  fclose(f);
}

#define SIS_EXTENSION ".sis.tsv"
void fb_stay_in_state_output(input_t *input) {
  fprintf(stderr,"Outputting Suyash-mode stay-in-state forward-backward results...  \n");
  int fname_length = strlen(rfmix_opts.output_basename) + strlen(SIS_EXTENSION) + 1;
  char fname[fname_length];

  sprintf(fname,"%s%s", rfmix_opts.output_basename, SIS_EXTENSION);
  FILE *f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open output file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  }
  fprintf(f,"#chm\tpos\tgpos");
  for(int j=0; j < input->n_samples; j++)
    fprintf(f,"\t%s.0\t%s.1", input->samples[j].sample_id, input->samples[j].sample_id);
  fprintf(f,"\n");
  
  for(int i=0; i < input->n_windows - 1; i++) {
    fprintf(f,"%s\t%d\t%1.2f", rfmix_opts.chromosome, input->snps[input->crf_windows[i].snp_idx].pos,
	    input->crf_windows[i].genetic_pos*100.);    
    for(int j=0; j < input->n_samples; j++) {
      if (input->samples[j].apriori_subpop != -1) continue;
      fprintf(f,"\t%1.5f\t%1.5f", input->samples[j].sis_p[0][i], input->samples[j].sis_p[1][i]);
    }
    fprintf(f,"\n");
  }
}
