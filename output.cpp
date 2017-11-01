/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

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
    if (samples[j].apriori_subpop != -1 || samples[j].s_sample == 1) continue;
    fprintf(f,"\t%d\t%d", samples[j].msp[0][w], samples[j].msp[1][w]);
  }
}

static int msp_compare(sample_t *samples, int n_samples, int n_windows, int a, int b) {
  for(int j=0; j < n_samples; j++) {
    if (samples[j].apriori_subpop != -1 || samples[j].s_sample == 1) continue;
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
  fprintf(f,"#");
  fprintf(f,"Subpopulation order/codes: %s=0", input->reference_subpops[0]);
  for(int i=1; i < input->n_subpops; i++) {
    fprintf(f,"\t%s=%d", input->reference_subpops[i], i);
  }
  fprintf(f,"\n");
  fprintf(f,"#chm\tspos\tepos\tsgpos\tegpos\tn snps");
  for(int j=0; j < input->n_samples; j++) {
    sample_t *sample = input->samples + j;
    if (sample->apriori_subpop != -1 || sample->s_sample == 1) continue;

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
    fprintf(f,"\t%1.5f",DF16(p[k]));
}

#define FB_EXTENSION ".fb.tsv"
void fb_output(input_t *input) {
/* 
the output is a tab separated file with the name \<output basename\>.fb.tsv

The first line is a comment line, that specifies the order of subpopulations:
eg:
  #reference_panel_population: golden_retriever  labrador_retriever  poodle  poodle_small
The second line specifies the column names, and every following lines gives 
data on a chunk of the genome, called a conditional random field (CRF) point.
The first few columns specifies the chromosome, genetic marker's physical 
  position in basepair units and genetic position in centiMorgans, and the 
  genetic marker's numerical index in the rfmix genetic map input file. 
  The remaining columns give the probabilities that the CRF point for a 
  genotype's haplotype was assigned to a specific reference panel population. 
  A genotype has two haplotypes, so the number of probabilities for a genotype 
  is 2*(number of reference panel populations). 
  The number of columns in the file is 
  4 + (number of genotypes) * 2 * (number of reference panel populations).

For example, for a rfmix run with 2 admixed genotype_ids run against 3 
reference panel populations, the columns would be:
  chromosome physical_position genetic_position genetic_marker_index 
  genotype_id1:::hap1:::subpop1
  genotype_id1:::hap1:::subpop2
  genotype_id1:::hap1:::subpop3
  genotype_id1:::hap2:::subpop1
  genotype_id1:::hap2:::subpop2
  genotype_id1:::hap2:::subpop3
  genotype_id2:::hap1:::subpop1
  genotype_id2:::hap1:::subpop2
  genotype_id2:::hap1:::subpop3
  genotype_id2:::hap2:::subpop1
  genotype_id2:::hap2:::subpop2
  genotype_id2:::hap2:::subpop3
  */

  fprintf(stderr,"Outputing forward-backward results.... \n");
  int fname_length = strlen(rfmix_opts.output_basename) + strlen(FB_EXTENSION) + 1;
  char fname[fname_length];

  sprintf(fname,"%s%s", rfmix_opts.output_basename, FB_EXTENSION);
  FILE *f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open output file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  }
  
  fprintf(f,"#");
  fprintf(f,"reference_panel_population:\t%s", input->reference_subpops[0]);
  for(int i=1; i < input->n_subpops; i++) {
    fprintf(f,"\t%s", input->reference_subpops[i]);
  }
  fprintf(f,"\n");
  fprintf(f,"chromosome\tphysical_position\tgenetic_position\tgenetic_marker_index");
  for(int j=0; j < input->n_samples; j++) {
    sample_t *sample = input->samples + j;
    if (sample->apriori_subpop != -1 || sample->s_sample == 1) continue;

    for(int k=0; k < input->n_subpops; k++) {
      fprintf(f, "\t%s:::hap1:::%s", sample->sample_id, input->reference_subpops[k]);
    }
    for(int k=0; k < input->n_subpops; k++) {
      fprintf(f, "\t%s:::hap2:::%s", sample->sample_id, input->reference_subpops[k]);
    }

  }
  fprintf(f,"\n");

  for(int i=0; i < input->n_windows; i++) {
    fprintf(f,"%s\t%d\t%1.5f\t%d", rfmix_opts.chromosome, input->snps[input->crf_windows[i].snp_idx].pos,
    	    input->crf_windows[i].genetic_pos*100.,input->crf_windows[i].snp_idx);
    for(int j=0; j < input->n_samples; j++) {
      if (input->samples[j].apriori_subpop != -1 || input->samples[j].s_sample == 1) continue;
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
  int fname_length = strlen(rfmix_opts.output_basename) + strlen(SIS_EXTENSION) + 1;
  char fname[fname_length];

  sprintf(fname,"%s%s", rfmix_opts.output_basename, SIS_EXTENSION);
  FILE *f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open output file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  }
  fprintf(f,"#chm\tpos\tgpos\tsnp idx");
  for(int j=0; j < input->n_samples; j++)
    fprintf(f,"\t%s.0\t%s.1", input->samples[j].sample_id, input->samples[j].sample_id);
  fprintf(f,"\n");
  
  for(int i=0; i < input->n_windows - 1; i++) {
    fprintf(f,"%s\t%d\t%1.5f\t%d", rfmix_opts.chromosome, input->snps[input->crf_windows[i].snp_idx].pos,
	    input->crf_windows[i].genetic_pos*100., input->crf_windows[i].snp_idx);    
    for(int j=0; j < input->n_samples; j++) {
      if (input->samples[j].apriori_subpop != -1 || input->samples[j].s_sample == 1) continue;
      fprintf(f,"\t%1.5f\t%1.5f", input->samples[j].sis_p[0][i], input->samples[j].sis_p[1][i]);
    }
    fprintf(f,"\n");
  }
}

#define Q_EXTENSION (".rfmix.Q")
void output_Q(input_t *input) {
  fprintf(stderr,"Outputing diploid global ancestry estimates.... \n");
  int fname_length = strlen(rfmix_opts.output_basename) + strlen(Q_EXTENSION) + 1;
  char fname[fname_length];

  sprintf(fname,"%s%s", rfmix_opts.output_basename, Q_EXTENSION);
  FILE *f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open output file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  }

  fprintf(f,"#rfmix diploid global ancestry .Q format output\n");
  fprintf(f,"#sample");
  for(int i=0; i < input->n_subpops; i++) {
    fprintf(f,"\t%s", input->reference_subpops[i]);
  }
  fprintf(f,"\n");

  double q[input->n_subpops];
  for(int i=0; i < input->n_samples; i++) {
    sample_t *sample = input->samples + i;
    if ((sample->apriori_subpop != -1 && rfmix_opts.em_iterations == 0) || sample->s_sample == 1) continue;
    
    for(int k=0; k < input->n_subpops; k++) q[k] = 0.;
    for(int j=0; j < input->n_windows; j++) {
      q[sample->msp[0][j]]++;
      q[sample->msp[1][j]]++;
    }
    for(int k=0; k < input->n_subpops; k++)
      q[k] = q[k]/(input->n_windows*2);
    fprintf(f,"%s", sample->sample_id);
    for(int k=0; k < input->n_subpops; k++)
      fprintf(f,"\t%1.5f", q[k]);
    fprintf(f,"\n");
  }

  fclose(f);
}
