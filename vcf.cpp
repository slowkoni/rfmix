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

#include "kmacros.h"
#include "inputline.h"
#include "rfmix.h"
#include "vcf.h"

int VCF::n_alleles = 0;
char **VCF::alleles = NULL;

char *VCF::find_allele_string(char *q) {

  for(int i=0; i < VCF::n_alleles; i++)
    if (strcmp(q, VCF::alleles[i]) == 0) return VCF::alleles[i];
  
  if (VCF::n_alleles % ALLELE_ALLOC_STEP == 0)
    RA(VCF::alleles, sizeof(char *)*(n_alleles + ALLELE_ALLOC_STEP), char *);
  VCF::alleles[VCF::n_alleles] = strdup(q);
  char *rval = alleles[VCF::n_alleles];
  
  VCF::n_alleles++;
  return rval;
}

static char *vcf_skip_headers(Inputline *vcf) {
  char *p;

  while((p = vcf->nextline(INPUTLINE_NOCOPY)) != NULL) {
    CHOMP(p);
    if (strncmp(p, "#CHROM", 6) == 0) break;
    if (p[0] == 0 || p[0] == '#') continue;
  }
  if (p == NULL) {
    fprintf(stderr,"\nNo genotype data found in VCF file %s.\n\n", vcf->fname);
    exit(-1);
  }
  return p;
}

static int8_t get_allele(char q) {
  switch(q) {
  case '0':
    return 0;
  case '1':
    return 1;
  case '.':
  default:
    return 2;
  }
}

static void skip_to_chromosome(Inputline *vcf, char *chm) {
  char *p, *q;
  
  while((p = vcf->nextline(INPUTLINE_NOCOPY)) != NULL) {
    q = strsep(&p, "\t");
    if (strcmp(q, chm) == 0) {
      vcf->pushback();
      return;
    }
  }

  fprintf(stderr,"\nCan't find chromosome %s in VCF file %s\n\n", chm, vcf->fname);
  exit(-1);
}

void VCF::parse_samples(char *sample_line) {
  char *p, *q;
  
  p = sample_line;
  for(int i=0; i < VCF_LEAD_COLS; i++) strsep(&p, "\t");

  samples = NULL;
  n_samples = 0;
  int col_idx = VCF_LEAD_COLS;
  while((q = strsep(&p, "\t")) != NULL) {
    if (n_samples % SAMPLE_ALLOC_STEP == 0)
      RA(samples, sizeof(sample_t)*(n_samples + SAMPLE_ALLOC_STEP), sample_t);

    samples[n_samples].sample_idx = n_samples;
    samples[n_samples].column_idx = col_idx;
    samples[n_samples].sample_id = strdup(q);
    samples[n_samples].haplotype[0] = NULL;
    samples[n_samples].haplotype[1] = NULL;
    
    n_samples++;
    col_idx++;
  }

  n_columns = col_idx;
  MA(column_map, sizeof(sample_t *)*(VCF_LEAD_COLS + n_columns), sample_t *);
  
  HashTable *sample_map = new HashTable(256);
  for(int i=0; i < n_samples; i++) {
    sample_map->insert(samples[i].sample_id, samples + i);
    column_map[samples[i].column_idx] = samples + i;
  }
}

void VCF::load_snps(char *chromosome, GeneticMap *genetic_map) {
  char *p, *q;

  snps = NULL;
  n_snps = 0;
  
  /* Note, the VCF file may need to be rewound/reopened to the beginning depending
     on what was done before this is called [again] */
  skip_to_chromosome(f, chromosome);
  while((p = f->nextline(INPUTLINE_NOCOPY)) != NULL) {
    CHOMP(p);
    if (p[0] == 0 || p[0] == '#') continue;
    
    char *chm = strsep(&p, "\t");
    if (strcmp(chm, chromosome) != 0) break;

    if (n_snps % SNP_ALLOC_STEP == 0)
      RA(snps, sizeof(snp_t)*(n_snps + SNP_ALLOC_STEP), snp_t);

    q = strsep(&p, "\t");
    snps[n_snps].pos = atoi(q);
    snps[n_snps].genetic_pos = genetic_map->translate_seqpos(snps[n_snps].pos);
    snps[n_snps].snp_id = strdup(strsep(&p,"\t"));
    q = strsep(&p, "\t");
    snps[n_snps].ref = VCF::find_allele_string(q);
    q = strsep(&p, "\t");
    snps[n_snps].alt = VCF::find_allele_string(q);
    
    n_snps++;
  }

}

void VCF::load_haplotypes(char *chromosome) {
  char *p, *q;

  delete f;
  f = new Inputline(fname, chromosome);
  vcf_skip_headers(f);
  skip_to_chromosome(f, chromosome);

  for(int i=0; i < n_samples; i++) {
    MA(samples[i].haplotype[0], sizeof(int8_t)*n_snps, int8_t);
    MA(samples[i].haplotype[1], sizeof(int8_t)*n_snps, int8_t);
    for(int j=0; j < n_samples; j++) {
      samples[i].haplotype[0][j] = 2;
      samples[i].haplotype[1][j] = 2;
    }
  }
  
  int s = 0;
  while(s < n_snps && (p = f->nextline()) != NULL) {
    CHOMP(p);
    if (p[0] == 0 || p[0] == '#') continue;

    q = strsep(&p, "\t");
    if (strcmp(chromosome, q) != 0) break;

    q = strsep(&p, "\t");
    int pos = atoi(q);

    if (snps[s].pos < pos) {
      fprintf(stderr,"VCF Error: A SNP in the snp array is not present in the VCF file\n");
      exit(-1);
    }
    if (snps[s].pos > pos) continue;

    int col_idx;
    for(col_idx=2; col_idx < VCF_LEAD_COLS; col_idx++) strsep(&p, "\t");
    while((q = strsep(&p, "\t")) != NULL) {
      sample_t *sample = column_map[col_idx];
      if (sample == NULL) continue;

      sample->haplotype[0][s] = get_allele(q[0]);
      sample->haplotype[1][s] = get_allele(q[2]);
      col_idx++;
    }
    s++;
  }

}

VCF::VCF(char *vcf_fname, char *chm) {

  sample_map = NULL;
  samples = NULL;
  n_samples = 0;

  n_columns = 0;
  column_map = NULL;

  n_snps = 0;
  snps = NULL;
  fname = strdup(vcf_fname);

  this->f = new Inputline(fname, chm);
  char *sample_line = vcf_skip_headers(f);
  parse_samples(sample_line);
}

void VCF::DropHaplotypes(void) {
  for(int i=0; i < n_samples; i++) {
    if (samples[i].haplotype[0]) free(samples[i].haplotype[0]);
    if (samples[i].haplotype[1]) free(samples[i].haplotype[1]);
    samples[i].haplotype[0] = NULL;
    samples[i].haplotype[1] = NULL;
  }
}

VCF::~VCF(void) {
  if (sample_map) delete sample_map;
  if (samples) {
    for(int i=0; i < n_samples; i++) {
      free(samples[i].sample_id);
      if (samples[i].haplotype[0]) free(samples[i].haplotype[0]);
      if (samples[i].haplotype[1]) free(samples[i].haplotype[1]);
    }
    free(samples);
  }
  n_samples = 0;
  /*  for(int i=0; i < n_alleles; i++) {
    free(alleles[i]);
  }
  free(alleles);*/
  
  if (column_map) free(column_map);
  n_columns = 0;
  
  if (snps) free(snps);
  n_snps = 0;

  if (f) delete f;
  if (fname) free(fname);
}
