#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <stdint.h>
#include <math.h>

#include "kmacros.h"
#include "s-sample.h"

Sample::Sample(char *sample_id, int subpop, snp_t *snps, int n_snps, int8_t *h1, int8_t *h2) {

  this->sample_id = strdup(sample_id);
  this->snps = snps;
  this->n_snps = n_snps;
  
  haplotype[0] = new int8_t[n_snps];
  haplotype[1] = new int8_t[n_snps];
  this->subpop[0] = new int8_t[n_snps];
  this->subpop[1] = new int8_t[n_snps];
  
  for(int i=0; i < n_snps; i++) {
    haplotype[0][i] = h1[i];
    haplotype[1][i] = h2[i];
    this->subpop[0][i] = subpop;
    this->subpop[1][i] = subpop;
  } 
}

Sample::Sample(Sample *p1, Sample *p2) {
  sample_id = new char[12];
  snprintf(sample_id, 12, "%08x", rand());

  if (p1->n_snps != p2->n_snps ||
      p1->snps != p2->snps) {
    fprintf(stderr,"Error: Trying to breed a new sample from two with incompatible SNPs\n");
    exit(-1);
  }
  
  snps = p1->snps;
  n_snps = p1->n_snps;

  haplotype[0] = new int8_t[n_snps];
  haplotype[1] = new int8_t[n_snps];
  subpop[0] = new int8_t[n_snps];
  subpop[1] = new int8_t[n_snps];
  
  p1->meiosis(haplotype[0], subpop[0]);
  p2->meiosis(haplotype[1], subpop[1]);
}

Sample::~Sample(void) {

  free(sample_id);
  delete[] haplotype[0];
  delete[] haplotype[1];
  delete[] subpop[0];
  delete[] subpop[1];
  snps = NULL;
  n_snps = 0;
}

void Sample::meiosis(int8_t *gamate, int8_t *gamate_subpop) {
  int sh = rand()/(RAND_MAX + 1.0) < 0.5 ? 0 : 1;

  int i = 0;
  while(i < n_snps) {
    double tmp = rand()/(RAND_MAX + 1.0);
    double dsb = -log(tmp)*100. + snps[i].genetic_pos;
    int j = i;
    while(j < n_snps && snps[j].genetic_pos < dsb) j++;

    for(int k=i; k < j; k++) {
      gamate[k] = this->haplotype[sh][k];
      gamate_subpop[k] = this->subpop[sh][k];
    }
    
    sh = sh ^ 0x1;
    i = j;
  }

}
	
	
  
