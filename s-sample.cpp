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

#include <stdint.h>
#include <math.h>

#include "kmacros.h"
#include "rfmix.h"
#include "s-sample.h"

HashTable *S_Sample::map = new HashTable(8192);

S_Sample::S_Sample(char *sample_id, int subpop, snp_t *snps, int n_snps, int8_t *h1, int8_t *h2) {

  this->sample_id = new char[strlen(sample_id)+1];
  strcpy(this->sample_id, sample_id);
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
  
  map->insert(sample_id, (void *) this);
}

S_Sample::S_Sample(S_Sample *p1, S_Sample *p2) {
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
  map->insert(sample_id, (void *) this);
}

S_Sample::~S_Sample(void) {

  map->remove(sample_id);
  delete[] sample_id;
  delete[] haplotype[0];
  delete[] haplotype[1];
  delete[] subpop[0];
  delete[] subpop[1];
  snps = NULL;
  n_snps = 0;
}

void S_Sample::meiosis(int8_t *gamate, int8_t *gamate_subpop) {
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
 
S_Sample *S_Sample::get_sample(char *sample_name) {
  return (S_Sample *) map->lookup(sample_name);
}
