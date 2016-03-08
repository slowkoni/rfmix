/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

#ifndef S_SAMPLE_H
#define S_SAMPLE_H
#include <stdint.h>
#include "vcf.h"
#include "hash-table.h"

class S_Sample {
 public:
  char *sample_id;
  S_Sample(char *sample_id, int subpop, snp_t *snps, int n_snps, int8_t *haplotype1, int8_t *haplotype2);
  S_Sample(S_Sample *p1, S_Sample *p2);
  ~S_Sample();
  static S_Sample *get_sample(char *sample_name);
  int8_t *haplotype[2];
  int8_t *subpop[2];
  
  
 private:
  int n_snps;

  /* Takes itself, a diploid sample, and produces a single haplotype and corresponding
     array tracking origin of each allele in the haplotype, simulating recombination
     between this sample's haplotypes as occurs during meiosis */
  void meiosis(int8_t *gamate, int8_t *gamate_subpop);
    
  /* Not owned by this class */
  snp_t *snps; 
  static HashTable *map;
  static int n_samples;
};
#endif
