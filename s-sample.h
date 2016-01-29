#ifndef S_SAMPLE_H
#define S_SAMPLE_H
#include <stdint.h>
#include "vcf.h"

class Sample {
 public:
  Sample(char *sample_id, int subpop, snp_t *snps, int n_snps, int8_t *haplotype1, int8_t *haplotype2);
  Sample(Sample *p1, Sample *p2);
  ~Sample();
  
 private:
  char *sample_id;
  int n_snps;
  int8_t *haplotype[2];
  int8_t *subpop[2];

  /* Takes itself, a diploid sample, and produces a single haplotype and corresponding
     array tracking origin of each allele in the haplotype, simulating recombination
     between this sample's haplotypes as occurs during meiosis */
  void meiosis(int8_t *gamate, int8_t *gamate_subpop);
    
  /* Not owned by this class */
  snp_t *snps;
};
#endif
