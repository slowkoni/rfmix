/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

#ifndef MD5RNG_H
#define MD5RNG_H
#include <stdint.h>

class md5rng {
 public:
  md5rng(uint32_t seed);
  ~md5rng();
  uint32_t uint32(uint32_t k2, uint32_t k3, uint32_t k4);
  int uniform_int(uint32_t k2, uint32_t k3, uint32_t k4, int l, int r);
  
 private:
  uint32_t key;
};

#endif
