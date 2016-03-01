/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

#include "md5rng.h"
#include <stdint.h>

inline uint32_t F(uint32_t X, uint32_t Y, uint32_t Z) {
  return (X & Y) | ((~X) & Z);
}

inline uint32_t G(uint32_t X, uint32_t Y, uint32_t Z) {
  return (X & Z) | (Y & (~Z));
}

inline uint32_t H(uint32_t X, uint32_t Y, uint32_t Z) {
  return X ^ Y ^ Z;
}

inline uint32_t I(uint32_t X, uint32_t Y, uint32_t Z) {
  return Y ^ (X | (~Z));
}

#define ROT32(x,s) (((x) << (s)) | ((x) >> (32 - (s))))
inline void FF(uint32_t &a, uint32_t b, uint32_t c, uint32_t d,
	    uint32_t M, uint32_t s, uint32_t t) {

  uint32_t tmp;
  tmp = a + F(b,c,d) + M + t;
  a = b + ROT32(tmp, s);
}

inline void GG(uint32_t &a, uint32_t b, uint32_t c, uint32_t d,
	    uint32_t M, uint32_t s, uint32_t t) {
  uint32_t tmp;

  tmp = a + G(b,c,d) + M + t;
  a = b + ROT32(tmp, s);
}

inline void HH(uint32_t &a, uint32_t b, uint32_t c, uint32_t d,
	    uint32_t M, uint32_t s, uint32_t t) {
  uint32_t tmp;

  tmp = a + H(b,c,d) + M + t;
  a = b + ROT32(tmp, s);
}

inline void II(uint32_t &a, uint32_t b, uint32_t c, uint32_t d,
	    uint32_t M, uint32_t s, uint32_t t) {
  uint32_t tmp;

  tmp = a + I(b,c,d) + M + t;
  a = b + ROT32(tmp, s);
}

uint32_t md5(uint32_t k1, uint32_t k2, uint32_t k3, uint32_t k4) {
  uint32_t a, b, c, d;
  
  a = 0x01234567;
  b = 0x89ABCDEF;
  c = 0xFEDCBA98;
  d = 0x76543210;
  
  FF(a, b, c, d, k1, 7, 0xD76AA478);
  FF(d, a, b, c, k2, 12, 0xE8C7B756);
  FF(c, d, a, b, k3, 17, 0x242070DB);
  FF(b, c, d, a, k4, 22, 0xC1BDCEEE);

  GG(a, b, c, d, k1, 5,  0xF61E2562);
  GG(d, a, b, c, k2, 9,  0xC040B340);
  GG(c, d, a, b, k3, 14, 0x265E5A51);
  GG(b, c, d, a, k4, 20, 0xE9B6C7AA);
  
  HH(a, b, c, d, k1, 4,  0xFFFA3942);
  HH(d, a, b, c, k2, 11, 0x8771F681);
  HH(c, d, a, b, k3, 16, 0x6D9D6122);
  HH(b, c, d, a, k4, 23, 0xFDE5380C);

  II(a, b, c, d, k1, 6,  0xF4292244);
  II(d, a, b, c, k2, 10, 0x432AFF97);
  II(c, d, a, b, k3, 15, 0xAB9423A7);
  II(b, c, d, a, k4, 21, 0xFC93A039);

  return a ^ b ^ c ^ d;
}

md5rng::md5rng(uint32_t key) {
  this->key = key;
}

md5rng::~md5rng() {
}

uint32_t md5rng::uint32(uint32_t k2, uint32_t k3, uint32_t k4) {
  return md5(key, k2, k3, k4);
}

int md5rng::uniform_int(uint32_t k2, uint32_t k3, uint32_t k4, int l, int r) {
  double tmp = md5(key, k2, k3, k4)/(double) 0x100000000L;

  return l + (int) ((r - l)*tmp);
}

#ifdef MAIN
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int k1 = 8;
  int k2 = 16;
  int k3 = 1024;
  uint32_t key = 0xDEADBEEF;
  
  if (argc > 1) {
    key = atoi(argv[1]);
  }
  
  if (argc > 2) {
    k1 = atoi(argv[2]);
  }
  if (argc > 3) {
    k2 = atoi(argv[3]);
  }
  if (argc > 4) {
    k3 = atoi(argv[4]);
  }

  md5rng *rng = new md5rng(key);
  
  for(int i = 0; i < k1; i++) {
    for(int j = 0; j < k2; j++) {
      for(int k = 0; k < k3; k++) {
	uint32_t tmp = rng->uint32(i, j, k);
	fwrite((void *) &tmp, sizeof(uint32_t), 1, stdout);
      }
    }
  }

  return 0;
}
#endif
