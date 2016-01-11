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
