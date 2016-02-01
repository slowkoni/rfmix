#ifndef S_SUBPOP_H
#define S_SUBPOP_H
#include "hash-table.h"

class Subpop {
 public:
  Subpop(char *name);
  ~Subpop(void);
  void add_sample(char *sample_name);
  char *subpop_name;
  int idx;

  static int n_subpops;
  static Subpop *lookup_sample_subpop(char *sample_name);
  static Subpop *lookup_subpop(char *subpop_name);
 private:

  static HashTable *subpop_map;
  static HashTable *sample_map;
};

#endif
