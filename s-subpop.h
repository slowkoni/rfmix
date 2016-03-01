/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

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
  static void set_ordering(void);
 private:

  static HashTable *subpop_map;
  static HashTable *sample_map;
};

#endif
