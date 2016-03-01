/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

#ifndef GENETIC_MAP_H
#define GENETIC_MAP_H

typedef struct {
  int seq_pos;
  double genetic_pos;
} map_pos_t;

class GeneticMap {
  char *chm;
  map_pos_t *map;
  int n_pos;

public:
  GeneticMap(void);
  ~GeneticMap(void);
  void load_map(char *fname, char *chm);
  double translate_seqpos(int seq_pos);

private:
  int binary_search(int seq_pos);
};
#endif
