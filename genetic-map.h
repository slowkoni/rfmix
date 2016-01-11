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
