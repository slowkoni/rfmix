/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>

#include "genetic-map.h"
#include "kmacros.h"

static int map_pos_compare(const void *ta, const void *tb) {
  map_pos_t *a = (map_pos_t *) ta;
  map_pos_t *b = (map_pos_t *) tb;

  if (a->seq_pos < b->seq_pos) return -1;
  if (a->seq_pos > b->seq_pos) return  1;
  return 0;
}

GeneticMap::GeneticMap(void) {
  chm = NULL;
  n_pos = 0;
  map = NULL;
}

#define POS_ALLOC_STEP (8192)
#define LINE_MAX (8192)
void GeneticMap::load_map(char *fname, char *chm) {
  FILE *f = fopen(fname, "r");
  if (f == NULL) {
    fprintf(stderr,"\nCan't open genetic map file %s (%s)\n\n", fname, strerror(errno));
    exit(-1);
  }

  map_pos_t *tmp_map = NULL;
  int n_pos = 0;
  char inputline[LINE_MAX];
  while(fgets(inputline, LINE_MAX, f) != NULL) {
    CHOMP(inputline);
    if (inputline[0] == 0 || inputline[0] == '#') continue;
    
    char *p, *q;
    p = inputline;
    q = strsep(&p," \t");
    if (strcmp(q, chm) != 0 && (strncasecmp(q, "chr", 3) != 0 || strcmp(q+3, chm) != 0)) 	
      continue;

    if (n_pos % POS_ALLOC_STEP == 0)
      RA(tmp_map, n_pos + POS_ALLOC_STEP, map_pos_t);
    
    q = strsep(&p, " \t");
    tmp_map[n_pos].seq_pos = atoi(q);
    q = strsep(&p, " \t");
    tmp_map[n_pos].genetic_pos = atof(q);
    n_pos++;
  }
  fclose(f);
  
  if (n_pos == 0) {
    fprintf(stderr,"\nSTOP: no genetic map positions for chromosome %s were found in %s\n\n",
	    chm, fname);
    exit(-1);
  }
  if (n_pos == 1) {
    fprintf(stderr,"\nSTOP: Number of genetic map positions must be 2 or more, only 1 found\n\n");
    exit(-1);
  }
  
  qsort(tmp_map, n_pos, sizeof(map_pos_t), map_pos_compare);
  for(int i = 1; i < n_pos; i++) {
    if (tmp_map[i].genetic_pos < tmp_map[i-1].genetic_pos) {
      fprintf(stderr,"\nSTOP: genetic map for chromosome %s is not strictly increasing.\n\n",
	      chm);
      exit(-1);
    }
  }

  this->map = tmp_map;
  this->chm = strdup(chm);
  this->n_pos = n_pos;
}

double GeneticMap::translate_seqpos(int seq_pos) {

  int i = binary_search(seq_pos);
  /* fprintf(stderr,"Searched for %d\t search returns %d: %1.2f\tpos %d\tnext pos %d\n",
     seq_pos, i, map[i].genetic_pos, map[i].seq_pos,
     i < this->n_pos - 1 ? map[i+1].seq_pos : map[i].seq_pos); */
  
  if (i == 0 && seq_pos < map[i].seq_pos)
    return map[0].genetic_pos * (seq_pos /(double) map[0].seq_pos);
  if (i >= this->n_pos-1) i = this->n_pos - 2;
  
  double rate = (map[i+1].genetic_pos - map[i].genetic_pos)/(double)
      (map[i+1].seq_pos - map[i].seq_pos);
  return map[i].genetic_pos + rate*(seq_pos - map[i].seq_pos);
}

int GeneticMap::binary_search(int seq_pos) {

  int i = 0;
  int j = n_pos;
  int m;
  for(;;) {
    m = (i+j)/2;
    if (m == i) return i;
    if (map[m].seq_pos == seq_pos) return m;
    if (map[m].seq_pos < seq_pos) i = m;
    if (map[m].seq_pos > seq_pos) j = m;
  }
}

GeneticMap::~GeneticMap(void) {
  if (map) free(map);
  if (chm) free(chm);
  n_pos = 0;
}

#ifdef UNIT_TEST

int main(int argc, char *argv[]) {
  char *tmp_fname = strdup("tmp-genetic-map-test-XXXXXX");
  int tmp_fd = mkstemp(tmp_fname);
  FILE *tmpf = fdopen(tmp_fd,"w");
  
  fprintf(tmpf,"1\t1000000\t1.0000\n");
  fprintf(tmpf,"1\t2000000\t3.0000\n");
  fprintf(tmpf,"1\t3000000\t4.0000\n");
  fprintf(tmpf,"2\t1000000\t20.000\n");
  fprintf(tmpf,"2\t2000000\t21.000\n");
  fclose(tmpf);

  GeneticMap *test_map = new GeneticMap;
  test_map->load_map(tmp_fname, "1");
  unlink(tmp_fname);

  for(int pos = 0; pos < 4000000; pos += 500000) {
    fprintf(stdout,"1\t%d\t%1.5f\n", pos, test_map->translate_seqpos(pos));
  }
  
  delete test_map;
  fprintf(stdout,"Test complete\n");

  return 0;
}

#endif
