/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "kmacros.h"
#include "hash-table.h"
#include "s-subpop.h"

HashTable *Subpop::sample_map = new HashTable(8192);
HashTable *Subpop::subpop_map = new HashTable(8192);
int Subpop::n_subpops = 0;

#define SUBPOP_ALLOC_STEP (8)
Subpop::Subpop(char *subpop_name) {
  if (subpop_map->lookup(subpop_name) != NULL) {
    fprintf(stderr,"Error: attempt to create a new subpop with the sample name as an existing one\n");
    exit(-1);
  }
  
  this->subpop_name = strdup(subpop_name);
  this->idx = Subpop::n_subpops++;
  subpop_map->insert(subpop_name, this);
}

Subpop::~Subpop(void) {
  subpop_map->remove(subpop_name);
}

Subpop *Subpop::lookup_subpop(char *name) {
  return (Subpop *) subpop_map->lookup(name);
}

Subpop *Subpop::lookup_sample_subpop(char *sample_name) {
  return (Subpop *) sample_map->lookup(sample_name);
}

void Subpop::add_sample(char *sample_name) {
  Subpop::sample_map->insert(sample_name, this);
}

static void selection_sort(char **a, int n) {

  for(int i=0; i < n; i++) {
    int m = i;
    for(int j=i+1; j < n; j++)
      if (strcmp(a[j],a[m]) < 0) m = j;
    if (m != i) {
      char *tmp = a[i];
      a[i] = a[m];
      a[m] = tmp;
    }
  }
}

void Subpop::set_ordering(void) {
  int n_subpops = subpop_map->n_items();
  char **subpop_names = new char*[n_subpops];
  
  int i = 0;
  char *key;
  subpop_map->reset();
  while((key = subpop_map->nextkey()) && i < n_subpops) {
    subpop_names[i++] = key;
  }

  selection_sort(subpop_names, n_subpops);
  for(i=0; i < n_subpops; i++) {
    Subpop *s = lookup_subpop(subpop_names[i]);
    s->idx = i;
  }
  
  
  delete[] subpop_names;
}
    
      


