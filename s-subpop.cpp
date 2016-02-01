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


