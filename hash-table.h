/* (c) 2008-2016 Mark Hamilton Wright */
#ifndef HASH_TABLE_H
#define HASH_TABLE_H
#include "kmacros.h"

typedef struct {
  char *key;
  int key_length;
  void  *item;
} hash_entry_t;

typedef struct {
  int size;
  hash_entry_t *table;
  int occupied;
  int tp;

  int n_ks;
  char **key_space;
  int *kp;
  int *free_keyspace;
} hash_table_t;

/*
hash_table_t *new_hash_table(int initial_size, int random_seed);
void ht_insert(hash_table_t *ht, char *key, int key_length, void *data);
void ht_delete(hash_table_t *ht, char *key, int key_length);
void *ht_lookup(hash_table_t *ht, char *key, int key_length);
char *ht_nextkey(hash_table_t *ht, int *key_length);
void *ht_nextitem(hash_table_t *ht, char **key, int *key_length);
double ht_collision_ratio(hash_table_t *ht);
void ht_reset(hash_table_t *ht);
void ht_free(hash_table_t *ht);
*/

class HashTable {
 public:
  HashTable();
  HashTable(int initial_size);
  HashTable(int initial_size, int random_seed);
  ~HashTable();

  void insert(char *key, void *data);
  void *lookup(char *key);
  void remove(char *key);
  double collision_ratio();
  void reset();
  char *nextkey();
  void *nextitem();
  int n_items();
  
 private:
  hash_table_t *ht;
};

#endif
