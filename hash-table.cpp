#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <time.h>
#include <stdint.h>

#include "hash-table.h"

/* HashTable - simple utility class for a hash table. No doubt any C++ programmer 
   looking at this file will wonder why the hell anyone would implement a class
   this way. What this really is is a C implementation of a hash table "class"
   used in many other C programs I have written over the years, that I adapted
   for C++ simply by leaving the C implementation mostly as it was but changing
   all of its exported functions to static, file-scope only functions, and then
   writing a C++ class that wraps the C implementation with public-facing 
   functions. It does not use fancy C++ things like containers and just returns
   void * pointers to whatever memory was passed as the data to store. As with
   use in C, it is the responsibility of the user to free/delete those objects.
   If the HashTable object is delete'd, it will free its internal memory usage
   but will not call destructors or free() on the data remaining in the hash
   table. */

static uint32_t random_tab[256];

/* This prime table was lifted from the compress.c LZW compression program,
   these days typically distributed as package called "ncompress" */
static int primetab[256] = {
  1013, -1061, 1109, -1181, 1231, -1291, 1361, -1429,
  1481, -1531, 1583, -1627, 1699, -1759, 1831, -1889,
  1973, -2017, 2083, -2137, 2213, -2273, 2339, -2383,
  2441, -2531, 2593, -2663, 2707, -2753, 2819, -2887,
  2957, -3023, 3089, -3181, 3251, -3313, 3361, -3449,
  3511, -3557, 3617, -3677, 3739, -3821, 3881, -3931,
  4013, -4079, 4139, -4219, 4271, -4349, 4423, -4493,
  4561, -4639, 4691, -4783, 4831, -4931, 4973, -5023,
  5101, -5179, 5261, -5333, 5413, -5471, 5521, -5591,
  5659, -5737, 5807, -5857, 5923, -6029, 6089, -6151,
  6221, -6287, 6343, -6397, 6491, -6571, 6659, -6709,
  6791, -6857, 6917, -6983, 7043, -7129, 7213, -7297,
  7369, -7477, 7529, -7577, 7643, -7703, 7789, -7873,
  7933, -8017, 8093, -8171, 8237, -8297, 8387, -8461,
  8543, -8627, 8689, -8741, 8819, -8867, 8963, -9029,
  9109, -9181, 9241, -9323, 9397, -9439, 9511, -9613,
  9677, -9743, 9811, -9871, 9941,-10061,10111,-10177,
  10259,-10321,10399,-10477,10567,-10639,10711,-10789,
  10867,-10949,11047,-11113,11173,-11261,11329,-11423,
  11491,-11587,11681,-11777,11827,-11903,11959,-12041,
  12109,-12197,12263,-12343,12413,-12487,12541,-12611,
  12671,-12757,12829,-12917,12979,-13043,13127,-13187,
  13291,-13367,13451,-13523,13619,-13691,13751,-13829,
  13901,-13967,14057,-14153,14249,-14341,14419,-14489,
  14557,-14633,14717,-14767,14831,-14897,14983,-15083,
  15149,-15233,15289,-15359,15427,-15497,15583,-15649,
  15733,-15791,15881,-15937,16057,-16097,16189,-16267,
  16363,-16447,16529,-16619,16691,-16763,16879,-16937,
  17021,-17093,17183,-17257,17341,-17401,17477,-17551,
  17623,-17713,17791,-17891,17957,-18041,18097,-18169,
  18233,-18307,18379,-18451,18523,-18637,18731,-18803,
  18919,-19031,19121,-19211,19273,-19381,19429,-19477
};

static int init_randomtab = 1;

#define get16bits(d) (*((const uint16_t *) (d)))
static uint32_t superfasthash (const char * data, int len) {
  uint32_t hash = len, tmp;
  int rem;

  if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}


#define F(X,Y,Z) (((X) & (Y)) ^ ((~(X)) & (Z)))
#define G(X,Y,Z) (((X) & (Y)) | ((Y) & (Z)) | ((X) & (Z)))
#define H(X,Y,Z) ((X) ^ (Y) ^ (Z))
#define ROTATE(x, s) (((x) << (s)) | ((x) >> (32 - (s))))
#define RD1(a,b,c,d,x,s) (a) = ROTATE(((a) + F((b),(c),(d)) + (x)), (s)) 
#define RD2(a,b,c,d,x,s) (a) = ROTATE(((a)+G((b),(c),(d))+(x)+0x5A827999), (s))
#define RD3(a,b,c,d,x,s) (a) = ROTATE(((a)+H((b),(c),(d))+(x)+0x6ED9EBA1), (s))
static uint32_t __attribute__((unused))md4(char *m, uint64_t length) {
  uint32_t a, b, c, d;
  uint64_t pad_length, i;
  
  a = 0x67452301;
  b = 0xEFCDAB89;
  c = 0x98BADCFE;
  d = 0x10325476;

  pad_length = length + 2 + (64 - (length + 2) % 64);
  char pad_m[pad_length];
  //  pad_m = alloca(pad_length);
  for(i=0;i<length;i++) {
    pad_m[i] = m[i];
  }
  pad_m[i++] = 0x80;
  while(i < pad_length - 2)
    pad_m[i++] = 0x00;

  memcpy(pad_m + i, (char *) &length, sizeof(uint64_t));

  for(i=0;i<length/16 + 1;i++) {
    uint32_t aa, bb, cc, dd;
    char *p;
    aa = a;
    bb = b;
    cc = c;
    dd = d;
    p = pad_m + i*16;

    RD1(a, b, c, d, p[0], 3);
    RD1(d, a, b, c, p[1], 7);
    RD1(c, d, a, b, p[2], 11);
    RD1(b, c, d, a, p[3], 19);

    RD1(a, b, c, d, p[4], 3);
    RD1(d, a, b, c, p[5], 7);
    RD1(c, d, a, b, p[6], 11);
    RD1(b, c, d, a, p[7], 19);

    RD1(a, b, c, d, p[8], 3);
    RD1(d, a, b, c, p[9], 7);
    RD1(c, d, a, b, p[10], 11);
    RD1(b, c, d, a, p[11], 19);

    RD1(a, b, c, d, p[12], 3);
    RD1(d, a, b, c, p[13], 7);
    RD1(c, d, a, b, p[14], 11);
    RD1(b, c, d, a, p[15], 19);

    RD2(a, b, c, d, p[0], 3);
    RD2(d, a, b, c, p[4], 5);
    RD2(c, d, a, b, p[8], 9);
    RD2(b, c, d, a, p[12], 13);

    RD2(a, b, c, d, p[1], 3);
    RD2(d, a, b, c, p[5], 5);
    RD2(c, d, a, b, p[9], 9);
    RD2(b, c, d, a, p[13], 13);

    RD2(a, b, c, d, p[2], 3);
    RD2(d, a, b, c, p[6], 5);
    RD2(c, d, a, b, p[10], 9);
    RD2(b, c, d, a, p[14], 13);

    RD2(a, b, c, d, p[3], 3);
    RD2(d, a, b, c, p[7], 5);
    RD2(c, d, a, b, p[11], 9);
    RD2(b, c, d, a, p[15], 13);

    RD3(a, b, c, d, p[0], 3);
    RD3(d, a, b, c, p[8], 9);
    RD3(c, d, a, b, p[4], 11);
    RD3(b, c, d, a, p[12], 15);

    RD3(a, b, c, d, p[2], 3);
    RD3(d, a, b, c, p[10], 9);
    RD3(c, d, a, b, p[6], 11);
    RD3(b, c, d, a, p[8], 15);

    RD3(a, b, c, d, p[1], 3);
    RD3(d, a, b, c, p[9], 9);
    RD3(c, d, a, b, p[5], 11);
    RD3(b, c, d, a, p[13], 15);

    RD3(a, b, c, d, p[3], 3);
    RD3(d, a, b, c, p[11], 9);
    RD3(c, d, a, b, p[7], 11);
    RD3(b, c, d, a, p[15], 15);

    a = a + aa;
    b = b + bb;
    c = c + cc;
    d = d + dd;
  }

  return a;
}

inline uint32_t simple_hash(char *key, int key_length, int mode) {
  uint32_t hash;
  int i, j;

  if (mode == 0) {
    hash = (key_length*5 << 24) | (key_length*3 << 16) | (key_length*3 << 8) 
      | key_length;
    for(i=0;i<key_length;i++) {
      j = random_tab[(unsigned) key[i]] & 0x3;
      hash |= (random_tab[((hash >> j*8) & 0xFF) ^ key[i]] << j*8);
      hash = (hash << 4) | (hash >> 28);
    }
  } else {
    hash = key_length;
    for(i=0;i<key_length;i++) {
      j = random_tab[(unsigned) key[i]] & 0x3;
      hash |= (random_tab[((hash >> j*8) & 0xFF) ^ key[i]] << j*8);
      hash = (hash >> 4) | (hash << 28);
    }
  }

  return hash;
}

inline uint32_t hashfunc(char *key, int key_length, int ht_size, int mode) {
  uint32_t hash;

  //hash = simple_hash(key, key_length, 0) & (ht_size - 1);
  //  hash = md4(key, key_length) & (ht_size - 1);
  hash = superfasthash(key, key_length) & (ht_size - 1);
  return hash;
}


static hash_table_t *new_hash_table(int initial_size, int random_seed) {
  hash_table_t *ht;
  int i;

  if (init_randomtab) {
    if (random_seed < 0) {
      srandom(time(NULL));
    } else if (random_seed > 0) {
      srandom(random_seed);
    }

    for(i=0;i<256;i++) {
      random_tab[i] = (random() >> 23);
    }
    init_randomtab = 0;
  }

  /* If initial_size is not a power of two, make it the next largest power
     of two */
  if ((initial_size & (initial_size - 1)) != 0) {
    i = 0x1;
    while(i < initial_size) {
      i <<= 1;
    }
    initial_size = i;
  }

  MA(ht, sizeof(hash_table_t), hash_table_t);
  ht->size = initial_size;
  MA(ht->table, sizeof(hash_entry_t)*initial_size, hash_entry_t);
  for(i=0;i<initial_size;i++) {
    ht->table[i].key = NULL;
    ht->table[i].key_length = 0;
    ht->table[i].item = NULL;
  }

  ht->n_ks = 0;
  ht->key_space = NULL;
  ht->free_keyspace = NULL;
  ht->kp = NULL;

  ht->occupied = 0;
  ht->tp = 0;
  return ht;
}

inline int key_compare(hash_entry_t *ent, char *key, int key_length) {
  int i;

  if (ent == NULL || key == NULL) return 1;
  if (ent->key_length == -1) return -1;
  if (ent->key_length == 0) return -2;

  if (ent->key_length != key_length) return 1;

  i = 0;
  while(i<key_length) {
    if (ent->key[i] != key[i]) return 1;
    i++;
  }
  return 0;
}

inline hash_entry_t *ht_search(hash_table_t *ht, char *key, int key_length,
			       int skip_deleted) {
  uint32_t k1, k2;
  int c;

  k1 = hashfunc(key, key_length, ht->size, 0);
  c = key_compare(ht->table + k1, key, key_length);
  if (c == 0 || (c == -1 && !skip_deleted) || c == -2) return ht->table + k1;

  k2 = primetab[hashfunc(key + 1, key_length - 1, ht->size, 0) & 0xFF];
  do {
    k1 = (k1 + k2) & (ht->size - 1);
    c = key_compare(ht->table + k1, key, key_length);
  } while(c > 0 || (skip_deleted && c==-1));

  return ht->table + k1;
}

#define KEYSPACE_SIZE (1024*1024)
static char *allocate_keyspace(hash_table_t *ht, char *key, int key_length) {
  int i;
  size_t size;
  char *p;

  i = 0;
  while(i<ht->n_ks && ht->free_keyspace[i] < key_length+1) i++;
  if (i == ht->n_ks) {
    if (i % 8 == 0) {
      RA(ht->key_space, sizeof(char *)*(ht->n_ks + 8), char *);
      RA(ht->kp, sizeof(int)*(ht->n_ks + 8), int);
      RA(ht->free_keyspace, sizeof(int)*(ht->n_ks + 8), int);
    }

    size = key_length > KEYSPACE_SIZE ? key_length + 1 : KEYSPACE_SIZE;
    MA(ht->key_space[i], sizeof(char)*size, char);
    ht->free_keyspace[i] = size;
    ht->kp[i] = 0;
    ht->n_ks++;
  }

  p = ht->key_space[i] + ht->kp[i];
  ht->kp[i] += key_length + 1;
  ht->free_keyspace[i] -= key_length + 1;
  return p;
}


static void ht_free(hash_table_t *ht) {
  int i;

  free(ht->table);
  for(i=0;i<ht->n_ks;i++) {
    free(ht->key_space[i]);
    ht->key_space[i] = NULL;
  }
  free(ht->key_space);
  free(ht->free_keyspace);
  free(ht->kp);
  free(ht);
}


static void rebuild_table(hash_table_t *ht, int new_size);
static void ht_insert(hash_table_t *ht, char *key, int key_length,
		      void *data) {
  hash_entry_t *ent;
  
  if (ht->occupied*2 >= ht->size) rebuild_table(ht, ht->size*2);

  ent = ht_search(ht, key, key_length, 0);
  if (ent->key) {
    ent->item = data;
    return;
  }
  ent->key = allocate_keyspace(ht, key, key_length);
  memcpy(ent->key, key, key_length);
  ent->key[key_length] = 0;
  ent->key_length = key_length;

  ent->item = data;

  ht->occupied ++;
  ht->tp = 0;
}

static void ht_delete(hash_table_t *ht, char *key, int key_length) {
  hash_entry_t *ent;

  ent = ht_search(ht, key, key_length, 1);
  if (ent->key_length == 0) return;

  ent->key = NULL;
  ent->key_length = -1;
  ent->item = NULL;

  ht->occupied--;
  if (ht->occupied * 4 < ht->size) rebuild_table(ht, ht->size/2);
}

void *ht_lookup(hash_table_t *ht, char *key, int key_length) {
  hash_entry_t *ent;

  ent = ht_search(ht, key, key_length, 1);
  if (ent->key_length <= 0) return NULL;
  return ent->item;
}

static void rebuild_table(hash_table_t *ht, int new_size) {
  int i;
  hash_table_t *new_ht;

  new_ht = new_hash_table(new_size, 0);
  for(i=0;i<ht->size;i++) {
    if (ht->table[i].key_length > 0) {
      ht_insert(new_ht, ht->table[i].key, ht->table[i].key_length, 
		ht->table[i].item);
    }
  }

  free(ht->table);
  for(i=0;i<ht->n_ks;i++) 
    free(ht->key_space[i]);
  free(ht->key_space);
  free(ht->kp);
  free(ht->free_keyspace);

  ht->table = new_ht->table;
  ht->size = new_size;

  ht->n_ks = new_ht->n_ks;
  ht->key_space = new_ht->key_space;
  ht->kp = new_ht->kp;
  ht->free_keyspace = new_ht->free_keyspace;

  free(new_ht);
}

static void __attribute__((unused))ht_reset(hash_table_t *ht) {
  ht->tp = 0;
}

static char __attribute__((unused))*ht_nextkey(hash_table_t *ht, int *key_length) {

  while(ht->tp < ht->size && ht->table[ht->tp].key_length <= 0) ht->tp++;

  if (ht->tp == ht->size) {
    ht->tp = 0;
    if (key_length != NULL) *key_length = 0;
    return NULL;
  }
  if (key_length != NULL) *key_length = ht->table[ht->tp].key_length;
  return ht->table[ht->tp++].key;
}

static void __attribute__((unused))*ht_nextitem(hash_table_t *ht, char **key, int *key_length) {
  while(ht->tp < ht->size && ht->table[ht->tp].key_length <= 0) ht->tp++;

  if (ht->tp == ht->size) {
    ht->tp = 0;
    if (key != NULL) *key = NULL;
    if (key_length != NULL) *key_length = 0;
    return NULL;
  }

  if (key != NULL) *key = ht->table[ht->tp].key;
  if (key_length != NULL) *key_length = ht->table[ht->tp].key_length;
  return ht->table[ht->tp++].item;
}

static double ht_collision_ratio(hash_table_t *ht) {
  int k1, k2;
  int collisions = 0;
  int i;
  hash_entry_t *ent;

  for(i=0;i<ht->size;i++) {
    if (ht->table[i].key_length > 0) {
      ent = ht->table + i;
      k1 = hashfunc(ent->key, ent->key_length, ht->size, 0);
      k2 = primetab[hashfunc(ent->key+1, ent->key_length-1, ht->size, 0)&0xFF];
      while(k1 != i) {
	collisions++;
	k1  = (k1 + k2) & (ht->size - 1);
      }
    }
  }
  /* average number of hash collisions per key in table */
  return collisions/(double) ht->occupied;
}

HashTable::HashTable(int initial_size, int random_seed) {
  ht = new_hash_table(initial_size, random_seed);
}

HashTable::HashTable(int initial_size) {
  ht = new_hash_table(initial_size, time(NULL));
}

HashTable::HashTable(void) {
  ht = new_hash_table(16, time(NULL));
}

HashTable::~HashTable(void) {
  ht_free(ht);
}

void HashTable::insert(char *key, void *data) {
  ht_insert(ht, key, strlen(key), data);
}

void *HashTable::lookup(char *key) {
  return ht_lookup(ht, key, strlen(key));
}

void HashTable::remove(char *key) {
  return ht_delete(ht, key, strlen(key));
}

double HashTable::collision_ratio(void) {
  return ht_collision_ratio(ht);
}

char *HashTable::nextkey(void) {
  int unused;
  return ht_nextkey(ht, &unused);
}

void *HashTable::nextitem(void) {
  char *key;
  int key_length;
  return ht_nextitem(ht, &key, &key_length);
}

void HashTable::reset(void) {
  ht_reset(ht);
}

int HashTable::n_items(void) {
  return ht->occupied;
}
