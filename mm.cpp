/* (c) 2007-2016 Mark Hamilton Wright */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <stdint.h>
#include <unistd.h>

#include "kmacros.h"
#include "mm.h"

/* implementation of a scatter/gather memory allocator -- see mm.h for 
   details */

#define BLOCK_LIST_ALLOC_STEP (32)
mm::mm(int block_size, WHEREARGS) {

  this->block_size = block_size*1024*1024;

  n_blocks = 1;
  MA(blocks, sizeof(void *)*BLOCK_LIST_ALLOC_STEP, void *);
  MA(blocks[0], sizeof(char)*this->block_size, char);
  current_block = 0;
  ptr = 0;
}

void *mm::allocate(int size, WHEREARGS) {
  void *p;
  size_t asize;

  if (size <= 0) {
    fprintf(stderr,"mm::allocate() called at %s():%d (%s) requesting zero or negative "
	    "size (%d)\n", callfunc, callline, callfile, size);
    return NULL;
  }
  /* keep alignment to 8 bytes for improved memory bus access speed on most
     architectures. For most architectures this may be required otherwise
     access to types allocated through this interface that are not of type
     char may generate bus errors. 8 bytes should work for any 64-bit 
     architecture, but x86_64 may work just as well with alignment to every
     4 bytes */
  asize = (size & 0x7) ? (size | 0x7) + 1 : size;

#ifdef DEBUG
  asize += sizeof(uint64_t)*2;
  if (asize > block_size) {
    fprintf(stderr,"Error: %s() (%s:%d) requesting %u bytes plus alignment "
	    "padding (%d bytes) exceeds mm block size (%d bytes)\n", callfunc,
	    callfile, callline, size, asize - size, block_size);
  }
#endif

  if (ptr + asize >= block_size) {
    current_block++;
    if (current_block == n_blocks) {
      if (n_blocks % BLOCK_LIST_ALLOC_STEP == 0)
	RA(blocks, sizeof(void *)*(BLOCK_LIST_ALLOC_STEP + n_blocks), void *);

      MA(blocks[n_blocks], sizeof(char)*block_size, char);
      n_blocks++;
    }

    ptr = 0;
  }

  p = (void *) ((uint64_t) blocks[current_block] + (uint64_t) ptr);
  ptr += asize;
  
#ifdef DEBUG 
# warning DEBUG path in mm.c enabled
  *((uint64_t *) p) = LEADING_SENTINAL;
  p += sizeof(uint64_t);
  
  p[np] = p;
  s[np] = size;
  
  /* place trailing sentinal directly after the allocated memory, which may
     not be uint64_t aligned */
  for(i=0;i<sizeof(uint64_t);i++)
    p[size+i] = (TRAILING_SENTINAL >> (i*8)) & 0xFF;  
#endif

  return p;
}

/* resets the allocation pointer for mm without actually freeing memory, unless
   the blocksize is adjusted or auto-tuned by repeated calling of this 
   function. After several calls, if being used again and again in the same 
   loop context, this function will simply reset pointers and not call free() */
void mm::recycle() {
  int __attribute__((unused))i;
  size_t __attribute__((unused))new_blocksize;

  /* if more than one block was allocated in last usage, increase the
     blocksize so that a single block would have been enough */
#if 1
  if (n_blocks > 1) {
    new_blocksize = block_size * n_blocks;
    //fprintf(stderr,"Set new blocksize to %1.1f Mb\n", new_blocksize/(double) (1024*1024));
    /* don't allow a blocksize increment to bring it over 16 Mb. Otherwise,
       make sure the new_blocksize is a multiple of the system pagesize */
    if (new_blocksize > 64*1024*1024) {
      new_blocksize = 64*1024*1024;
    } else {
      new_blocksize = (new_blocksize | (getpagesize()-1)) + 1;
    }
    //fprintf(stderr,"Set new blocksize to %1.1f Mb\n", new_blocksize/(double) (1024*1024));

    if (n_blocks > 1) {
      for(i=1;i<n_blocks;i++) free(blocks[i]);
      
      RA(blocks[0], sizeof(char)*new_blocksize, char);
      block_size = new_blocksize;
      n_blocks = 1;
    }
  }
#endif
  
  current_block = 0;
  ptr = 0;
}

/* permanently done with using mm for thread specific memory allocation -- 
   release all allocated memory at once. The main purpose of this system is
   so that many many small objects can be allocated as pointers to within 
   large pre-allocated blocks, and then *all* of the pointers can be "freed"
   with a single (or very few) calls to free() */
mm::~mm() {
  int i;

  for(i=0; i < n_blocks;i++) {
    free(blocks[i]);
    blocks[i] = NULL;
  }
  free(blocks);
  
  n_blocks = -1;
  current_block = -1;
  ptr = -1;
}
