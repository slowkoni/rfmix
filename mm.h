/* (c) 2007-2016 Mark Hamilton Wright */
#ifndef KONI_MM_H
#define KONI_MM_H

#include <kmacros.h>

/* this is a simple memory management interface for common cases where 
   sequential operations are expected to want to call malloc() or something
   to allocate small chunks of memory, but would allocate very many of them
   but how many can not be known in advance, and the memory allocated will
   either never be freed until the program is finished or it will be freed
   all at once.

   Memory is allocated via mmap() in much larger blocks, and calls to mm_a() 
   are made to allocate the small chunks as needed. mm_a() return pointers 
   within blocks that are not freed individually. As needed, more blocks
   are allocated but not via realloc() so as to avoid large amounts of 
   memory being. When the program is finished with *all* of these little
   bits of memory, allocated blocks are freed and all pointers that were
   allocated via mm_a() are now invalid. */
#if 0
typedef struct {
  int n_blocks;
  void **blocks;
  size_t block_size;
  int current_block;
  int ptr;
} mm_block_t;


/* NOTE!! blocksize is in MB */
mm_block_t *mm_block_init(int block_size, WHEREARGS);

void *mm_a(mm_block_t *mm, int size, WHEREARGS);

void mm_recycle(mm_block_t *mm);
void mm_free(mm_block_t *mm);


#endif

#ifndef WHEREARGS
#define WHEREARGS const char *callfunc, const char *callfile, int callline
#define WHEREFROM __func__, __FILE__, __LINE__
#endif

class mm {
 public:
  mm(int mb_blocksize, WHEREARGS);
  ~mm();
  void *allocate(int size, WHEREARGS);
  void recycle();

 private:
  int n_blocks;
  void **blocks;
  size_t block_size;
  int current_block;
  size_t ptr;
};

#endif

