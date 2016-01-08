#ifndef INPUTLINE_H
#define INPUTLINE_H

#include <stdio.h>
#include <pthread.h>

typedef struct {
  FILE *f;

  char *input_buf; /* verbatim copy of input line that was read */
  char *mod_buf;   /* copy of input_buf returned that the user may modify */

  int alloc_length;
  int line_stored;
  int line_no;
  pthread_mutex_t lock;
} inputline_t;

enum { INPUTLINE_NOCOPY = 0, INPUTLINE_RETURN_COPY = 1 };

#ifdef __cplusplus
extern "C" {
#endif

inputline_t *init_inputline(FILE *f);
void inputline_free(inputline_t *obj);

char *inputline_nextline(inputline_t *obj, int return_copy);

/* the user may "unread" a line by calling this function. It will cause the
   next call to inputline_nextline() to copy the unmodified version of the
   last read line in input_buf[] to mod_buf[] and return mod_buf[] (or a copy
   if return_copy is set) to the user without reading anything new from the
   file */
void inputline_pushback(inputline_t *obj);
#ifdef __cplusplus
}
#endif
#endif
