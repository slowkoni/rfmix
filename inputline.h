#ifndef INPUTLINE_H
#define INPUTLINE_H

#include <stdio.h>
#include <pthread.h>

class Inputline {
 public:
  char *fname;
  int line_no;
  
  Inputline(char *fname);
  Inputline(char *fname, char *chm);
  ~Inputline();
  char *nextline();
  char *nextline(int return_copy);
  void pushback();
  
  
 private:
  FILE *f;
  int child_pid;
  char *input_buf; /* verbatim copy of input line that was read */
  char *mod_buf;   /* copy of input_buf returned that the user may modify */

  int alloc_length;
  int line_stored;
  pthread_mutex_t lock;

  FILE *open_gzip_read(char *fname);
  FILE *open_bcftools_read(char *fname, char *chm);
};

enum { INPUTLINE_NOCOPY = 0, INPUTLINE_RETURN_COPY = 1 };

#endif
