#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <kmacros.h>
#include "inputline.h"

#define INPUTLINE_CHUNK (8192)
Inputline::Inputline(char *fname) {

  f = fopen(fname, "r");
  if (f == NULL) {
    fprintf(stderr,"\nCan't open input file %s (%s)\n\n", fname, strerror(errno));
    exit(-1);
  }

  MA(input_buf, INPUTLINE_CHUNK, char);
  MA(mod_buf, INPUTLINE_CHUNK, char);
  alloc_length = INPUTLINE_CHUNK;
  line_stored = 0;
  line_no = 0;

  pthread_mutex_init(&lock, NULL);
}

char *Inputline::nextline(int return_copy) {
  int l, input_l;
  char *rval;

  pthread_mutex_lock(&lock);
  if (line_stored) {
    line_stored = 0;

    if (return_copy) 
      rval = strdup(input_buf);
    else {
      strcpy(mod_buf, input_buf);
      rval = mod_buf;
    }

    pthread_mutex_unlock(&lock);
    return rval;
  }

  l = 0;
  input_l = 0;
  input_buf[0] = 0;
  line_no++;
  for(;;) {
    if (alloc_length - l < INPUTLINE_CHUNK) {
      alloc_length += INPUTLINE_CHUNK;
      RA(input_buf, alloc_length, char);
      RA(mod_buf, alloc_length, char);
    }

    if (fgets(input_buf + l, INPUTLINE_CHUNK, f) == NULL) break;

    input_l = strlen(input_buf + l);
    if (input_l < INPUTLINE_CHUNK - 1) break;
    l += input_l;
  }

  if (input_l == 0) {
    rval = NULL;
  } else {
    if (return_copy) 
      rval = strdup(input_buf);
    else {
      strcpy(mod_buf, input_buf);
      rval = mod_buf;
    }
  }
  
  pthread_mutex_unlock(&lock);
  return rval;
}

void Inputline::pushback() {
  /* it probably does not make sense for inputline_pushback() to be used in
     a multithreaded context */
  pthread_mutex_lock(&lock);
  line_stored = 1;
  pthread_mutex_unlock(&lock);
}

Inputline::~Inputline() {

  if (pthread_mutex_destroy(&lock) == EBUSY)
    fprintf(stderr,"Warning: Call to inputline_free() with object locked\n");

  fclose(f);
  f = NULL;
  
  free(input_buf);
  input_buf = NULL;

  free(mod_buf);
  mod_buf = NULL;

  alloc_length = 0;
  line_stored = 0;
}
