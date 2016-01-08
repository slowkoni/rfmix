#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <kmacros.h>
#include "inputline.h"

#define INPUTLINE_CHUNK (8192)
inputline_t *init_inputline(FILE *f) {
  inputline_t *obj;

  MA(obj, 1, inputline_t);
  obj->f = f;

  MA(obj->input_buf, INPUTLINE_CHUNK, char);
  MA(obj->mod_buf, INPUTLINE_CHUNK, char);
  obj->alloc_length = INPUTLINE_CHUNK;
  obj->line_stored = 0;
  obj->line_no = 0;

  pthread_mutex_init(&obj->lock, NULL);
  return obj;
}

char *inputline_nextline(inputline_t *obj, int return_copy) {
  int l, input_l;
  char *rval;

  pthread_mutex_lock(&obj->lock);
  if (obj->line_stored) {
    obj->line_stored = 0;

    if (return_copy) 
      rval = strdup(obj->input_buf);
    else {
      strcpy(obj->mod_buf, obj->input_buf);
      rval = obj->mod_buf;
    }

    pthread_mutex_unlock(&obj->lock);
    return rval;
  }

  l = 0;
  input_l = 0;
  obj->input_buf[0] = 0;
  obj->line_no++;
  for(;;) {
    if (obj->alloc_length - l < INPUTLINE_CHUNK) {
      obj->alloc_length += INPUTLINE_CHUNK;
      RA(obj->input_buf, obj->alloc_length, char);
      RA(obj->mod_buf, obj->alloc_length, char);
    }

    if (fgets(obj->input_buf + l, INPUTLINE_CHUNK, obj->f) == NULL) break;

    input_l = strlen(obj->input_buf + l);
    if (input_l < INPUTLINE_CHUNK - 1) break;
    l += input_l;
  }

  if (input_l == 0) {
    rval = NULL;
  } else {
    if (return_copy) 
      rval = strdup(obj->input_buf);
    else {
      strcpy(obj->mod_buf, obj->input_buf);
      rval = obj->mod_buf;
    }
  }
  
  pthread_mutex_unlock(&obj->lock);
  return rval;
}

void inputline_pushback(inputline_t *obj) {
  /* it probably does not make sense for inputline_pushback() to be used in
     a multithreaded context */
  pthread_mutex_lock(&obj->lock);
  obj->line_stored = 1;
  pthread_mutex_unlock(&obj->lock);
}

void inputline_free(inputline_t *obj) {

  if (pthread_mutex_destroy(&obj->lock) == EBUSY)
    fprintf(stderr,"Warning: Call to inputline_free() with object locked\n");

  free(obj->input_buf);
  free(obj->mod_buf);
  obj->input_buf = NULL;
  obj->mod_buf = NULL;
  obj->alloc_length = 0;
  obj->line_stored = 0;

  free(obj);  
}
