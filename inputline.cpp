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

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>

#include "kmacros.h"
#include "inputline.h"
#include "rfmix.h"

#define INPUTLINE_CHUNK (8192)
FILE *Inputline::open_gzip_read(char *fname) {
  int fds[2];
  int child_pid;
  FILE *f;
  
  if (pipe(fds) != 0) {
    fprintf(stderr,"Can't create pipe to read output from gzip -dc (%s)\n", strerror(errno));
    exit(-1);
  }

  child_pid = fork();
  if (child_pid == 0) {
    /* New child process executes this code, which just dups the write end of the pipe to
       standard output, so that gzip which writes to stdout will really be outputing to 
       the pipe for the parent to read */
    close(0);
    close(1);
    close(fds[0]);
    dup2(fds[1], 1);
    close(fds[1]);
    execlp("gzip","gzip","-dc",fname, NULL);

    /* If gzip is successfully started, execlp() never returns and these lines
       are never reached. Only if execlp() fails will the process continue here */
    fprintf(stderr,"Can't execute gzip to decompress input file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  } else {
    /* parent process. If the child_pid is -1, then fork() failed, which pretty much
       never happens. */
    if (child_pid == -1) {
      fprintf(stderr,"Can't fork() a child process to start gzip for decompression (%s)\n",
	      strerror(errno));
      exit(-1);
    }

    /* Close the write end of the pipe because we are only going to read. Then open a normal
       stdio FILE stream from the pipe read end file descriptor. We set the child_pid field
       of the object so that the destructor knows to call wait() to reap the child process.
       Otherwise, it is a zombie hanging around */
    close(fds[1]);
    f = fdopen(fds[0], "r");
    if (f == NULL) {
      fprintf(stderr,"Can't open file handle from file descriptor to read gzip output (%s)\n",
	      strerror(errno));
      exit(-1);
    }
    this->child_pid = child_pid;
  }
  
  return f;
}

FILE *Inputline::open_bcftools_read(char *fname, char *chm) {
  int fds[2];
  int child_pid;
  FILE *f;
  
  if (pipe(fds) != 0) {
    fprintf(stderr,"Can't create pipe to read output from bcftools (%s)\n", strerror(errno));
    exit(-1);
  }

  child_pid = fork();
  if (child_pid == 0) {
    /* New child process executes this code, which just dups the write end of the pipe to
       standard output, so that bcftools which writes to stdout will really be outputing to 
       the pipe for the parent to read */
    close(0);
    close(1);
    close(fds[0]);
    dup2(fds[1], 1);
    close(fds[1]);
    if (chm != NULL) {
      execlp("bcftools","bcftools","view","--regions", chm, fname, NULL);
    } else {
      execlp("bcftools","bcftools","view", fname, NULL);
    }
    
    /* If bcftools is successfully started, execlp() never returns and these lines
       are never reached. Only if execlp() fails will the process continue here */
    fprintf(stderr,"Can't execute bcftools to read BCF file %s (%s) - is bcftools installed?\n",
	    fname, strerror(errno));
    exit(-1);
  } else {
    /* parent process. If the child_pid is -1, then fork() failed, which pretty much
       never happens. */
    if (child_pid == -1) {
      fprintf(stderr,"Can't fork() a child process to start bcftools (%s)\n",
	      strerror(errno));
      exit(-1);
    }

    /* Close the write end of the pipe because we are only going to read. Then open a normal
       stdio FILE stream from the pipe read end file descriptor. We set the child_pid field
       of the object so that the destructor knows to call wait() to reap the child process.
       Otherwise, it is a zombie hanging around */
    close(fds[1]);
    f = fdopen(fds[0], "r");
    if (f == NULL) {
      fprintf(stderr,"Can't open file handle from file descriptor to read bcftools output (%s)\n",
	      strerror(errno));
      exit(-1);
    }
    this->child_pid = child_pid;
  }
  
  return f;
}

Inputline::Inputline(char *fname, char *chm) {
  /* use the filename extension to determine if we need a helper program to read 
     the file. Functions above take care of starting the relevant program and
     piping its output back to a FILE handle we can read like any other file */
  int l = strlen(fname);
  if ((l > 7 && strcmp(".bcf.gz", fname + l - 7) == 0) ||
      (l > 4 && strcmp(".bcf", fname +l - 4) == 0)) {
    f = open_bcftools_read(fname, chm);
  }
  else if (l > 3 && strcmp(".gz", fname + l - 3) == 0) {
    f = open_gzip_read(fname);
  }
  else {
    f = fopen(fname, "r");
    child_pid = 0; // There is no helper program in a child process if we are reading the file directly
  }
  
  if (f == NULL) {
    fprintf(stderr,"\nCan't open input file %s (%s)\n\n", fname, strerror(errno));
    exit(-1);
  }
  this->fname = strdup(fname);

  MA(input_buf, INPUTLINE_CHUNK, char);
  MA(mod_buf, INPUTLINE_CHUNK, char);
  alloc_length = INPUTLINE_CHUNK;
  line_stored = 0;
  line_no = 0;

  pthread_mutex_init(&lock, NULL);
}

Inputline::Inputline(char *fname) {
  Inputline(fname, NULL);
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

char *Inputline::nextline(void) {
  return nextline(INPUTLINE_NOCOPY);
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
  free(fname);
  if (child_pid) {
    kill(child_pid, SIGINT);
    waitpid(child_pid, NULL, 0);
  }
  
  free(input_buf);
  input_buf = NULL;

  free(mod_buf);
  mod_buf = NULL;

  alloc_length = 0;
  line_stored = 0;
}
