/* (c) 2007-2016 Mark Hamilton Wright */
#ifndef KMACROS_H
#define KMACROS_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <stdint.h>

#include <sys/time.h>

#ifndef MA
#define MA(p,s,t)       {			\
    (p) = (t*) malloc(sizeof(t)*(s));		\
    if ((p) == NULL) {						   \
      fprintf(stderr,"Failed allocating memory at %s:%d (%1.1f Mb)\n",	\
	  __FILE__,__LINE__, (s)/1e6);				    \
      abort();							    \
    }								    \
}

#define CA(p,s)	 {			\
    (p) = malloc(s);			\
if ((p) == NULL) {  \
  fprintf(stderr,"Failed allocating memory at %s:%d (%1.1f Mb)\n", \
	  __FILE__,__LINE__, (s)/1e6);				    \
  abort();							    \
  }\
    memset((p), 0x0, (s)); \
}

#define RA(p,s,t)           {			\
    (p) = (t *) realloc((p),(s)*sizeof(t));	\
if ((p) == NULL) {  \
  fprintf(stderr,"Failed allocating memory at %s:%d (%1.1f Mb)\n", \
	  __FILE__,__LINE__, (s)/1e6);				    \
  abort();							    \
 }                                \
}
#endif

#define CHOMP(s) { int l; char *s1; s1 = (s); l = strlen(s1)-1; while(l>=0 && (s1[l]=='\n' || s1[l]=='\r')) { s1[l--] = 0; } }

#define STREQ(s1, s2) (strcmp((s1), (s2))==0)
#define STRNEQ(s1, s2) (strncmp((s1), (s2))==0)

#define FREAD(p, s, n, f) if (fread((p), (s), (n), (f)) < (n)) { fprintf(stderr,"Warning: fread() returned short read at %s:%d\n",__FILE__,__LINE__); }

#define FWRITE(p, s, n, f) if (fwrite((p), (s), (n), (f)) < (n)) { fprintf(stderr,"Warning: fwrite() returned short write at %s:%d\n",__FILE__,__LINE__); }


#define WHEREARGS const char *callfunc, const char *callfile, int callline
#define WHEREFROM __func__, __FILE__, __LINE__

#ifndef FOPEN
#define FOPEN
#define FOPEN_EXIT_ON_ERROR     (1)
#define FOPEN_CONTINUE_ON_ERROR (0)
static FILE *k_fopen(char *fname, char *mode, int exit_on_error) __attribute__((unused));
static FILE *k_fopen(char *fname, char *mode, int exit_on_error) {
  FILE *f;

  if (fname == NULL || fname[0] == 0) return NULL;
 
  f = fopen(fname, mode); 
  if (f == NULL) { 
    fprintf(stderr,"Can't open file %s for %s (%s)\n", fname, 
	    (mode[0] == 'r') ? "input" : "output", strerror(errno)); 
    if (exit_on_error == FOPEN_EXIT_ON_ERROR) exit(-1);
  } 

  return f;
}
#endif

#if 0
typedef union {
  char ch;
  char *str;
  int32_t int32;
  int64_t int64;
  uint32_t uint32;
  uint64_t uint64;
  double dbl;
} splitfield_type_t;

#define SPLITFIELD_INITIAL_ALLOC (16)
typedef struct {
  int n_fields;
  char **fields;
  splitfield_type_t *tf;

  char *s_copy;
  int n_alloc;
} splitfield_t;

#define SPLIT_NOLIMIT (-1)
enum { SPLIT_INLINE = 0, SPLIT_MAKECOPY };
static splitfield_t *split(splitfield_t *sf, char *s, char *d, char *types,
			   int limit, int use_copy) __attribute__((unused));

/* NOTE: if use_copy is set, then input string (s) is not modified and fields
   extracted point to the copy string, so the input can be free()'d or 
   altered. However, the field strings themselves are not copies, so if they
   are free()'d a segfault will probably happen immediately or soon, and
   free()ing the first field will free the copy that all other fields
   are pointing to */
static splitfield_t *split(splitfield_t *sf, char *s, char *d, char *types,
			   int limit, int use_copy) {
  char *p, *q;
  int n_types;

  if (sf == NULL) {
    MA(sf, sizeof(splitfield_t));
    sf->s_copy = NULL;

    MA(sf->fields, sizeof(char *)*SPLITFIELD_INITIAL_ALLOC);
    MA(sf->tf, sizeof(splitfield_type_t)*SPLITFIELD_INITIAL_ALLOC);
    sf->n_alloc = SPLITFIELD_INITIAL_ALLOC;
  }

  sf->n_fields = 0;
  if (sf->s_copy) {
    free(sf->s_copy);
    sf->s_copy = NULL;
  }

  if (use_copy) sf->s_copy = s = strdup(s);
  n_types = (types == NULL) ? 0 : strlen(types);

  p = s;
  while((q = strsep(&p, d)) != NULL && (limit == -1 || sf->n_fields < limit)) {
    if (sf->n_fields >= sf->n_alloc) {
      sf->n_alloc <<= 1;
      RA(sf->fields, sizeof(char *)*sf->n_alloc);
      RA(sf->tf, sizeof(splitfield_type_t)*sf->n_alloc);
    }

    sf->fields[sf->n_fields] = q;
    if (types && sf->n_fields < n_types) {
      switch(types[sf->n_fields]) {
      case 'c':
	sf->tf[sf->n_fields].ch = q[0];
	break;
      case 'd':
	sf->tf[sf->n_fields].int32 = strtol(q, NULL, 10);
	break;
      case 'u':
	sf->tf[sf->n_fields].uint32 = strtoul(q, NULL, 10);
	break;
      case 'D':
	sf->tf[sf->n_fields].int64 = strtol(q, NULL, 10);
	break;
      case 'U':
	sf->tf[sf->n_fields].uint64 = strtoul(q, NULL, 10);
	break;
      case 's':
	sf->tf[sf->n_fields].str = q;
	break;
      case 'f':
	sf->tf[sf->n_fields].dbl = atof(q);
	break;
      case ' ':
	/* don't convert field */
	break;
      default:
	fprintf(stderr,"split(): unknown conversion type %c\n", 
		types[sf->n_fields]);
	break;
      }
    }

    sf->n_fields++;
  }

  return sf;
}

static splitfield_t *sf_init(void) __attribute__((unused));
static splitfield_t *sf_init(void) {
  splitfield_t *sf;

  MA(sf, sizeof(splitfield_t));
  sf->s_copy = NULL;
  
  MA(sf->fields, sizeof(char *)*SPLITFIELD_INITIAL_ALLOC);
  MA(sf->tf, sizeof(splitfield_type_t)*SPLITFIELD_INITIAL_ALLOC);
  sf->n_alloc = SPLITFIELD_INITIAL_ALLOC;
  
  return sf;
}

static void sf_destroy(splitfield_t *sf) __attribute__((unused));
static void sf_destroy(splitfield_t *sf) {
  if (sf->s_copy) free(sf->s_copy);
  free(sf->fields);
  free(sf->tf);
  free(sf);	    
}
#endif
 
static double __attribute__((unused))elapsed_time_ms(struct timeval *ref_time) {
  struct timeval now;

  gettimeofday(&now, NULL);
  return (now.tv_sec - ref_time->tv_sec) * 1000. +
    (now.tv_usec - ref_time->tv_usec)/1000.;
}

static double __attribute__((unused))elapsed_time(struct timeval *ref_time) {
  struct timeval now;

  gettimeofday(&now, NULL);
  return (now.tv_sec - ref_time->tv_sec) +
    (now.tv_usec - ref_time->tv_usec)/1000000.;
}

#define SEARCH_STRS_NULL_TERM (-1)
#define SEARCH_STRS_NOT_FOUND (-1)
static int __attribute__((unused))search_strs(char **strs, char *q, int n) {
  int i;

  if (strs == NULL) return SEARCH_STRS_NOT_FOUND;

  i = 0;
  while(i<n || (n == SEARCH_STRS_NULL_TERM && strs[i] != NULL)) {
    if (strs[i] != NULL && strcmp(strs[i], q) == 0) return i;
    i++;
  }

  return SEARCH_STRS_NOT_FOUND;
}

static void __attribute__((unused))join(FILE *f, char *d, char **strs, int n) {
  int i;

  fprintf(f,"%s", strs[0]);
  for(i=1;i<n;i++) {
    fprintf(f,"%s%s", d, strs[i]);
  }
}

#define ADD_ITEM(list, item, n_items, alloc_step, type) {	\
    type x;							\
    int as;                                                     \
    int n;                                                      \
    x = (item);							\
    n = (n_items);                                              \
    as = (alloc_step);                                          \
    if (n_items % as == 0) {					\
      RA((list), sizeof(type)*(n + alloc_step));		\
    };								\
    (list)[n] = x;                                              \
    n_items = n+1;						\
  }


#endif
