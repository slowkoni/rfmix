/* (c) 2007-2016 Mark Hamilton Wright */
#ifndef _CMDLINE_UTILS_H
#define _CMDLINE_UTILS_H

#ifndef UCHAR
typedef unsigned char uchar;
#define UCHAR
#endif

/* Use these more brief type indicators in new code -- older ones below are
   annoying and take up more line space */
enum { OPT_STR = 0, OPT_INT, OPT_DBL, OPT_FLAG, OPT_CHAR, N_OPT };

/* Older, longer and annoying option types enumeration -- for compatibility
   with a large amount of existing code using them */
enum { OPTION_STRINGTYPE = 0, OPTION_INTTYPE, OPTION_DOUBLETYPE, OPTION_FLAGTYPE, N_OPTION_TYPE };

typedef struct {
  char shortopt;
  char const *longopt;
  void *var;
  int type;
  int required;
  int argument_required;
  char const *docstring;
} option_t;

#ifdef __cplusplus
extern "C" 
#endif

int cmdline_getoptions(option_t *options, int argc, char *argv[]);
void cmdline_printoptions(option_t *options);

#endif
