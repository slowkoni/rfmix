/* (c) 2007-2016 Mark Hamilton Wright */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>

#include <unistd.h>

#include "cmdline-utils.h"

static char const* typestring[N_OPTION_TYPE] = { "string", "int", "float", "flag" };

static void setoption(option_t *option, char *optarg) {

  switch(option->type) {
  case OPTION_STRINGTYPE:
    *((char **) option->var) = strdup(optarg);
    break;
  case OPTION_INTTYPE:
    *((int *) option->var) = atoi(optarg);
    break;
  case OPTION_DOUBLETYPE:
    *((double *) option->var) = atof(optarg);
    break;
  case OPTION_FLAGTYPE:
    *((int *) option->var) ^= 0x1;
    break;
  case OPT_CHAR:
    *((char *) option->var) = optarg[0];
    break;
  }
}

int cmdline_getoptions(option_t *options, int argc, char *argv[]) {
  int max_optsize;
  char *optname, *optarg;
  int longopt;
  int i,j,k;
  int error;

  if (argc == 1 || strcmp(argv[1],"--help") == 0) {
    cmdline_printoptions(options);
    exit(0);
  }

  max_optsize = 0;
  for(i=1;i<argc;i++) {
    if (strlen(argv[i]) > max_optsize) {
      max_optsize = strlen(argv[i]);
    }
  }
  optname = alloca(max_optsize+1);


  error = 0;
  i = 1;
  while(i<argc) {
    if (argv[i][0] == '-') {
      if (argv[i][1] == '-') {
	longopt = 1;
	k = 2;
	while(argv[i][k] != 0 && argv[i][k]!='=') {
	  optname[k-2] = argv[i][k];
	  k++;
	}
	optname[k-2] = 0;
	if (argv[i][k] == 0) {
	  optarg = NULL;
	} else {
	  optarg = argv[i] + k + 1;
	}
      } else {
	longopt = 0;
	optname[0] = argv[i][1];
	optarg = argv[i+1];
      }
    } else {
      /* Not an option then... */
      i++;
      continue;
    }

    j = 0;
    while(options[j].var != NULL) {
      if (longopt) {
	if (strcmp(optname, options[j].longopt)==0) {
	  setoption(options + j, optarg);
	  break;
	}	
      } else {
	if (options[j].shortopt == optname[0]) {
	  setoption(options + j, optarg);
	  break;
	}
      }
      j++;
    }
    if (options[j].var == NULL) {
      fprintf(stderr,"Unrecognized option \"%s\"\n",argv[i]);
      error = 1;
    }
    if (longopt || options[j].type == OPT_FLAG) {
      i++;
    } else {
      i+=2;
    }
  }

  return error;
}


void cmdline_printoptions(option_t *options) {
  FILE *f;
  int i;

  f = isatty(1) ? stdout : stderr;

  fprintf(f,"Summary of command line options - see manual for details\n\n");
  i = 0;
  while(options[i].var != NULL) {  
    fprintf(f,"   ");
    if (options[i].shortopt) {
      fprintf(f,"-%c", options[i].shortopt);
      if (options[i].argument_required)
	fprintf(f," <%s>", typestring[options[i].type]);

      if (options[i].longopt) fprintf(f,", ");
    }
    
    if (options[i].longopt) {
      fprintf(f,"--%s", options[i].longopt);
      if (options[i].argument_required)
	fprintf(f,"=<%s>", typestring[options[i].type]);
    }
    if (options[i].required)
      fprintf(f,"\t(required)");
    fprintf(f,"\n");
    
    if (options[i].docstring)
      fprintf(f,"      %s\n", options[i].docstring);

    i++;
  }
  fprintf(f,"\n\n");
}
