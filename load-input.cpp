#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "kmacros.h"
#include "rfmix.h"
#include "genetic-map.h"
#include "load-input.h"
#include "inputline.h"
#include "hash-table.h"

extern rfmix_opts_t rfmix_opts;

#define SUBPOP_ALLOC_STEP (8)
#define SAMPLE_ALLOC_STEP (256)
#define SNP_ALLOC_STEP (16384)

/* Reads through and ignores all VCF header lines and returns the sample header line */
static char *vcf_skip_headers(Inputline *vcf) {
  char *p;

  while((p = vcf->nextline(INPUTLINE_NOCOPY)) != NULL) {
    CHOMP(p);
    if (strncmp(p, "#CHROM", 6) == 0) break;
    if (p[0] == 0 || p[0] == '#') continue;
  }
  if (p == NULL) {
    fprintf(stderr,"\nNo genotype data found in VCF file %s.\n\n", vcf->fname);
    exit(-1);
  }
  return p;
}

static void load_samples(input_t *input) {
  sample_t *samples;
  int n_samples, i, ref_idx, *tmp;
  char *sample_id, *reference_pop, *p;

  MA(input->reference_subpops, sizeof(char *)*SUBPOP_ALLOC_STEP, char *);
  input->reference_subpops[0] = strdup("query");
  input->n_subpops = 1;
  
  /* This will be used to rapidly locate a sample when matching the reference VCF file 
     to the sample ids loaded from the sample map file. Any sample_id found in the VCF
     sample header line that is not defined in the sample map file will be excluded.
     The need for the hash table is perceived in the case the reference file has a very
     large number of samples (thousands). The integer index to the samples array is
     the data value stored in the hash table, because the samples array may be
     copied to a new location when extended by RA() (realloc). */
  HashTable *sample_hash = new HashTable(256);
  samples = NULL;
  n_samples = 0;
  
  /* All samples in the query VCF file will be analyzed, and are not expected to be
     named in the seperate sample map file. Grab them from the VCF header and add
     them to the sample array first */
  Inputline *qvcf = new Inputline(rfmix_opts.qvcf_fname);
  p = vcf_skip_headers(qvcf);

  CHOMP(p);
  for(i=0; i < 9; i++) strsep(&p, "\t");
  while((sample_id = strsep(&p, "\t")) != NULL) {
    if (n_samples % SAMPLE_ALLOC_STEP == 0)
      RA(samples, sizeof(sample_t)*(SAMPLE_ALLOC_STEP + n_samples), sample_t);
	
    samples[n_samples].sample_id = strdup(sample_id);
    samples[n_samples].apriori_subpop = 0;

    MA(tmp, sizeof(int), int);
    *tmp = n_samples;
    sample_hash->insert(sample_id, tmp);
    n_samples++;
  }
  delete qvcf; 

  /* Now load sample ids from the sample map file */
  Inputline *f = new Inputline(rfmix_opts.class_fname);

  /* Now scan the sample map file and determine the reference subpops and sample mapping to them */
  while((p = f->nextline(INPUTLINE_NOCOPY)) != NULL) {
    CHOMP(p);
    if (p[0] == 0 || p[0] == '#') continue;

    sample_id = strsep(&p, "\t");
    reference_pop = strsep(&p, "\t");

    /* Search for this reference subpop in the already known list */
    for(i=0; i < input->n_subpops; i++)
      if (strcasecmp(input->reference_subpops[i], reference_pop) == 0) break;
      
    /* Add the subpop name to the list of reference subpops if it is not found */
    if (i == input->n_subpops) {
      if (input->n_subpops % SUBPOP_ALLOC_STEP == 0)
	RA(input->reference_subpops, sizeof(char *)*(input->n_subpops + SUBPOP_ALLOC_STEP), char *);
      input->reference_subpops[input->n_subpops] = strdup(reference_pop);
      input->n_subpops++;
    }
    /* Whether just added or found, the reference subpop index is i */
    ref_idx = i;

    if (n_samples % SAMPLE_ALLOC_STEP == 0)
      RA(samples, sizeof(sample_t)*(SAMPLE_ALLOC_STEP + n_samples), sample_t);
	
    samples[n_samples].sample_id = strdup(sample_id);
    samples[n_samples].apriori_subpop = ref_idx;

    MA(tmp, sizeof(int), int);
    *tmp = n_samples;
    sample_hash->insert(sample_id, tmp);
    n_samples++;
  }
  delete f;
  
  /* initialize to empty/null values all other sample struct fields */
  for(i=0; i < n_samples; i++) {
    samples[i].haplotype[0] = NULL;
    samples[i].haplotype[1] = NULL;
    samples[i].current_p[0] = NULL;
    samples[i].current_p[1] = NULL;
    samples[i].est_p[0] = NULL;
    samples[i].est_p[1] = NULL;
    samples[i].est_p[2] = NULL;
    samples[i].est_p[3] = NULL;
  }

  /* All work of this function is stored into the input_t struct and made available
     essentially everywhere through that */
  input->samples = samples;
  input->n_samples = n_samples;
  input->sample_hash = sample_hash;  
}

static void skip_to_chromosome(Inputline *vcf, char *chm) {
  char *p, *q;
  
  while((p = vcf->nextline(INPUTLINE_NOCOPY)) != NULL) {
    q = strsep(&p, "\t");
    if (strcmp(q, chm) == 0) {
      vcf->pushback();
      return;
    }
  }

  fprintf(stderr,"\nCan't find chromosome %s in VCF file %s\n\n", chm, vcf->fname);
  exit(-1);
}

/* Grabs the next SNP line from a VCF and returns the chromosome, position, and the snp line
   pointing to the next field past the position */
static char *get_next_snp(Inputline *vcf, char **chm, int *pos) {
  char *p = vcf->nextline();
  if (p == NULL) {
    *chm = NULL;
    *pos = -1;
    return NULL;
  }
  CHOMP(p);
  
  *chm = strsep(&p, "\t");
  char *q = strsep(&p, "\t");

  *pos = atoi(q);
  return p;
}

static void identify_common_snps(input_t *input) {
  char *qvcf_sample_header, *rvcf_sample_header;

  Inputline *qvcf = new Inputline(rfmix_opts.qvcf_fname);
  vcf_skip_headers(qvcf);
  skip_to_chromosome(qvcf, rfmix_opts.chromosome);
  
  Inputline *rvcf = new Inputline(rfmix_opts.rvcf_fname);
  vcf_skip_headers(rvcf);
  skip_to_chromosome(rvcf, rfmix_opts.chromosome);
  
  snp_t *snps = NULL;
  int n_snps = 0;
  char *q_chm, *r_chm;
  int q_pos, r_pos;
  
  get_next_snp(qvcf, &q_chm, &q_pos);
  get_next_snp(rvcf, &r_chm, &r_pos);

  for(;;) {
    while(q_pos != -1 && strcmp(q_chm, rfmix_opts.chromosome) == 0 &&
	  q_pos < r_pos)
      get_next_snp(qvcf, &q_chm, &q_pos);
    if (q_pos == -1 || strcmp(q_chm, rfmix_opts.chromosome) != 0) break;

    while(r_pos != -1 && strcmp(r_chm, rfmix_opts.chromosome) == 0 &&
	  r_pos < q_pos)
      get_next_snp(rvcf, &r_chm, &r_pos);
    if (r_pos == -1 || strcmp(r_chm, rfmix_opts.chromosome) != 0) break;

    if (q_pos == r_pos) {
      if (n_snps % SNP_ALLOC_STEP == 0)
	RA(snps, sizeof(snp_t)*(n_snps + SNP_ALLOC_STEP), snp_t);
      snps[n_snps].pos = q_pos;
      snps[n_snps].genetic_pos = input->genetic_map->translate_seqpos(q_pos);
      snps[n_snps].crf_index = -1;
      n_snps++;

      get_next_snp(qvcf, &q_chm, &q_pos);      
      get_next_snp(rvcf, &r_chm, &r_pos);
    }
  }

  input->snps = snps;
  input->n_snps = n_snps;

  delete qvcf;
  delete rvcf;
}

#define VCF_LEAD_COLS (9)
typedef struct {
  int col;
  char *sample_id; // NOTE: not a sample_id if leading VCF cols
  int sample_idx;
} vcf_column_map_t;

static int vcf_parse_column_header(vcf_column_map_t **rcolumn_map, char *column_header,
				   input_t *input) {
  vcf_column_map_t *column_map = NULL;
  int n_cols = 0;
  char *p = column_header;
  char *q;
  while((q = strsep(&p,"\t")) != NULL) {
    if (n_cols % SAMPLE_ALLOC_STEP == 0)
      RA(column_map, sizeof(vcf_column_map_t)*(n_cols + SAMPLE_ALLOC_STEP), vcf_column_map_t);
    
    column_map[n_cols].col = n_cols;
    column_map[n_cols].sample_id = strdup(q);
    int *tmp = (int *) input->sample_hash->lookup(q);
    if (tmp != NULL) {
      column_map[n_cols].sample_idx = *tmp;
    } else {
      column_map[n_cols].sample_idx = -1;
    }
    n_cols++;
  }

  *rcolumn_map = column_map;
  return n_cols;
}

static int8_t get_allele(char q) {
  switch(q) {
  case '0':
    return 0;
  case '1':
    return 1;
  case '.':
  default:
    return 2;
  }
}

static void parse_alleles(input_t *input, Inputline *vcf, vcf_column_map_t *column_map,
			  int n_cols) {
  char *p, *q;
  char *chm;
  int pos;

  int snp_idx = 0;
  while(snp_idx < input->n_snps &&
	(p = get_next_snp(vcf, &chm, &pos)) != NULL &&
	strcmp(chm, rfmix_opts.chromosome) == 0) {
    if (input->snps[snp_idx].pos != pos) continue;

    int col_idx = 2;
    while(col_idx < VCF_LEAD_COLS) { strsep(&p, "\t"); col_idx++; }

    while(col_idx < n_cols && (q = strsep(&p, "\t")) != NULL) {
      if (column_map[col_idx].sample_idx == -1) {
	col_idx++;
	continue;
      }
      sample_t *sample = input->samples + column_map[col_idx].sample_idx;
      
      if (strlen(q) < 2) {
	fprintf(stderr,"VCF parsing error - valid genotype not detected on line %d of %s\n",
		vcf->line_no, vcf->fname);
	exit(-1);
      }
      if (q[1] != '|' && q[0] != '.' && q[2] != '.') {
	fprintf(stderr,"Warning: unphased genotype detected on line %d of %s\n",
		vcf->line_no, vcf->fname);
      }
      sample->haplotype[0][snp_idx] = get_allele(q[0]);
      sample->haplotype[1][snp_idx] = get_allele(q[2]);
      col_idx++;
    }
    snp_idx++;
  }
}

static void load_alleles(input_t *input) {
  Inputline *qvcf = new Inputline(rfmix_opts.qvcf_fname);
  char *sample_header = vcf_skip_headers(qvcf);

  vcf_column_map_t *column_map;
  int n_cols = vcf_parse_column_header(&column_map, sample_header, input);

  skip_to_chromosome(qvcf, rfmix_opts.chromosome);
  parse_alleles(input, qvcf, column_map, n_cols);

  delete qvcf;
  for(int i=0; i < n_cols; i++) {
    if (column_map[i].sample_id) free(column_map[i].sample_id);
  }
  free(column_map);
  n_cols = 0;

  Inputline *rvcf = new Inputline(rfmix_opts.rvcf_fname);
  sample_header = vcf_skip_headers(rvcf);
  n_cols = vcf_parse_column_header(&column_map, sample_header, input);

  skip_to_chromosome(rvcf, rfmix_opts.chromosome);
  parse_alleles(input, rvcf, column_map, n_cols);

  delete rvcf;
  for(int i=0; i < n_cols; i++) {
    if (column_map[i].sample_id) free(column_map[i].sample_id);
  }
  free(column_map);
  n_cols = 0;
}

static void set_crf_points(input_t *input) {
  /* Determine the number of defined CRF points we have (CRF windows), the
     central SNP that defines each one, and the boundaries of the larger
     window used to source SNPs for the random forest classification */
  fprintf(stderr,"\n\tsetting up CRF points and random forest windows... ");
  input->n_windows = input->n_snps / rfmix_opts.crf_spacing;
  MA(input->crf_windows, sizeof(crf_window_t)*input->n_windows, crf_window_t);

  int rf_start, rf_end;
  int j = rfmix_opts.crf_spacing / 2;
  for(int i=0; i < input->n_windows; i++,j += rfmix_opts.crf_spacing) {
    input->crf_windows[i].snp_idx = j;

    rf_start = j - rfmix_opts.rf_window_size / 2;
    if (rf_start < 0) rf_start = 0;

    rf_end = rf_start + rfmix_opts.rf_window_size - 1;
    if (rf_end >= input->n_snps) {
      rf_end = input->n_snps - 1;
      rf_start = rf_end - rfmix_opts.rf_window_size + 1;
      if (rf_start < 0) rf_start = 0;
    }

    input->crf_windows[i].rf_start_idx = rf_start;
    input->crf_windows[i].rf_end_idx = rf_end;
  }
  fprintf(stderr,"done\n");
  
  /* Set up and initialize the current (starting) marginal probabilities for subpop
     assignment for each haplotype at each CRF window. These values start as 100%
     probability the haplotypes are from the apriori subpopulation for reference
     individuals, and just initialized to zero for all query individuals. These are
     calculated at each EM iteration by the Forward-Backward algorithm in the
     conditional random field code */
  fprintf(stderr,"\tinitializing apriori reference subpop across CRF... ");
  for(int k=0; k < input->n_samples; k++) {
    sample_t *sample = input->samples + k;

    for(int h=0; h < 2; h++) {
      MA(sample->current_p[h], sizeof(AF_TYPE *)*input->n_windows, AF_TYPE *);

      for(int i=0; i < input->n_windows; i++) {
	MA(sample->current_p[h][i], sizeof(AF_TYPE)*input->n_subpops, AF_TYPE);

	for(int s=0; s < input->n_subpops; s++) 
	  sample->current_p[h][i][s] = 0.;
	  
	if (sample->apriori_subpop != 0) sample->current_p[h][i][sample->apriori_subpop - 1] = 1.;
      }
    }
  }
  fprintf(stderr,"done\n");

  fprintf(stderr,"\tsetting up random forest probability estimation arrays... ");
  for(int k=0; k < input->n_samples; k++) {
    sample_t *sample = input->samples + k;

    for(int h=0; h < 4; h++) {
      MA(sample->est_p[h], sizeof(AF_TYPE *)*input->n_windows, AF_TYPE *);

      for(int i=0; i < input->n_windows; i++)
	MA(sample->est_p[h][i], sizeof(AF_TYPE)*input->n_subpops, AF_TYPE);
    }
      /*    MA(sample->est_p, sizeof(AF_TYPE *)*input->n_windows, AF_TYPE *);
    for(int i=0; i < input->n_windows; i++)
    MA(sample->est_p[i], sizeof(AF_TYPE)*n_states, AF_TYPE);*/
  }
  fprintf(stderr,"done\n");
}

input_t *load_input(void) {
  input_t *input;
  MA(input, sizeof(input_t), input_t);

  fprintf(stderr,"Loading genetic map for chromosome %s ...  ", rfmix_opts.chromosome);
  input->genetic_map = new GeneticMap();
  input->genetic_map->load_map(rfmix_opts.genetic_fname, rfmix_opts.chromosome);
  fprintf(stderr,"done\n");
  
  fprintf(stderr,"Scanning input VCFs for common SNPs on chromosome %s ...   ", rfmix_opts.chromosome);
  identify_common_snps(input);
  fprintf(stderr,"%d SNPs\n", input->n_snps);

  /* Find and map out all the samples that we will be loading */
  fprintf(stderr,"Mapping samples ... ");
  load_samples(input);
  fprintf(stderr,"%d samples combined\n", input->n_samples);
  
  /* Now we know all the samples that we will be loading, and all the SNPs,
     allocate the space to store the haplotypes */
  for(int i=0; i < input->n_samples; i++) {
    MA(input->samples[i].haplotype[0], sizeof(int8_t)*input->n_snps, int8_t);
    MA(input->samples[i].haplotype[1], sizeof(int8_t)*input->n_snps, int8_t);
  }

  fprintf(stderr,"Loading haplotypes... ");
  load_alleles(input);
  fprintf(stderr,"done\n");

  fprintf(stderr,"Defining and initializing conditional random field...  ");
  set_crf_points(input);
  fprintf(stderr,"done\n");
  
  return input;
}


void free_input(input_t *input) {

  free(input->snps);
  input->n_snps = 0;
  
  for(int i=0; i < input->n_samples; i++) {
    sample_t *sample = input->samples + i;

    for(int h = 0; h < 2; h++) {
      free(sample->haplotype[h]);
      for(int w=0; w < input->n_windows; w++)
	free(sample->current_p[h][w]);
      free(sample->current_p[h]);
    }

    for(int h=0; h < 4; h++) {
      for(int w=0; w < input->n_windows; w++)
	free(sample->est_p[h][w]);
      free(sample->est_p[h]);
    }
  }

  free(input->crf_windows);

  delete input->sample_hash;
  free(input->samples);
  input->n_samples = 0;
  
  delete input->genetic_map;
  free(input);
}
