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

#include <float.h>

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

#define VCF_LEAD_COLS (9)

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

static int8_t get_allele(char q) {
  switch(q) {
  case '0':
    return 0;
  case '1':
    return 1;
  case '.':
    return 2;
  default:
    return 2;
  }
}

static void load_samples(input_t *input) {
  sample_t *samples;
  int n_samples, i, ref_idx, *tmp;
  char *sample_id, *reference_pop, *p;

  input->reference_subpops = NULL;
  input->n_subpops = 0;
  
  /* This will be used to rapidly locate a sample when matching the reference VCF file 
     to the sample ids loaded from the sample map file. Any sample_id found in the VCF
     sample header line that is not defined in the sample map file will be excluded.
     The need for the hash table is perceived in the case the reference file has a very
     large number of samples (thousands). The integer index to the samples array is
     the data value stored in the hash table, because the samples array may be
     copied to a new location when extended by RA() (realloc). */
  HashTable *sample_hash = new HashTable(256);
  HashTable *tmp_hash = new HashTable(256); // used to check if sample defined in sample map is in reference VCF
  samples = NULL;
  n_samples = 0;
  
  /* All samples in the query VCF file will be analyzed, and are not expected to be
     named in the seperate sample map file. Grab them from the VCF header and add
     them to the sample array first */
  Inputline *qvcf = new Inputline(rfmix_opts.qvcf_fname, rfmix_opts.chromosome);
  p = vcf_skip_headers(qvcf);

  CHOMP(p);
  for(i=0; i < 9; i++) strsep(&p, "\t");
  while((sample_id = strsep(&p, "\t")) != NULL) {
    if (n_samples % SAMPLE_ALLOC_STEP == 0)
      RA(samples, sizeof(sample_t)*(SAMPLE_ALLOC_STEP + n_samples), sample_t);
	
    samples[n_samples].sample_id = strdup(sample_id);
    samples[n_samples].apriori_subpop = -1;
    samples[n_samples].s_parent = 0;
    samples[n_samples].s_sample = 0;

    if (sample_hash->lookup(sample_id) != NULL) {
      fprintf(stderr,"Error: Sample id %s occurs twice or more in input - samples must have unique identifiers both within and across query and reference\n", sample_id);
      exit(-1);
    }
    
    MA(tmp, sizeof(int), int);
    *tmp = n_samples;
    sample_hash->insert(sample_id, tmp);
    
    n_samples++;
  }
  delete qvcf; 

  /* Parse up the column header from the reference VCF to see what reference
     samples are actually in the file. Then, we will only define a sample in
     the program's sample array if in fact that sample is present in the
     reference. Otherwise, extraneous samples defined in the sample map cause
     havoc if they are never actually loaded from the reference. This does happen
     if a common sample map is used for different reference VCFs where certain
     subpopulations have simply not been included. */
  Inputline *rvcf = new Inputline(rfmix_opts.rvcf_fname, rfmix_opts.chromosome);
  p = vcf_skip_headers(rvcf);

  CHOMP(p);
  for(i=0; i < 9; i++) strsep(&p, "\t");
  while((sample_id = strsep(&p, "\t")) != NULL) {
    tmp_hash->insert(sample_id, sample_id);
  }

  delete rvcf;

  /* This is an embarassing afterthought hack - this insists the subpopulation name
     to index number goes in alphabetical order, always, for consistency even if the
     sample map is reordered but otherwise the same. What is more likely to happen
     is someone comments out a few samples causing the reference subpop names to
     first appear to the program in a different order, resulting in different
     indexing. To add this, I'm just scanning the entire sample map file first and
     only recording the subpop names, sorting them, then letting the original loop
     below scan the sample map for sample names and their subpops, find the subpop
     is already defined, and use the array index of the sorted array. Just hold your
     noses... this whole damn file needs a rewrite, not just this ugly hack */
  Inputline *f = new Inputline(rfmix_opts.class_fname, rfmix_opts.chromosome);

  while((p = f->nextline(INPUTLINE_NOCOPY)) != NULL) {
    
    CHOMP(p);
    if (p[0] == 0 || p[0] == '#' || p[0]=='^') continue;

    sample_id = strsep(&p, "\t");
    if (sample_id[0] == 0) continue;
    /* Do not define a sample entry if the sample is not in the reference vcf */
    if (tmp_hash->lookup(sample_id) == NULL) continue;
    
    reference_pop = strsep(&p, "\t");
    if (reference_pop[0] == 0) continue;

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
  }
  delete f;

  for(int i=0; i < input->n_subpops; i++) {
    int m = i;
    for(int j=i+1; j < input->n_subpops; j++)
      if (strcmp(input->reference_subpops[j], input->reference_subpops[m]) < 0) m = j;
    if (m != i) {
      char *tmp;
      tmp = input->reference_subpops[i];
      input->reference_subpops[i] = input->reference_subpops[m];
      input->reference_subpops[m] = tmp;
    }
  }
  
  /* Reopen the file and start over... sorry. */
  f = new Inputline(rfmix_opts.class_fname, rfmix_opts.chromosome);

  /* Now scan the sample map file and determine the reference subpops and sample mapping to them */
  while((p = f->nextline(INPUTLINE_NOCOPY)) != NULL) {
    CHOMP(p);
    if (p[0] == 0 || p[0] == '#' || p[0] == '^') continue;

    sample_id = strsep(&p, "\t");
    if (sample_id[0] == 0) continue;
    /* Do not define a sample entry if the sample is not in the reference vcf */
    if (tmp_hash->lookup(sample_id) == NULL) continue;
    
    reference_pop = strsep(&p, "\t");
    if (reference_pop[0] == 0) continue;
    
    /* Search for this reference subpop in the already known list */
    for(i=0; i < input->n_subpops; i++)
      if (strcasecmp(input->reference_subpops[i], reference_pop) == 0) break;
      
    /* Add the subpop name to the list of reference subpops if it is not found */
    if (i == input->n_subpops) {
      fprintf(stderr,"Warning: subpopulation %s should already be defined but it is not, "
	      "subpopulation order\nin output not guaranteed to be alphabetical.\n", reference_pop);
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

    if (sample_hash->lookup(sample_id) != NULL) {
      fprintf(stderr,"Error: Sample id %s occurs twice or more in input - samples must have unique identifiers both within and across query and reference\n", sample_id);
      exit(-1);
    }
    MA(tmp, sizeof(int), int);
    *tmp = n_samples;
    sample_hash->insert(sample_id, tmp);
    n_samples++;
  }
  delete f;
  
  /* initialize to empty/null values all other sample struct fields */
  for(i=0; i < n_samples; i++) {
    samples[i].s_sample = 0;
    samples[i].s_parent = 0;
    samples[i].haplotype[0] = NULL;
    samples[i].haplotype[1] = NULL;
    samples[i].current_p[0] = NULL;
    samples[i].current_p[1] = NULL;
    samples[i].est_p[0] = NULL;
    samples[i].est_p[1] = NULL;
    samples[i].est_p[2] = NULL;
    samples[i].est_p[3] = NULL;
    samples[i].sis_p[0] = NULL;
    samples[i].sis_p[1] = NULL;
    samples[i].logl[0] = -DBL_MAX;
    samples[i].logl[1] = -DBL_MAX;
    samples[i].logl[2] = -DBL_MAX;
    samples[i].logl[3] = -DBL_MAX;
    samples[i].msp[0] = NULL;
    samples[i].msp[1] = NULL;
    samples[i].msp[2] = NULL;
    samples[i].msp[3] = NULL;
    samples[i].ksp[0] = NULL;
    samples[i].ksp[1] = NULL;
  }
  
  /* All work of this function is stored into the input_t struct and made available
     essentially everywhere through that */
  input->samples = samples;
  input->n_samples = n_samples;
  input->sample_hash = sample_hash;

  delete tmp_hash;
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

/* Parses up a VCF genotype line and computes the minor allele frequency as well as 
   the proportion of missing data. This will be used to determine whether a SNP 
   should be included in the input.

   Afterthoughts: This functionality was added out of necessity after discovering that
                  large amounts of missing data in the input at SNPs caused serious
		  ambiguity in the estimation of ancestry at that CRF window. The
		  function duplicates code with parse_alleles() and a refactoring 
		  should consider combining the functions, loading the alleles
		  into a temporary array with this function just counting the alleles
		  in the temporary array, and parse_alleles() copying the alleles
		  into the haplotypes arrays, both calling the same function that
		  actually parses up the genotypes. */
static double vcf_snp_maf(int *r_mac, double *r_miss, char *p) {
  char *q;
  int n_total = 0;
  int n_obs = 0;
  double ref_freq = 0.;
  
  while((q = strsep(&p, "\t")) != NULL) {

    /* We are going to skip errors at this stage, they will be reported later when
       loading alleles. Errors will count as missing data since alleles can not be
       parsed from the genotype. If missing data is high the snp may simply be 
       excluded in which case we need not care about the errors later */
    if (strlen(q) < 2) {
      n_total += 2;
      continue;
    }

    /* NOTE: perhaps unphased genotypes should be treated as missing data and we
             should detect and count them as such here. Currently the allele 
	     loading code just loads the alleles if they were phased but prints
	     a warning */
    if (get_allele(q[0]) == 0) {
      ref_freq += 1.0;
      n_obs++;
    } else if (get_allele(q[0]) == 1) {
      n_obs++;
    }
    n_total++;

    if (get_allele(q[2]) == 0) {
      ref_freq += 1.0;
      n_obs++;
    } else if (get_allele(q[2]) == 1) {
      n_obs++;
    }
    n_total++;

  }
  
  *r_mac = ref_freq < n_obs - ref_freq ? ref_freq : n_obs - ref_freq;
  ref_freq /= n_obs;
  *r_miss = (n_total - n_obs) / (double) n_total;

  return *r_mac /(double) n_obs;
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

  int col_idx = 2;
  while(col_idx < VCF_LEAD_COLS) { strsep(&p, "\t"); col_idx++; }

  return p;
}

static void identify_common_snps(input_t *input) {
  char *pq, *pr;
  
  Inputline *qvcf = new Inputline(rfmix_opts.qvcf_fname, rfmix_opts.chromosome);
  vcf_skip_headers(qvcf);
  skip_to_chromosome(qvcf, rfmix_opts.chromosome);
  
  Inputline *rvcf = new Inputline(rfmix_opts.rvcf_fname, rfmix_opts.chromosome);
  vcf_skip_headers(rvcf);
  skip_to_chromosome(rvcf, rfmix_opts.chromosome);
  
  snp_t *snps = NULL;
  int n_snps = 0;
  char *q_chm, *r_chm;
  int q_pos, r_pos;
  
  pq = get_next_snp(qvcf, &q_chm, &q_pos);
  pr = get_next_snp(rvcf, &r_chm, &r_pos);

  double maf, miss;
  int mac;
  for(;;) {
    while(q_pos != -1 && strcmp(q_chm, rfmix_opts.chromosome) == 0 &&
	  q_pos < r_pos)
      pq = get_next_snp(qvcf, &q_chm, &q_pos);
    if (q_pos == -1 || strcmp(q_chm, rfmix_opts.chromosome) != 0) break;

    while(r_pos != -1 && strcmp(r_chm, rfmix_opts.chromosome) == 0 &&
	  r_pos < q_pos)
      pr = get_next_snp(rvcf, &r_chm, &r_pos);
    if (r_pos == -1 || strcmp(r_chm, rfmix_opts.chromosome) != 0) break;

    if (q_pos == r_pos) {
      if (q_pos < rfmix_opts.analyze_range[0] || q_pos > rfmix_opts.analyze_range[1]) {
	pq = get_next_snp(qvcf, &q_chm, &q_pos);      
	pr = get_next_snp(rvcf, &r_chm, &r_pos);
	continue;
      }
      
      /* Discard SNPs with too much missing data in either the query or the
	 reference files. If desired, insert minor allele frequency or minor
         allele count filters here */
      maf = vcf_snp_maf(&mac, &miss, pq);
      if (miss > rfmix_opts.maximum_missing_data_freq) {
	pq = get_next_snp(qvcf, &q_chm, &q_pos);      
	pr = get_next_snp(rvcf, &r_chm, &r_pos);
	continue;
      }

      maf = vcf_snp_maf(&mac, &miss, pr);
      if (miss > rfmix_opts.maximum_missing_data_freq) {
	pq = get_next_snp(qvcf, &q_chm, &q_pos);      
	pr = get_next_snp(rvcf, &r_chm, &r_pos);
	continue;
      }
      
      if (n_snps % SNP_ALLOC_STEP == 0)
	RA(snps, sizeof(snp_t)*(n_snps + SNP_ALLOC_STEP), snp_t);
      snps[n_snps].pos = q_pos;
      snps[n_snps].genetic_pos = input->genetic_map->translate_seqpos(q_pos);
      snps[n_snps].crf_index = -1;
      n_snps++;

      pq = get_next_snp(qvcf, &q_chm, &q_pos);      
      pr = get_next_snp(rvcf, &r_chm, &r_pos);
    }
  }

  input->snps = snps;
  input->n_snps = n_snps;

  delete qvcf;
  delete rvcf;
}

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

static void parse_alleles(input_t *input, Inputline *vcf, vcf_column_map_t *column_map,
			  int n_cols) {
  char *p, *q;
  char *chm;
  int pos;

  int snp_idx = 0;
  int n_unphased = 0;
  while(snp_idx < input->n_snps &&
	(p = get_next_snp(vcf, &chm, &pos)) != NULL &&
	strcmp(chm, rfmix_opts.chromosome) == 0) {
    if (input->snps[snp_idx].pos != pos) continue;

    int col_idx = VCF_LEAD_COLS;
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
	n_unphased++;
      }
      sample->haplotype[0][snp_idx] = get_allele(q[0]);
      sample->haplotype[1][snp_idx] = get_allele(q[2]);
      col_idx++;
    }
    snp_idx++;
  }

  if (n_unphased > 0) {
    fprintf(stderr,"\nWarning: %s - %d unphased genotypes treated as phased\n", vcf->fname, n_unphased);
  }
}

static void load_alleles(input_t *input) {
  Inputline *qvcf = new Inputline(rfmix_opts.qvcf_fname, rfmix_opts.chromosome);
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

  Inputline *rvcf = new Inputline(rfmix_opts.rvcf_fname, rfmix_opts.chromosome);
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

#define WINDOW_ALLOC_STEP 128
static void layout_random_forest(crf_window_t *crf, int n_crf, snp_t *snps, int n_snps,
				 double rf_spacing) {
  int minimum_snps = rfmix_opts.minimum_snps;
  int rf_start, rf_end;
  double s_gain, e_gain;

  /* spacing size less than 2.0 (limit could be larger) is intepreted as meaning a genetic window
     size in cM. Otherwise, it means a fixed number of SNPs. In either case, we widen the RF window
     by one SNP to the left or the right, which ever side would increase the genetic size of the 
     window the least. Thus, SNPs in the window are closest to the CRF window position taken as 
     the center point */
  if (rf_spacing <= 2.0) {  
    for(int w=0; w < n_crf; w++) {
      rf_start = crf[w].snp_idx;
      rf_end = crf[w].snp_idx;
      while((rf_start > 0 || rf_end < n_snps - 1) &&
	    (snps[rf_end].genetic_pos - snps[rf_start].genetic_pos < rf_spacing ||
	     rf_end - rf_start < rfmix_opts.minimum_snps)) {
	double last_window_size = snps[rf_end].genetic_pos - snps[rf_start].genetic_pos;

	s_gain = rf_start > 0 ? snps[rf_start].genetic_pos - snps[rf_start-1].genetic_pos : DBL_MAX;
	e_gain = rf_end < n_snps - 1 ? snps[rf_end+1].genetic_pos - snps[rf_end].genetic_pos : DBL_MAX;
	if (s_gain < e_gain) {
	  if (rf_spacing - last_window_size < s_gain && rf_end - rf_start >= minimum_snps) break; 
	  rf_start--;
	} else {
	  if (rf_spacing - last_window_size < e_gain && rf_end - rf_start >= minimum_snps) break;
	  rf_end++;
	}

	if (rf_end - rf_start >= 500) break;
      }

      if (w > 0 && rf_start > crf[w-1].rf_end_idx + 1) rf_start = crf[w-1].rf_end_idx + 1;
      crf[w].rf_start_idx = rf_start;
      crf[w].rf_end_idx = rf_end;
    }
  } else {
    for(int w=0; w < n_crf; w++) {
      rf_start = crf[w].snp_idx - rf_spacing / 2;
      if (rf_start < 0) rf_start = 0;
      rf_end = crf[w].snp_idx + rf_spacing / 2;
      if (rf_end >= n_snps) rf_end = n_snps - 1;

      if (w > 0 && rf_start > crf[w-1].rf_end_idx + 1) rf_start = crf[w-1].rf_end_idx + 1;
      /*      while((rf_start > 0 || rf_end < n_snps - 1) && rf_end - rf_start < rf_spacing) {

	s_gain = rf_start > 0 ? snps[rf_start].genetic_pos - snps[rf_start-1].genetic_pos : DBL_MAX;
	e_gain = rf_end < n_snps - 1 ? snps[rf_end+1].genetic_pos - snps[rf_end].genetic_pos : DBL_MAX;
	if (s_gain < e_gain) {
	  rf_start--;
	} else {
	  rf_end++;
	}
	}*/

      crf[w].rf_start_idx = rf_start;
      crf[w].rf_end_idx = rf_end;
    }
  }
  
}
static void set_crf_points(input_t *input) {

  fprintf(stderr,"\n   setting up CRF points and random forest windows... ");

  /* Local variable is needed for IDX(window,subpop) macro */
  int n_subpops = input->n_subpops;
  snp_t *snps = input->snps;
  int n_snps = input->n_snps;
  
  /* Determine the number of defined CRF points we have (CRF windows), the
     central SNP that defines each one, and the boundaries of the larger
     window used to source SNPs for the random forest classification */
  int w = 0;
  input->n_windows = 0;
  input->crf_windows = NULL;
  
  if (rfmix_opts.crf_spacing < 1.0) {
    
    MA(input->crf_windows, sizeof(crf_window_t)*(WINDOW_ALLOC_STEP), crf_window_t);
    input->crf_windows[0].genetic_pos = snps[0].genetic_pos;
    input->crf_windows[0].snp_idx = 0;

    w = 1;
    int s = 0;    
    while(s < n_snps - 1) {
      double target_pos = snps[s].genetic_pos + rfmix_opts.crf_spacing;
      if (snps[n_snps-1].genetic_pos - target_pos < rfmix_opts.crf_spacing / 2.)
	target_pos = snps[n_snps-1].genetic_pos;

      int t = s + 1;
      while(t < n_snps - 1 && snps[t].genetic_pos < target_pos) t++;
	    
      if (target_pos < snps[n_snps-1].genetic_pos && s - t > 1 && t < n_snps - 1 &&
	  snps[t].genetic_pos - target_pos < target_pos - snps[t-1].genetic_pos) t--;

      if (w % WINDOW_ALLOC_STEP == 0)
	RA(input->crf_windows, sizeof(crf_window_t)*(w + WINDOW_ALLOC_STEP), crf_window_t);

      input->crf_windows[w].snp_idx = t;
      input->crf_windows[w].genetic_pos = snps[t].genetic_pos;
      w++;
      s = t;
    }
    input->n_windows = w;
  } else {
    int s = 0;
    while(s < n_snps) {
      if (w % WINDOW_ALLOC_STEP == 0)
	RA(input->crf_windows, sizeof(crf_window_t)*(w + WINDOW_ALLOC_STEP), crf_window_t);

      input->crf_windows[w].snp_idx = s;
      input->crf_windows[w].genetic_pos = snps[s].genetic_pos;
      w++;
      
      s += rfmix_opts.crf_spacing;
    }
  }
  input->n_windows = w;

  fprintf(stderr,"\n   computing random forest window spacing overlay... ");
  layout_random_forest(input->crf_windows, input->n_windows, snps, n_snps, rfmix_opts.rf_window_size);
  /* Convert cM to M as we will always need in M in the CRF */
  for(w=0; w < input->n_windows; w++) input->crf_windows[w].genetic_pos /= 100.;
  
  /* Set up and initialize the current (starting) marginal probabilities for subpop
     assignment for each haplotype at each CRF window. These values start as 100%
     probability the haplotypes are from the apriori subpopulation for reference
     individuals, and just initialized to zero for all query individuals. These are
     calculated at each EM iteration by the Forward-Backward algorithm in the
     conditional random field code */
  fprintf(stderr,"\n   initializing apriori reference subpop across CRF... ");
  for(int k=0; k < input->n_samples; k++) {
    sample_t *sample = input->samples + k;
    
    for(int h=0; h < 2; h++) {
      MA(sample->current_p[h], sizeof(int16_t)*input->n_windows*n_subpops, int16_t);
      MA(sample->sis_p[h], sizeof(float)*input->n_windows*n_subpops, float);

      for(int i=0; i < input->n_windows; i++) {	
	for(int s=0; s < n_subpops; s++)
	  sample->current_p[h][ IDX(i,s) ] = ef16(0.0001/(n_subpops-1.));
	
	if (sample->apriori_subpop != -1)
	  sample->current_p[h][ IDX(i,sample->apriori_subpop) ] = ef16(0.9999);
      }
    }
  }

  fprintf(stderr,"\n   setting up random forest probability estimation arrays... ");
  for(int k=0; k < input->n_samples; k++) {
    sample_t *sample = input->samples + k;

    for(int h=0; h < 4; h++) {
      MA(sample->msp[h], sizeof(int8_t)*input->n_windows, int8_t);
      for(int i=0; i < input->n_windows; i++)
	sample->msp[h][i] = sample->apriori_subpop;

      MA(sample->est_p[h], sizeof(int16_t)*input->n_windows*n_subpops, int16_t);
    }
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

    /* Find and map out all the samples that we will be loading */
  fprintf(stderr,"Mapping samples ... ");
  load_samples(input);
  fprintf(stderr,"%d samples combined\n", input->n_samples);

  fprintf(stderr,"Scanning input VCFs for common SNPs on chromosome %s ...   ", rfmix_opts.chromosome);
  identify_common_snps(input);
  fprintf(stderr,"%d SNPs\n", input->n_snps);
  
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
  fprintf(stderr,"Defining and initializing conditional random field...   done\n");

  int64_t n_variant = 0;
  int64_t n_missing = 0;
  for(int i=0; i < input->n_samples; i++) {
    for(int h=0; h < 2; h++) {
      for(int s=0; s < input->n_snps; s++) {
	if (input->samples[i].haplotype[h][s] == 1) n_variant++;
	if (input->samples[i].haplotype[h][s] == 2) n_missing++;
      }
    }
  }

  fprintf(stderr,"%ld (%1.1f%%) variant alleles\t%ld (%1.1f%%) missing alleles\n",
	  n_variant, n_variant/((double) input->n_snps*2*input->n_samples)*100.,
	  n_missing, n_missing/((double) input->n_snps*2*input->n_samples)*100.);

  return input;
}


void free_input(input_t *input) {

  free(input->crf_windows);

  free(input->snps);
  input->n_snps = 0;
  
  for(int i=0; i < input->n_samples; i++) {
    sample_t *sample = input->samples + i;

    for(int h = 0; h < 2; h++) {
      free(sample->haplotype[h]);
      free(sample->current_p[h]);
      if (sample->ksp[h]) free(sample->ksp[h]);
      if (sample->sis_p[h]) free(sample->sis_p[h]);
    }

    for(int h=0; h < 4; h++) {
       free(sample->est_p[h]);
       free(sample->msp[h]);
    }   

    int *tmp = (int *) input->sample_hash->lookup(sample->sample_id);
    free(tmp);
    free(sample->sample_id);
  }


  delete input->sample_hash;
  free(input->samples);
  input->n_samples = 0;

  for(int i=0; i < input->n_subpops; i++)
    free(input->reference_subpops[i]);
  free(input->reference_subpops);

  delete input->genetic_map;
  free(input);

}
