#ifndef GENSAMPLES_H
#define GENSAMPLES_H

#include "rfmix.h"

void generate_simulated_samples(input_t *);
double score_msp(double ***r_m, double *r_ha, input_t *input);
void print_simulation_scoring_matrix(double **m, int n);
void free_simulation_scoring_matrix(double **m, int n);

#endif
