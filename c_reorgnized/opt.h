#ifndef PALD_OPT_H
#define PALD_OPT_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

void pald_opt(double *D, double beta, int n, double *C, const int b);
void pald_opt_par(double *D, double beta, int n, double *C, const int b, int num_threads);
#endif
