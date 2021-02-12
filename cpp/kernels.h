#ifndef PALD_KERNELS_H
#define PALD_KERNELS_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

// pald algorithms
void pald_orig(double *D, double beta, int n, double *C);
void pald_opt(double *D, double beta, int n, double *C, const int b);
void pald_opt_par(double *D, double beta, int n, double *C, const int b, int num_threads);

#endif //PALD_KERNELS_H
