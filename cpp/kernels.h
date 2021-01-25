#ifndef PALD_KERNELS_H
#define PALD_KERNELS_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// pald algorithms
void pald_orig(double *D, double beta, int n, double *C);
void pald_opt(double *D, double beta, int n, double *C, const int b);

#endif //PALD_KERNELS_H
