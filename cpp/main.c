//
// Created by Yixin Zhang on 1/10/2021.
//
#include <stdio.h>
#include <stdlib.h>
#include "pald_opt.h"
#include "pald_orig.h"
#include "utils.h"

void print_out(int n, double *C) {
    printf("\n");
    int i, j;
    register int temp;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            temp = i * n + j;
            C[temp] /= (n - 1);
            printf("%.5f ", C[temp]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {

    //creating dist matrix
    int n = 10, cache_size = 2, i;
    unsigned int num_gen = n * n;
    double *C = calloc(num_gen, sizeof(double));

    if (n < 1)
        fprintf(stderr, "Matrix edge length must be positive");

    double *D = malloc(sizeof(double) * num_gen);
    dist_mat_gen(D, n, 30, 0);

    //print out dist matrix
    for (i = 0; i < num_gen; i++) {

        if (i % n == 0) {
            printf("\n");
        }
        printf("%2.0f ", D[i]);
    }

    //computing C with optimal block algorithm
    pald_opt(D, 1, n, C, cache_size);

    //print out block algorithm result
    print_out(n, C);

    free(C);

    //computing C with original algorithm
    C = calloc(num_gen, sizeof(double));
    pald_orig(D, 1, n, C);

    //print out result of original algorithm
    print_out(n, C);
    // print out for error checking

    free(D);
    free(C);
}