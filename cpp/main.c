//
// Created by Yixin Zhang on 1/10/2021.
//
#include <stdio.h>
#include <stdlib.h>
#include "pald_opt.c"
#include "pald_orig.c"
#include "utils.c"

void print_out(int n, double *C) {
    printf("\n");
    int i, j;
    register int temp;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            temp = j * n + i;
            C[temp] /= (n - 1);
            printf("%.7f ", C[temp]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {

    //initializing testing environment spec
    int n, cache_size, i;
    if ((argc != 2 && argc != 3) || !(n = atoi(argv[1]))) {
        fprintf(stderr, "Usage: ./name distance_mat_size block_size\n");
        exit(-1);
    }

    cache_size = argc == 2 ? 2 : atoi(argv[2]);

    unsigned int num_gen = n * n;

    double *C1 = calloc(num_gen, sizeof(double));
    double *C2 = calloc(num_gen, sizeof(double));
    double *D = malloc(sizeof(double) * num_gen);
    dist_mat_gen2D(D, n, 1, 5, 12345, '1');

    //print out dist matrix
    /*for (i = 0; i < num_gen; i++) {

        if (i % n == 0) {
            printf("\n");
        }
        printf("%.2f ", D[i]);
    }*/

    //computing C with optimal block algorithm
    pald_opt(D, 1, n, C1, cache_size);

    //print out block algorithm result
    //print_out(n, C1);


    //computing C with original algorithm
    pald_orig(D, 1, n, C2);

    //print out result of original algorithm
    //print_out(n, C2);
    // print out for error checking

    // compute max norm error between two cohesion matrices
    double diff, maxdiff = 0.;
    for (i = 0; i < num_gen; i++) {
        diff = fabs(C1[i]-C2[i]);
        maxdiff = diff > maxdiff ? diff : maxdiff;
    }
    printf("Maximum difference: %1.1e \n", maxdiff);

    free(D);
    free(C2);
    free(C1);
}
