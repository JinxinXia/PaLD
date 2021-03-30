#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "kernels.h"
#include "utils.h"

int main(int argc, char **argv) {

    //initializing testing environment spec
    int n, b, t, i;

    // initialize timers
    double start, /*time_seq,*/ time_par;
    
    if ((argc != 4) || !(n = atoi(argv[1])) || !(b = atoi(argv[2])) || !(t = atoi(argv[3])) ) {
        fprintf(stderr, "Usage: ./name mat_dim block_size num_threads\n");
        exit(-1);
    }
    unsigned int num_gen = n * n;

    float *C1 = calloc(num_gen, sizeof(float));
    float *D = malloc(sizeof(float) * num_gen);
    dist_mat_gen2D(D, n, 1, 10*n, 12345, '2');

    //computing C with parallel algorithm
    start =  omp_get_wtime();
    pald_opt_par(D, 1, n, C1, b, t);
    time_par =  omp_get_wtime() - start;

    //computing C with sequential alg 
    /* 
    start = omp_get_wtime();
    pald_opt(D, 1, n, C2, b);
    time_seq = omp_get_wtime() - start;
    */
    // compute max norm error between two cohesion matrices
    /*
    float d, maxdiff = 0.;
    for (i = 0; i < num_gen; i++) {
        d = fabs(C1[i]-C2[i]);
        maxdiff = d > maxdiff ? d : maxdiff;
    }
    printf("Maximum difference: %1.1e \n", maxdiff);
    */
    printf("%d Par time: %.3fs\n",n ,time_par);
    //printf("Parallel speedup is %f, efficiency is %2.1f\n",time_seq/time_par,time_seq/time_par/t*100);

    free(D);
    free(C1);
}
