//
// Created by Yixin Zhang on 1/10/2021.
//

#ifndef PALD_UTILS_H
#define PALD_UTILS_H

#include <stdlib.h>
#include <stdio.h>

/*
This is a distance matrix generator

params
D           out empty pointer for the distance matrix:
                D(x,y) is distance between x and y (symmetric,
                but assumed to be stored in full)
edge_len    in  the edge length of a square matrix (total size edge_len^2)
upper_limit in  maximum distance to be generated
seed        in  random seed for matrix generation
*/

void dist_mat_gen(double *D, const int edge_len, const int upper_limit, const unsigned int seed) {

    srand(seed);
    int i, j;
    for (i = 0; i < edge_len - 1; i++) {
        for (j = i + 1; j < edge_len; j++) {
            int temp = rand() % upper_limit + 1;
            D[i * edge_len + j] = temp;
            D[j * edge_len + i] = temp;
        }
    }
    for (i = 0; i < edge_len; i++)
        D[i * edge_len + i] = 0;

}

#endif //PALD_UTILS_H
