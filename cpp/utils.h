//
// Created by Yixin Zhang on 1/10/2021.
//

#ifndef PALD_UTILS_H
#define PALD_UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

void sym_mat_gen(double *D, const int edge_len, const int upper_limit, const unsigned int seed) {

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

//upper_limit is a place holder
void L1_dist_mat_gen2D(double *D, const int edge_len, int upper_limit, const unsigned int seed) {
    upper_limit = 35000;
    srand(seed);
    int i, j;
    int x[edge_len];
    int y[edge_len];

    for (i = 0; i < edge_len; i++)
        x[i] = rand() % upper_limit - 17500;

    for (i = 0; i < edge_len; i++)
        y[i] = rand() % upper_limit - 17500;

    for (i = 0; i < edge_len - 1; i++) {
        D[i * edge_len + i] = 0;
        for (j = i + 1; j < edge_len; j++) {
            register int x_dist = x[i] - x[j];
            register int y_dist = y[i] - y[j];
            register double temp = abs(x_dist) + abs((y_dist));
            //using L1 distance as example
            D[i * edge_len + j] = temp;
            D[j * edge_len + i] = temp;
        }
    }
    D[edge_len * edge_len - 1] = 0;


}

void L2_dist_mat_gen2D(double *D, const int edge_len, int upper_limit, const unsigned int seed) {
    upper_limit = 35000;
    srand(seed);
    int i, j;
    int x[edge_len];
    int y[edge_len];

    for (i = 0; i < edge_len; i++)
        x[i] = rand() % upper_limit - 17500;

    for (i = 0; i < edge_len; i++)
        y[i] = rand() % upper_limit - 17500;

    for (i = 0; i < edge_len - 1; i++) {
        D[i * edge_len + i] = 0;
        for (j = i + 1; j < edge_len; j++) {
            register int x_dist = x[i] - x[j];
            register int y_dist = y[i] - y[j];
            register double temp = sqrt(x_dist * x_dist + y_dist * y_dist);
            //using L1 distance as example
            D[i * edge_len + j] = temp;
            D[j * edge_len + i] = temp;
        }
    }
    D[edge_len * edge_len - 1] = 0;


}

#endif //PALD_UTILS_H
