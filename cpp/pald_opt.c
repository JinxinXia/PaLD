//
// Created by Yixin Zhang on 1/10/2021.
//

/*
params
D    in  distance matrix: D(x,y) is distance between x and y (symmetric,
         but assumed to be stored in full)
beta in  conflict focus parameter: z is in focus of (x,y) if
         min(d(z,x),d(z,y)) <= beta * d(x,y)
n    in  number of points
C    out cohesion matrix: C(x,z) is z's support for x
b    in  blocking parameter for cache efficiency
*/

#include <stdlib.h>
#include <string.h>

void pald_opt(double *D, double beta, int n, double *C, const int b) {
    // declare indices
    int x, y, z, i, j, xb, yb, ib;

    // pre-allocate conflict focus and distance cache blocks
    int *UXY = (int *) malloc(b * b * sizeof(int));
    double *DXY = (double *) malloc(b * b * sizeof(double));

    // initialize pointers for cache-block subcolumn vectors
    double *DXz, *DYz, *CXz, *CYz;

    // loop over blocks of points Y = (y,...,y+b-1)
    for (y = 0; y < n; y += b) {
        // define actual block size (for corner cases)
        yb = (b < n - y ? b : n - y);

        // loop over blocks of points X = (x,...,x+b-1)
        for (x = 0; x <= y; x += b) {
            // define actual block size (for corner cases)
            xb = (b < n - x ? b : n - x);

            // copy distances into cache block one column at a time
            for (j = 0; j < yb; j++) {
                // DXY(:,j) = D(x:x+xb,y+j) in off-diagonal case
                ib = (x == y ? j : xb); // handle diagonal blocks
                memcpy(DXY + j * xb, D + x + (y + j) * n, ib * sizeof(double));
            }

            // DEBUG: print out DXY cache block
            /*printf("x %d y %d xb %d yb %d\n",x,y,xb,yb);
            printf("\nDXY\n");
            for (int i = 0; i < xb; i++)
            {
                for (int j = 0; j < yb; j++)
                {
                    printf("%f ", DXY[i+j*xb]);
                }
                printf("\n");
            }
            printf("\n");*/

            // compute block's conflict focus sizes by looping over all points z
            memset(UXY, 0, b * b * sizeof(int)); // clear old values
            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D
            for (z = 0; z < n; z++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < yb; j++) {
                    ib = (x == y ? j : xb); // handle diagonal blocks
                    for (i = 0; i < ib; i++)
                        // DXY[i+j*xb] is distance between x+i and y+j
                        // DXz[i] is distance between x+i and z
                        // DYz[j] is distance between y+j and z

                        // determine if z is in conflict focus of x+i and y+j
                        if (DYz[j] <= beta * DXY[i + j * xb] || DXz[i] <= beta * DXY[i + j * xb])
                            UXY[i + j * xb]++;
                }

                // update pointers to subcolumns of D
                DXz += n;
                DYz += n;
            }

            // DEBUG: print out UXY cache block
            /*for (int i = 0; i < xb; i++)
            {
                for (int j = 0; j < yb; j++)
                {
                    printf("%d ", UXY[i+j*xb]);
                }
                printf("\n");
            }
            printf("\n");*/

            // update cohesion values according to conflicts between X and Y
            // by looping over all points z
            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D
            CXz = C + x;
            CYz = C + y; // init pointers to subcolumns of C
            for (z = 0; z < n; z++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < yb; j++) {
                    ib = (x == y ? j : xb); // handle diagonal blocks
                    for (i = 0; i < ib; i++) {
                        // DXY[i+j*xb] is distance between x+i and y+j
                        // DXz[i] is distance between x+i and z
                        // DYz[j] is distance between y+j and z

                        // check if z is in conflict of (x+i,y+j)
                        if (DYz[j] <= beta * DXY[i + j * xb] || DXz[i] <= beta * DXY[i + j * xb]) {
                            // z supports x+i
                            if (DXz[i] < DYz[j])
                                CXz[i] += 1.0 / UXY[i + j * xb];
                                // z supports y+j
                            else if (DYz[j] < DXz[i])
                                CYz[j] += 1.0 / UXY[i + j * xb];
                                // z splits its support
                            else {
                                CXz[i] += 0.5 / UXY[i + j * xb];
                                CYz[j] += 0.5 / UXY[i + j * xb];
                            }
                        }
                    }
                }

                // update pointers to subcolumns of D and C
                DXz += n;
                DYz += n;
                CXz += n;
                CYz += n;
            }
        }
    }

    // free up cache blocks
    free(DXY);
    free(UXY);
}


