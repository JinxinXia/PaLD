#include "opt_partial.h"

void pald_opt_new(double *D, double beta, int n, double *C, const int b) {
    // declare indices
    int x, y, z, i, j, xb, yb, ib;

    // pre-allocate conflict focus and distance cache blocks
    double *UXY = (double *) _mm_malloc(BLOCKSIZE * BLOCKSIZE * sizeof(double),64);
    double *DXY = (double *) _mm_malloc(BLOCKSIZE * BLOCKSIZE * sizeof(double),64);
    __assume_aligned(UXY,64);
    __assume_aligned(DXY,64);

    // initialize pointers for cache-block subcolumn vectors
    double *CXz, *CYz, *DXz, *DYz;
    __assume_aligned(DXz,64);
    __assume_aligned(DYz,64);
    __assume_aligned(CXz,64);
    __assume_aligned(CYz,64);

    //up_left main block
    for (y = 0; y < n; y += BLOCKSIZE) {

        // loop over blocks of points X = (x,...,x+b-1)
        for (x = 0; x <= y; x += BLOCKSIZE) {

            for (j = 0; j < BLOCKSIZE; j++) {
                // DXY(:,j) = D(x:x+xb,y+j) in off-diagonal case
                ib = (x == y ? j : BLOCKSIZE); // handle diagonal blocks
                memcpy(DXY + j * BLOCKSIZE, D + x + (y + j) * n, ib * sizeof(double));
            }

            /*
            if (x == y){
                for (j = 0; j < BLOCKSIZE; j++)
                    memcpy(DXY + j * BLOCKSIZE, D + x + (y + j) * n, j * sizeof(double));

            }
            else{
                for (j = 0; j < BLOCKSIZE; j++)
                    memcpy(DXY + j * BLOCKSIZE, D + x + (y + j) * n, BLOCKSIZE * sizeof(double));

            }*/

            // compute block's conflict focus sizes by looping over all points z
            memset(UXY, 0, BLOCKSIZE * BLOCKSIZE * sizeof(double)); // clear old values
            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D

            for (z = 0; z < n; z ++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < BLOCKSIZE; j++) {
                    ib = (x == y ? j : BLOCKSIZE); 
                    for (i = 0; i < ib; i++)
     
                        // determine if z is in conflict focus of x+i and y+j
                        if (DYz[j] <= beta * DXY[i + j * BLOCKSIZE] || DXz[i] <= beta * DXY[i + j * BLOCKSIZE])
                            UXY[i + j * BLOCKSIZE]++;
                }

                // update pointers to subcolumns of D
                DXz += n;
                DYz += n;
            }

            // update cohesion values according to conflicts between X and Y
            // by looping over all points z
            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D
            CXz = C + x;
            CYz = C + y; // init pointers to subcolumns of C

            for (z = 0; z < n; z++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < BLOCKSIZE; j++) {
                    ib = (x == y ? j : BLOCKSIZE); 
                    for (i = 0; i < ib; i++) {
                        
                        // check if z is in conflict of (x+i,y+j)
                        if (DYz[j] <= beta * DXY[i + j * BLOCKSIZE] || DXz[i] <= beta * DXY[i + j * BLOCKSIZE]) {
                            // z supports x+i
                            if (DXz[i] < DYz[j])
                                CXz[i] += 1.0 / UXY[i + j * BLOCKSIZE];
                                // z supports y+j
                            else if (DYz[j] < DXz[i])
                                CYz[j] += 1.0 / UXY[i + j * BLOCKSIZE];
                                // z splits its support
                            else {
                                CXz[i] += 0.5 / UXY[i + j * BLOCKSIZE];
                                CYz[j] += 0.5 / UXY[i + j * BLOCKSIZE];
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
    _mm_free(DXY);
    _mm_free(UXY);
}