#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

// linear indexing function assuming column major
int lin(int i, int j, int n) { return i + j * n; }

/*
params
D    in  distance matrix: D(x,y) is distance between x and y (symmetric)
beta in  conflict focus parameter: z is in focus of (x,y) if
         min(d(z,x),d(z,y)) <= beta * d(x,y)
n    in  number of points
C    out cohesion matrix: C(x,z) is z's support for x
*/
void pald_orig(double *D, double beta, int n, double *C) {
    // input checking
    if (beta < 0)
        fprintf(stderr, "beta must be positive\n");

    // loop over pairs of points x and y (only for x < y)
    for (int x = 0; x < n - 1; x++)
        for (int y = x + 1; y < n; y++) {
            int cfs = 0;                // conflict focus size of x,y
            double dxy = D[lin(x, y, n)]; // distance between x and y

            // loop over all points z to determine conflict focus size
            for (int z = 0; z < n; z++) {
                if (D[lin(z, x, n)] <= beta * dxy || D[lin(z, y, n)] <= beta * dxy)
                    cfs++;
            }

            // loop over all points z to determine contributions to x or y
            for (int z = 0; z < n; z++) {
                double dzx = D[lin(z, x, n)]; // dist between z and x
                double dzy = D[lin(z, y, n)]; // dist between z and y

                // z contributes to x or y only if in conflict focus
                if (dzx <= beta * dxy || dzy <= beta * dxy) {
                    if (dzx < dzy)
                        C[lin(x, z, n)] += 1.0 / cfs; // z closer to x than y
                    else if (dzy < dzx)
                        C[lin(y, z, n)] += 1.0 / cfs; // z closer to y than x

                    else {
                        // z equidistant to x and y
                        C[lin(x, z, n)] += 0.5 / cfs;
                        C[lin(y, z, n)] += 0.5 / cfs;
                    }
                }
            }
        }
}