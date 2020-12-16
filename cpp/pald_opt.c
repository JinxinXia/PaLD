#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/*
params
D    in  distance matrix: D(x,y) is distance between x and y (symmetric)
beta in  conflict focus parameter: z is in focus of (x,y) if 
         min(d(z,x),d(z,y)) <= beta * d(x,y)
n    in  number of points
C    out cohesion matrix: C(x,z) is z's support for x
b    in  blocking parameter for cache efficiency
*/
void pald_opt(double *D, double beta, int n, double* C, const int b)
{
    // pre-allocate conflict focus block
    double *Uij = calloc(b * b, sizeof(double));

    // loop over blocks of points I = (i,...,i+b-1)
    for (int i = 0; i < n; i += b)
    {
        int Ib = min(i + b - 1, n);
        // loop over blocks of points J = (j,...,j+b-1)
        for (j = i; j < n; j += b)
        {
            int Jb = min(j + b - 1, n);

            int dim_row = min(b, n - i + 1);
            int dim_col = min(b, n - j + 1);
            double Dij = D[Ib * n + Jb];
            int k;

            for (k = 0; k < n; k++)
            {
                if (fmin(D[Ib * n + k], D[k * n + Jb]) < Dij)
                {
                    int others;
                    for (others = 0; others < dim_row * dim_col; others++)
                        Uij[others]++;
                }
            }
            if (i == j)
            {
                int others;
                for (others = 0; others < dim_row; others++)
                {
                    Uij[others * dim_row + others] = DBL_MAX;
                }
                int k;

                //FIXME
                /*for (k = 0; i < n; k++)
                {
                    C[Ib * dim_row + k] +=
                }*/
            }
            else
            {
                /* code */
            }
        }
    }

    // free up cache block
    free(Uij);
}

int main(int argc, char **argv)
{

    int n = 4;
    double *C = calloc(n * n, sizeof(double));

    double D[] = {0, 1, 2, 3, 1, 0, 4, 5, 2, 4, 0, 6, 3, 5, 6, 0};
    pald_opt(D, 1, n, C, 3);

    free(C);
}
