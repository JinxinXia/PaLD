#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

void optimal_matmul_block(double *C, double *D, double beta, int edge_len, const int b)
{
    int i, j;
    for (i = 0; i < edge_len; i += b)
    {
        int Ib = min(i + b - 1, edge_len);
        for (j = i; j < edge_len; j += b)
        {
            int Jb = min(j + b - 1, edge_len);

            int dim_row = min(b, edge_len - i + 1);
            int dim_col = min(b, edge_len - j + 1);
            double *Uij = calloc(dim_col * dim_row, sizeof(double));
            double Dij = D[Ib * edge_len + Jb];
            int k;

            for (k = 0; k < edge_len; k++)
            {
                if (fmin(D[Ib * edge_len + k], D[k * edge_len + Jb]) < Dij)
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
                /*for (k = 0; i < edge_len; k++)
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
}

int main(int argc, char **argv)
{

    int edge_len = 4;
    double *C = calloc(edge_len * edge_len, sizeof(double));

    double D[] = {0, 1, 2, 3, 1, 0, 4, 5, 2, 4, 0, 6, 3, 5, 6, 0};
    orig_contribute(C, D, 1, edge_len, 3);

    free(C);
}
