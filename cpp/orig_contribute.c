#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

void orig_contribute(double *C, double *D, int *F, double beta, int edge_len)
{
    int x, y;
    if (beta < 0)
        fprintf(stderr, "beta must be positive\n");

    for (x = 0; x < edge_len - 1; x++)
    {
        for (y = x + 1; y < edge_len; y++)
        {
            register int num_involve = 0;
            register int other = 0;
            register int px, py;
            register int x_offset = x * edge_len;
            register int y_offset = y * edge_len;

            for (other = 0; other < edge_len; other++)
            {
                px = x_offset + other;
                py = y_offset + other;

                if (D[px] <= beta * D[x_offset + y] || D[py] <= beta * D[x_offset + y])
                    num_involve++;
                //calculate local depth
            }
            F[x_offset + y] = num_involve;

            for (other = 0; other < edge_len; other++)
            {
                px = x_offset + other;
                py = y_offset + other;

                if (D[px] < D[py])
                    C[px] += 1.0 / num_involve;
                else if (D[px] > D[py])
                    C[py] += 1.0 / num_involve;
                else
                {
                    C[px] += 0.5 / num_involve;
                    C[py] += 0.5 / num_involve;
                }
            }
        }
    }
    printf("\n");
    int i, j;
    register int temp;
    for (i = 0; i < edge_len; i++)
    {
        for (j = 0; j < edge_len; j++)
        {
            temp = i * edge_len + j;
            C[temp] /= (edge_len - 1);
            C[temp] *= 3;
            printf("%.5f ", C[temp]);
        }
        printf("\n");
    }

    for (i = 0; i < edge_len; i++)
    {
        for (j = 0; j < edge_len; j++)
        {
            temp = i * edge_len + j;

            printf("%d ", F[temp]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv)
{

    int edge_len = 4;
    double *C = calloc(edge_len * edge_len, sizeof(double));
    int *F = calloc(edge_len * edge_len, sizeof(int));

    double D[] = {0, 1, 2, 3, 1, 0, 4, 5, 2, 4, 0, 6, 3, 5, 6, 0};
    orig_contribute(C, D, F, 1, edge_len);

    free(F);
    free(C);
}
/*int main(int argc, char **argv)
{
    FILE *fp;
    float beta;
    unsigned int edge_length;
    if (argc != 4 || !(fp = fopen(argv[1], 'r')) || !(beta = atof(argv[2])) || !(edge_length = atoi(argv[3])))
    {
        fprintf(stderr, "Usage: ./orig_contribute Distance_matrix.txt beta edge_length_of_distance_mat\n");
        if (beta < 0)
            fprintf(stderr, "beta must be positive\n");
        exit(-1);
    }
    int mat_size = edge_length * edge_length;
    float *D = malloc(sizeof(float) * mat_size);

    read_D(fp, D);
    if (-1 * is_symmetric(D, edge_length))
    {
        fprintf(stderr, 'The distance matrix D must be symmetric!\n');
        exit(-1);
    }
    float *C = malloc(sizeof(float) * mat_size);
    float *F = malloc(sizeof(float) * mat_size);

    return 0;
}

int is_symmetric(double *D, int edge_len)
{
    int i, j;
    for (i = 0; i < edge_len - 1; i++)
        for (j = i + 1; j < edge_len; j++)
            if (D[i * edge_len + j] != D[j * edge_len + i])
                return -1;

    return 1;
}
void read_D(FILE *fp, double *D)
{
    register unsigned char temp;
    int j = 0;
    int i = 0;

    while (temp = fgetc(fp))
    {
        unsigned char *temp_char[128];

        if (isdigit(temp) || temp == '.')
        {
            temp_char[i] = temp;
            i++;
        }

        else
        {
            D[j] = atof(temp_char);
            i = 0;
            j++;
        }
    }
}*/