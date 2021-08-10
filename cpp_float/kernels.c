#include "kernels.h"

// linear indexing function assuming column major
inline int lin(int i, int j, int n) { return i + j * n; }

/*
params
D    in  distance matrix: D(x,y) is distance between x and y (symmetric)
beta in  conflict focus parameter: z is in focus of (x,y) if
         min(d(z,x),d(z,y)) <= beta * d(x,y)
n    in  number of points
C    out cohesion matrix: C(x,z) is z's support for x
*/
void pald_orig(float *D, float beta, int n, float *C) {
    // input checking
    if (beta < 0)
        fprintf(stderr, "beta must be positive\n");

    // loop over pairs of points x and y (only for x < y)
    for (int x = 0; x < n - 1; x++)
        for (int y = x + 1; y < n; y++) {
            int cfs = 0;                // conflict focus size of x,y
            float dxy = D[lin(x, y, n)]; // distance between x and y

            // loop over all points z to determine conflict focus size
            for (int z = 0; z < n; z++) {
                if (D[lin(z, x, n)] <= beta * dxy || D[lin(z, y, n)] <= beta * dxy)
                    cfs++;
            }

            // loop over all points z to determine contributions to x or y
            for (int z = 0; z < n; z++) {
                float dzx = D[lin(z, x, n)]; // dist between z and x
                float dzy = D[lin(z, y, n)]; // dist between z and y

                // z contributes to x or y only if in conflict focus
                if (dzx <= beta * dxy || dzy <= beta * dxy) {
                    if (dzx < dzy)
                        C[lin(x, z, n)] += 1.0f / cfs; // z closer to x than y
                    else if (dzy < dzx)
                        C[lin(y, z, n)] += 1.0f / cfs; // z closer to y than x

                    else {
                        // z equidistant to x and y
                        C[lin(x, z, n)] += 0.5f / cfs;
                        C[lin(y, z, n)] += 0.5f / cfs;
                    }
                }
            }
        }
}

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

void pald_opt_new(float* restrict D, float beta, int n, float* restrict C) {
    // declare indices
    int x, y, z, i, j, k, xb, yb, ib, divisible, remainder;

    remainder = n % BLOCKSIZE;
    divisible = n - remainder;

    float contains_one;
    // pre-allocate conflict focus and distance cache blocks
    float *UXY = (float *) _mm_malloc(BLOCKSIZE * BLOCKSIZE * sizeof(float),64);
    float *DXY = (float *) _mm_malloc(BLOCKSIZE * BLOCKSIZE * sizeof(float),64);
    float *in_range = (float *) _mm_malloc(BLOCKSIZE * sizeof(float),64);
    //handeling cases of unequal distance
    float *in_logic = (float *) _mm_malloc(BLOCKSIZE  * sizeof(float),64);
    float *in_logic_2 = (float *) _mm_malloc(BLOCKSIZE  * sizeof(float),64);
    float CYz_reduction = 0.f;

    __assume_aligned(C,64);
    __assume_aligned(D,64);
    __assume_aligned(UXY,64);
    __assume_aligned(DXY,64);
    __assume_aligned(in_range,64);
    __assume_aligned(in_logic,64);
    __assume_aligned(in_logic_2,64);


    // initialize pointers for cache-block subcolumn vectors
    float *CXz, *CYz, *DXz, *DYz;
    __assume_aligned(DXz,64);
    __assume_aligned(DYz,64);
    __assume_aligned(CXz,64);
    __assume_aligned(CYz,64);
    
    float dist_cutoff = 0., dist_cutoff_tmp = 0.;
    double temp =0., sum =0.;
    //up_left main block
    for (y = 0; y < divisible; y += BLOCKSIZE) {

        // loop over blocks of points X = (x,...,x+b-1)
        for (x = 0; x <= y; x += BLOCKSIZE) {

            for (j = 0; j < BLOCKSIZE; j++) {
                ib = (x == y ? j : BLOCKSIZE); // handle diagonal blocks
		        memcpy(DXY + j * BLOCKSIZE, D + x + (y + j) * n, ib * sizeof(float));
            }

            // compute block's conflict focus sizes by looping over all points z
            memset(UXY, 0, BLOCKSIZE * BLOCKSIZE * sizeof(float)); // clear old values
            // init pointers to subcolumns of D
            DXz = D + x;
            DYz = D + y; 

            //calculating conflict focus size
           
            for (z = 0; z < n; z ++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < BLOCKSIZE; j++) {
                    ib = (x == y ? j : BLOCKSIZE); 
                    for (i = 0; i < ib; i++) {
			            // determine if z is in conflict focus of x+i and y+j
                        if (DYz[j] <= beta * DXY[i + j * BLOCKSIZE] ||
                            DXz[i] <= beta * DXY[i + j * BLOCKSIZE]){
                            UXY[i + j * BLOCKSIZE]++;
			    
			            }
		            }
                }
                // update pointers to subcolumns of D
                DXz += n;
                DYz += n;
            }

            // pre-compute support from comflict focus size
            for (k=0;k<BLOCKSIZE*BLOCKSIZE;k++)
                UXY[k] = 1.0f/UXY[k];
            
            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D
            CXz = C + x;
            CYz = C + y; // init pointers to subcolumns of C

            // update cohesion values according to conflicts between X and Y
            
            for (z = 0; z < n; z++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < BLOCKSIZE; j++) {
                    ib = (x == y ? j : BLOCKSIZE); 
		            // masking z that supports y + j
                    temp = omp_get_wtime();
		            for (i = 0; i < ib; ++i) 
                        in_logic[i]= DXz[i] < DYz[j]? 1.0f:0.0f;
		            sum += omp_get_wtime()-temp;
                    // masking z that supports both
                    contains_one = 0.0f;
                    for (i = 0; i < ib; ++i) {
                        in_logic_2[i]= DXz[i] == DYz[j]? 1.0f:0.0f;
                        contains_one += in_logic_2[i];
                    }

                    // masking z in conflict focus
		            for (i = 0; i < ib; i++) {
                        if (DYz[j] <= beta * DXY[i + j * BLOCKSIZE] || 
                            DXz[i] <= beta * DXY[i + j * BLOCKSIZE]) 
                            in_range[i] = 1.0f;
                        else
                            in_range[i] = 0.0f;
		            }

                    // calculating support to x + j
		            for (i =0;i<ib ;++i)
                        CXz[i] +=  UXY[i + j * BLOCKSIZE] * in_range[i] * in_logic[i];
		            
		            // calculating support to x + i
		            CYz_reduction = 0;
		            for (i =0;i<ib;++i)
                        CYz_reduction +=  UXY[i + j * BLOCKSIZE] 
                        * in_range[i] * (1 - in_logic[i]);
		            
		            CYz[j] += CYz_reduction;
                    
		            // compensating float count in ties
                    if (contains_one > 0.5f){                      
                        CYz_reduction = CYz[j];
			            for (i = 0;i < ib; i++){
                            CXz[i] += 0.5f * UXY[i + j * BLOCKSIZE]*in_range[i]*in_logic_2[i];
                            CYz_reduction -= 0.5f * UXY[i + j * BLOCKSIZE]*in_range[i]*in_logic_2[i];
                        }
			            CYz[j] = CYz_reduction;
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
    //computing cases for non-even division of BLOCKSIZE
    printf("The time for parallel region is: %f\n", sum);
    if (remainder == 0){
        _mm_free(in_range);
        _mm_free(in_logic);
        _mm_free(in_logic_2);
        _mm_free(DXY);
        _mm_free(UXY);    
        return;
    }
    else{
        // blocks in upper-right rectangle, will not encounter diagonal
        for (x = 0; x < y; x += BLOCKSIZE) {
            for (j = 0; j < remainder; j++) {
		        memcpy(DXY + j * BLOCKSIZE, D + x + (y + j) * n, BLOCKSIZE * sizeof(float));
            }

            // compute block's conflict focus sizes by looping over all points z
            memset(UXY, 0, BLOCKSIZE * BLOCKSIZE * sizeof(float)); // clear old values
            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D
            for (z = 0; z < n; z ++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < remainder; j++) {
                    ib = (x == y ? j : BLOCKSIZE); 
                    for (i = 0; i < ib; i++) {
			// determine if z is in conflict focus of x+i and y+j
                        if (DYz[j] <= beta * DXY[i + j * BLOCKSIZE] || DXz[i] <= beta * DXY[i + j * BLOCKSIZE]){
                            UXY[i + j * BLOCKSIZE]++;
			    //in_range[i] = 1.0f;
			}
			
		   }
		    	
                }

                // update pointers to subcolumns of D
                DXz += n;
                DYz += n;
            }
            for (k=0;k<BLOCKSIZE*BLOCKSIZE;k++)
                UXY[k] = 1.0f/UXY[k];


            // update cohesion values according to conflicts between X and Y
            // by looping over all points z
            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D
            CXz = C + x;
            CYz = C + y; // init pointers to subcolumns of C

            for (z = 0; z < n; z++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < remainder; j++) {
                    ib = (x == y ? j : BLOCKSIZE); 
		    
		   // z supports y+j
		    for (i = 0; i < ib; ++i) {
                        in_logic[i]= DXz[i] < DYz[j]? 1.0f:0.0f;
		    }

                    contains_one = 0.0f;
		    for (i = 0; i < ib; ++i) {
                        in_logic_2[i]= DXz[i] == DYz[j]? 1.0f:0.0f;
			contains_one += in_logic_2[i];
		    }

		    for (i = 0; i < ib; i++) {
			//dist_cutoff = beta * DXY[i + j * BLOCKSIZE]; 
                        //if (DYz[j] <= dist_cutoff || DXz[i] <= dist_cutoff)
                        if (DYz[j] <= beta * DXY[i + j * BLOCKSIZE] || DXz[i] <= beta * DXY[i + j * BLOCKSIZE]) 
		      	    in_range[i] = 1.0f;
                        else
                            in_range[i] = 0.0f;
		    }

		    for (i =0;i<ib ;++i){
                        CXz[i] +=  UXY[i + j * BLOCKSIZE]*in_range[i]*(in_logic[i]);
                        //CYz[j] +=  UXY[i + j * BLOCKSIZE]*in_range[i]*(1 - in_logic[i]);
		    }

		    // z supports x+i
		    CYz_reduction = 0;
		    for (i =0;i<ib;++i){
                        CYz_reduction +=  UXY[i + j * BLOCKSIZE]*in_range[i]*(1 - in_logic[i]);
		    }
		    CYz[j] += CYz_reduction;
                    
		    // z is evenly divided
                    if (contains_one > 0.5f){                      
                        CYz_reduction = CYz[j];
			for (i = 0;i < ib; i++){
                            CXz[i] += 0.5f * UXY[i + j * BLOCKSIZE]*in_range[i]*in_logic_2[i];
                            CYz_reduction -= 0.5f * UXY[i + j * BLOCKSIZE]*in_range[i]*in_logic_2[i];
                        }
			CYz[j] = CYz_reduction;
                    }
                    
                }

                // update pointers to subcolumns of D and C
                DXz += n;
                DYz += n;
                CXz += n;
                CYz += n;
            }
        }


        for (j = 0; j < remainder; j++) {
            // DXY(:,j) = D(x:x+xb,y+j) in off-diagonal case
            ib = (x == y ? j : remainder); // handle diagonal blocks
	        memcpy(DXY + j * remainder, D + x + (y + j) * n, ib * sizeof(float));
        }

            // compute block's conflict focus sizes by looping over all points z
        memset(UXY, 0, BLOCKSIZE * BLOCKSIZE * sizeof(float)); // clear old values
        DXz = D + x;
        DYz = D + y; // init pointers to subcolumns of D

        for (z = 0; z < n; z ++) {
            // loop over all (i,j) pairs in block
            for (j = 0; j < remainder; j++) {
                ib = (x == y ? j : remainder); 
                for (i = 0; i < ib; i++) {
		            // determine if z is in conflict focus of x+i and y+j
                    if (DYz[j] <= beta * DXY[i + j * remainder] || 
                        DXz[i] <= beta * DXY[i + j * remainder])
                        UXY[i + j * remainder]++;    
		        }
		    	
            }

                // update pointers to subcolumns of D
                DXz += n;
                DYz += n;
        }
            for (k=0;k<BLOCKSIZE*BLOCKSIZE;k++)
                UXY[k] = 1.0f/UXY[k];


            // update cohesion values according to conflicts between X and Y
            // by looping over all points z
            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D
            CXz = C + x;
            CYz = C + y; // init pointers to subcolumns of C

            for (z = 0; z < n; z++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < remainder; j++) {
                    ib = (x == y ? j : remainder); 
		    
		   // z supports y+j
		    for (i = 0; i < ib; ++i) {
                        in_logic[i]= DXz[i] < DYz[j]? 1.0f:0.0f;
		    }

                    contains_one = 0.0f;
		    for (i = 0; i < ib; ++i) {
                        in_logic_2[i]= DXz[i] == DYz[j]? 1.0f:0.0f;
			contains_one += in_logic_2[i];
		    }

		    for (i = 0; i < ib; i++) {
			//dist_cutoff = beta * DXY[i + j * BLOCKSIZE]; 
                        //if (DYz[j] <= dist_cutoff || DXz[i] <= dist_cutoff)
                        if (DYz[j] <= beta * DXY[i + j * remainder] || DXz[i] <= beta * DXY[i + j * remainder]) 
		      	    in_range[i] = 1.0f;
                        else
                            in_range[i] = 0.0f;
		    }

		    for (i =0;i<ib ;++i){
                        CXz[i] +=  UXY[i + j * remainder]*in_range[i]*(in_logic[i]);
                        //CYz[j] +=  UXY[i + j * BLOCKSIZE]*in_range[i]*(1 - in_logic[i]);
		    }

		    // z supports x+i
		    CYz_reduction = 0;
		    for (i =0;i<ib;++i){
                        CYz_reduction +=  UXY[i + j * remainder]*in_range[i]*(1 - in_logic[i]);
		    }
		    CYz[j] += CYz_reduction;
                    
		    // z is evenly divided
                    if (contains_one > 0.5f){                      
                        CYz_reduction = CYz[j];
			for (i = 0;i < ib; i++){
                            CXz[i] += 0.5f * UXY[i + j * remainder]*in_range[i]*in_logic_2[i];
                            CYz_reduction -= 0.5f * UXY[i + j * remainder]*in_range[i]*in_logic_2[i];
                        }
			CYz[j] = CYz_reduction;
                    }
                    
                }
                // update pointers to subcolumns of D and C
                DXz += n;
                DYz += n;
                CXz += n;
                CYz += n;
            }
    }
    //printf("x = %d, y = %d \n", x,  y);
    
    // free up cache blocks
    _mm_free(in_range);
    _mm_free(in_logic);
    _mm_free(in_logic_2);
    _mm_free(DXY);
    _mm_free(UXY);
}

void pald_opt_new_par(float *D, float beta, int n, float *C, int t){
    // declare indices
    
       // declare indices
    int x, y, z, j, k, xb, yb, ib, divisible, remainder;

    remainder = n % BLOCKSIZE;
    divisible = n - remainder;

    float contains_one;
    // pre-allocate conflict focus and distance cache blocks
    float *UXY = (float *) _mm_malloc(BLOCKSIZE * BLOCKSIZE * sizeof(float),64);
    float *DXY = (float *) _mm_malloc(BLOCKSIZE * BLOCKSIZE * sizeof(float),64);
    float *in_range = (float *) _mm_malloc(BLOCKSIZE * sizeof(float),64);
    //handeling cases of unequal distance
    float *in_logic = (float *) _mm_malloc(BLOCKSIZE  * sizeof(float),64);
    float *in_logic_2 = (float *) _mm_malloc(BLOCKSIZE  * sizeof(float),64);
    float CYz_reduction = 0.f;

    __assume_aligned(C,64);
    __assume_aligned(D,64);
    __assume_aligned(UXY,64);
    __assume_aligned(DXY,64);
    __assume_aligned(in_range,64);
    __assume_aligned(in_logic,64);
    __assume_aligned(in_logic_2,64);


    // initialize pointers for cache-block subcolumn vectors
    float *CXz, *CYz, *DXz, *DYz;
    __assume_aligned(DXz,64);
    __assume_aligned(DYz,64);
    __assume_aligned(CXz,64);
    __assume_aligned(CYz,64);
    
    float dist_cutoff = 0., dist_cutoff_tmp = 0.;
    double temp =0., sum =0.;
    //up_left main block
        for (y = 0; y < divisible; y += PBLOCKSIZE) {

        // loop over blocks of points X = (x,...,x+b-1)
        
        for (x = 0; x <= y; x += PBLOCKSIZE) {

            for (int j = 0; j < PBLOCKSIZE; j++) {
                int ib = (x == y ? j : PBLOCKSIZE); // handle diagonal blocks
		        memcpy(DXY + j * PBLOCKSIZE, D + x + (y + j) * n, ib * sizeof(float));
            }

            // compute block's conflict focus sizes by looping over all points z
            memset(UXY, 0, PBLOCKSIZE * PBLOCKSIZE * sizeof(float)); // clear old values
            // init pointers to subcolumns of D
         
            //calculating conflict focus size
            

            #pragma omp parallel for num_threads(t) reduction(+:UXY[:PBLOCKSIZE*PBLOCKSIZE]) schedule(static,n/t)
            for (int z = 0; z < n; z ++) {
                float* DXz1 = D + x + z*n;
                float* DYz1 = D + y + z*n;
                __assume_aligned(DXz1,64);
                __assume_aligned(DYz1,64);

                // loop over all (i,j) pairs in block
                for (int j = 0; j < PBLOCKSIZE; j++) {
                    int ib = (x == y ? j : PBLOCKSIZE); 
                    for (int i = 0; i < ib; i++) {
			            // determine if z is in conflict focus of x+i and y+j
                        if (DYz1[j] <= beta * DXY[i + j * PBLOCKSIZE] ||
                            DXz1[i] <= beta * DXY[i + j * PBLOCKSIZE]){
                            UXY[i + j * PBLOCKSIZE]++;
			    
			            }
		            }
                }
            }

            // pre-compute support from comflict focus size
            #pragma omp parallel for num_threads(t) schedule(static,n/t)
            for (int k = 0;k<PBLOCKSIZE*PBLOCKSIZE;k++)
                UXY[k] = 1.0f/UXY[k];

            DXz = D + x;
            DYz = D + y; // init pointers to subcolumns of D
            CXz = C + x;
            CYz = C + y; // init pointers to subcolumns of C

            // update cohesion values according to conflicts between X and Y
            for (z = 0; z < n; z++) {
                // loop over all (i,j) pairs in block
                for (j = 0; j < BLOCKSIZE; j++) {
                    ib = (x == y ? j : BLOCKSIZE); 
		            // masking z that supports y + j
                    temp = omp_get_wtime();
                    //#pragma omp parallel for num_threads(t) schedule(static,n/t)
		            for (int i = 0; i < ib; i++) 
                        in_logic[i]= DXz[i] < DYz[j]? 1.0f:0.0f;
		            
                    sum += omp_get_wtime() - temp;
                    // masking z that supports both
                    contains_one = 0.0f;
                    for (int i = 0; i < ib; i++) {
                        in_logic_2[i]= DXz[i] == DYz[j]? 1.0f:0.0f;
                        contains_one += in_logic_2[i];
                    }

                    // masking z in conflict focus
		            for (int i = 0; i < ib; i++) {
                        if (DYz[j] <= beta * DXY[i + j * BLOCKSIZE] || 
                            DXz[i] <= beta * DXY[i + j * BLOCKSIZE]) 
                            in_range[i] = 1.0f;
                        else
                            in_range[i] = 0.0f;
		            }

                    // calculating support to x + j
		            for (int i = 0;i < ib; i++)
                        CXz[i] +=  UXY[i + j * BLOCKSIZE] * in_range[i] * in_logic[i];
		            
		            // calculating support to x + i
		            CYz_reduction = 0;
		            for (int i =0; i < ib; i++)
                        CYz_reduction +=  UXY[i + j * BLOCKSIZE] 
                        * in_range[i] * (1 - in_logic[i]);
		            
		            CYz[j] += CYz_reduction;
                    
		            // compensating float count in ties
                    if (contains_one > 0.5f){                      
                        CYz_reduction = CYz[j];
			            for (int i = 0;i < ib; i++){
                            CXz[i] += 0.5f * UXY[i + j * BLOCKSIZE]*in_range[i]*in_logic_2[i];
                            CYz_reduction -= 0.5f * UXY[i + j * BLOCKSIZE]*in_range[i]*in_logic_2[i];
                        }
			            CYz[j] = CYz_reduction;
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
    //computing cases for non-even division of BLOCKSIZE
    // free up cache blocks
    printf("The time for parallel region (Par) is: %f\n", sum);
    _mm_free(in_range);
    _mm_free(in_logic);
    _mm_free(in_logic_2);
    _mm_free(DXY);
    _mm_free(UXY);
}

void pald_opt(float *D, float beta, int n, float *C, const int b) {
    // declare indices
    int x, y, z, i, j, xb, yb, ib;

    // pre-allocate conflict focus and distance cache blocks
    int *UXY = (int *) malloc(b * b * sizeof(int));
    float *DXY = (float *) malloc(b * b * sizeof(float));

    // initialize pointers for cache-block subcolumn vectors
    float *DXz, *DYz, *CXz, *CYz;

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
                memcpy(DXY + j * xb, D + x + (y + j) * n, ib * sizeof(float));
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
                                CXz[i] += 1.0f / UXY[i + j * xb];
                                // z supports y+j
                            else if (DYz[j] < DXz[i])
                                CYz[j] += 1.0f / UXY[i + j * xb];
                                // z splits its support
                            else {
                                CXz[i] += 0.5f / UXY[i + j * xb];
                                CYz[j] += 0.5f / UXY[i + j * xb];
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

/*
params
D    in  distance matrix: D(x,y) is distance between x and y (symmetric,
         but assumed to be stored in full)
beta in  conflict focus parameter: z is in focus of (x,y) if
         min(d(z,x),d(z,y)) <= beta * d(x,y)
n    in  number of points
C    out cohesion matrix: C(x,z) is z's support for x
b    in  blocking parameter for cache efficiency
t    in  number of OMP threads to use
*/
void pald_opt_par(float *D, float beta, int n, float *C, const int b, int t) {

    // pre-allocate conflict focus and distance cache blocks
    int *UXY = (int *) malloc(b * b * sizeof(int));
    float *DXY = (float *) malloc(b * b * sizeof(float));

    // loop over blocks of points Y = (y,...,y+b-1)
    for (int y = 0; y < n; y += b) {
        // define actual block size (for corner cases)
        int yb = (b < n - y ? b : n - y);

        // loop over blocks of points X = (x,...,x+b-1)
        for (int x = 0; x <= y; x += b) {
            // define actual block size (for corner cases)
            int xb = (b < n - x ? b : n - x);

            // copy distances into cache block one column at a time
            for (int j = 0; j < yb; j++) {
                // DXY(:,j) = D(x:x+xb,y+j) in off-diagonal case
                int ib = (x == y ? j : xb); // handle diagonal blocks
                memcpy(DXY + j * xb, D + x + (y + j) * n, ib * sizeof(float));
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
            #pragma omp parallel for num_threads(t) reduction(+:UXY[:b*b])
            for (int z = 0; z < n; z++) {
                // set pointers to subcolumns of D
                float* DXz = D + x + z*n;
                float* DYz = D + y + z*n; 
                // loop over all (i,j) pairs in block
                for (int j = 0; j < yb; j++) {
                    int ib = (x == y ? j : xb); // handle diagonal blocks
                    for (int i = 0; i < ib; i++)
                        // DXY[i+j*xb] is distance between x+i and y+j
                        // DXz[i] is distance between x+i and z
                        // DYz[j] is distance between y+j and z

                        // determine if z is in conflict focus of x+i and y+j
                        if (DYz[j] <= beta * DXY[i + j * xb] || DXz[i] <= beta * DXY[i + j * xb])
                            UXY[i + j * xb]++;
                }
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
            #pragma omp parallel for num_threads(t)
            for (int z = 0; z < n; z++) {
                // set pointers to subcolumns of D
                float* DXz = D + x + z*n;
                float* DYz = D + y + z*n; 
                // set pointers to subcolumns of C
                float* CXz = C + x + z*n;
                float* CYz = C + y + z*n; 
                // loop over all (i,j) pairs in block
                for (int j = 0; j < yb; j++) {
                    int ib = (x == y ? j : xb); // handle diagonal blocks
                    for (int i = 0; i < ib; i++) {
                        // DXY[i+j*xb] is distance between x+i and y+j
                        // DXz[i] is distance between x+i and z
                        // DYz[j] is distance between y+j and z

                        // check if z is in conflict of (x+i,y+j)
                        if (DYz[j] <= beta * DXY[i + j * xb] || DXz[i] <= beta * DXY[i + j * xb]) {
                            // z supports x+i
                            if (DXz[i] < DYz[j])
                                CXz[i] += 1.0f / UXY[i + j * xb];
                                // z supports y+j
                            else if (DYz[j] < DXz[i])
                                CYz[j] += 1.0f / UXY[i + j * xb];
                                // z splits its support
                            else {
                                CXz[i] += 0.5f / UXY[i + j * xb];
                                CYz[j] += 0.5f / UXY[i + j * xb];
                            }
                        }
                    }
                }
            }
        }
    }

    // free up cache blocks
    free(DXY);
    free(UXY);

    // print out timing results before returning
}
