#include "opt_partial.h"

void pald_opt_new(double* restrict D, double beta, int n, double* restrict C) {
    // declare indices
    int x, y, z, i, j, k, xb, yb, ib, divisible, remainder;

    remainder = n % BLOCKSIZE;
    divisible = n - remainder;

    double contains_one;
    // pre-allocate conflict focus and distance cache blocks
    double *UXY = (double *) _mm_malloc(BLOCKSIZE * BLOCKSIZE * sizeof(double),64);
    double *DXY = (double *) _mm_malloc(BLOCKSIZE * BLOCKSIZE * sizeof(double),64);
    double *in_range = (double *) _mm_malloc(BLOCKSIZE * sizeof(double),64);
    //handeling cases of unequal distance
    double *in_logic = (double *) _mm_malloc(BLOCKSIZE  * sizeof(double),64);
    double *in_logic_2 = (double *) _mm_malloc(BLOCKSIZE  * sizeof(double),64);
    double CYz_reduction = 0.f;

    __assume_aligned(C,64);
    __assume_aligned(D,64);
    __assume_aligned(UXY,64);
    __assume_aligned(DXY,64);
    __assume_aligned(in_range,64);
    __assume_aligned(in_logic,64);
    __assume_aligned(in_logic_2,64);


    // initialize pointers for cache-block subcolumn vectors
    double *CXz, *CYz, *DXz, *DYz;
    __assume_aligned(DXz,64);
    __assume_aligned(DYz,64);
    __assume_aligned(CXz,64);
    __assume_aligned(CYz,64);
    
    double dist_cutoff = 0., dist_cutoff_tmp = 0.;
    
    //up_left main block
    for (y = 0; y < divisible; y += BLOCKSIZE) {

        // loop over blocks of points X = (x,...,x+b-1)
        for (x = 0; x <= y; x += BLOCKSIZE) {

            for (j = 0; j < BLOCKSIZE; j++) {
                ib = (x == y ? j : BLOCKSIZE); // handle diagonal blocks
		        memcpy(DXY + j * BLOCKSIZE, D + x + (y + j) * n, ib * sizeof(double));
            }

            // compute block's conflict focus sizes by looping over all points z
            memset(UXY, 0, BLOCKSIZE * BLOCKSIZE * sizeof(double)); // clear old values
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
		            for (i = 0; i < ib; ++i) 
                        in_logic[i]= DXz[i] < DYz[j]? 1.0f:0.0f;
		            
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
                    
		            // compensating double count in ties
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
		        memcpy(DXY + j * BLOCKSIZE, D + x + (y + j) * n, BLOCKSIZE * sizeof(double));
            }

            // compute block's conflict focus sizes by looping over all points z
            memset(UXY, 0, BLOCKSIZE * BLOCKSIZE * sizeof(double)); // clear old values
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
	        memcpy(DXY + j * remainder, D + x + (y + j) * n, ib * sizeof(double));
        }

            // compute block's conflict focus sizes by looping over all points z
        memset(UXY, 0, BLOCKSIZE * BLOCKSIZE * sizeof(double)); // clear old values
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
