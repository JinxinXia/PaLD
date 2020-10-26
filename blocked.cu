#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#define TILE_WIDTH 32


__global__ void matmul(float* A, float* B, float* C, int width){
	//using shared memory method
	__shared__ float submatrix_A[TILE_WIDTH][TILE_WIDTH];
	__shared__ float submatrix_B[TILE_WIDTH][TILE_WIDTH];
	
	int bx = blockIdx.x; int by = blockIdx.y;
	int tx = threadIdx.x; int ty = threadIdx.y;
	int row = bx*TILE_WIDTH + tx;
	int col = by*TILE_WIDTH + ty;
	float temp=0;
	
	
	for(int phase=0; phase < ceil(1.0*width/TILE_WIDTH); phase++){
		
		//Prepare shared memory at each phase
		if ( row < width && phase*TILE_WIDTH+ty < width)
			submatrix_A[tx][ty] = A[row*width + phase*TILE_WIDTH+ty];
		else 
			submatrix_A[tx][ty] = 0;

		if (col < width &&  phase*TILE_WIDTH +tx < width)
			submatrix_B[tx][ty] = B[(phase*TILE_WIDTH +tx)*width + col]; 
		else
			submatrix_B[tx][ty] = 0;
	
		__syncthreads();

		// calculation within a phase
		if (row < width && col < width)
			for(int i = 0; i<TILE_WIDTH; i++)
				temp += submatrix_A[tx][i] * submatrix_B[i][ty];

		__syncthreads();
	}
	//write back values
	if (row<width && col <width)
		C[row*width + col] = temp;
}

int main(int argc, char** argv){
	int width, i;
	
	if (argc!=2 || !(width=atoi(argv[1])) || width<=0){
		fprintf(stderr,"Please enter a valid width!\n");
		exit(-1);
	}

	//allocating space
	u_int32_t num_elements = width*width;
	u_int32_t matrix_size = sizeof(float)*num_elements;
	float* temp=(float*)malloc(matrix_size);
	float *d_A, *d_B, *d_C;
	cudaMalloc(&d_A, matrix_size);
	cudaMalloc(&d_B, matrix_size);
	cudaMalloc(&d_C, matrix_size);

	srand(0);
	//create and transfer random numbers
	for(i=0;i<num_elements;i++)
		temp[i] = 20.*rand()/RAND_MAX-10;
	cudaMemcpy(d_A, temp, matrix_size, cudaMemcpyHostToDevice);
	
	for(i=0;i<num_elements;i++)
		temp[i] = 20.*rand()/RAND_MAX-10;
	cudaMemcpy(d_B, temp, matrix_size, cudaMemcpyHostToDevice);
	
	int block_len = ceil(1.*width/TILE_WIDTH);
	dim3 dimGrid(block_len, block_len);
	dim3 dimBlock(TILE_WIDTH, TILE_WIDTH);

	//warm_up
	matmul<<<dimGrid,dimBlock>>>(d_A, d_B, d_C, width);
	cudaDeviceSynchronize();

	//initializing timing tools
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//timing
	cudaEventRecord(start);
	matmul<<<dimGrid,dimBlock>>>(d_A, d_B, d_C, width);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	//Transfer back and clean up
	cudaMemcpy(temp, d_C, matrix_size, cudaMemcpyDeviceToHost);
	float msec;
	cudaEventElapsedTime(&msec, start, stop);
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
	
	
	// outputing
	
	FILE* fp = fopen("product.dat","w+");
	if(fp == NULL){
		fprintf(stderr, "Failed to load output file");
		exit(-1);
	}
	int j;
	for(i=0;i<width;i++){
		for(j=0;j<width-1;j++)
			fprintf(fp,"%-5.2f\t",temp[i * width + j]);
		fprintf(fp,"%-5.2f\n",temp[i * width + width -1]);
	}
	//close files and free electrons
	fclose(fp);
	free(temp);
	printf("Computation finished in %.5f seconds with N = %d \n", msec/1000, width);
	
}
