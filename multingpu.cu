#include <iostream>
#include <stdlib.h>
#include <set>
#include <fstream>

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <math.h>
#include <random>
#include <chrono>
#include <ratio>
#include <thread>
#include <mutex>

//#define MODULUS_PRIME  1073741827// 30 bit prime
#define MODULUS_PRIME 536870909 //29 bit prime
//#define MODULUS_PRIME 9973 //14 bit prime
//#define MODULUS_PRIME 11 //4 bit prime

static const int numberOfThreads = 250;
static const int NumberOfTrials = 1;

void generateQVector(unsigned *queryVector, int p);
void printMatrix(int **a,int r, int c);
void printVector(unsigned *a,int c);
void printVector2(unsigned long long int *a,int c);
int checkIfEqual(unsigned long long int  *resultSparseVectorValueArray, unsigned long long int  *resultSparseVectorValueArrayDD, int length);
int checkIfEqual2(unsigned  *resultSparseVectorValueArray, unsigned  *resultSparseVectorValueArrayDD, int length);
__global__ void	dvsmMultKernel(int lengthOfResultVectorReduced, int *columnPtr,int *rowIndex,unsigned *valueArray,
		unsigned *queryVector,unsigned long long int *resultSparseVectorValueArray);

int main()
{
	std::ofstream myfile;
	myfile.open ("TestCase-CUDA100trialsR20.txt");
	uint64_t seed = time(NULL);
	myfile << "Seed: " << seed ;
//	std::cout << "Seed: " << seed ;
	srand(seed);
	const long max_u =16, r = 1L << 20;
	for (long p = 2; p <= r; p <<= 1)
	{
		for (long u = 1; u <= max_u; u <<=1)
		{
			//			std::cout << "************************************************************************\n";
//			std::cout << "\n\np: " << p << "\tr: " << r << "\tu: " << u ;

			//			myfile << "************************************************************************\n";
			myfile << "\n\np: " << p << "\tr: " << r << "\tu: " << u ;

			for(int trials=0;trials<NumberOfTrials;trials++){

				//top:


				std::chrono::high_resolution_clock::time_point pruTotalTimeStart = std::chrono::high_resolution_clock::now();
				
				int numberOfNonZeroElements = p*u;
				int lengthOfColumnPtr = r+1;
				unsigned *valueArray = (unsigned*)malloc(sizeof(unsigned)*numberOfNonZeroElements);
				int *rowIndex = (int*)malloc(sizeof(int)*numberOfNonZeroElements);
				int *columnPtr = (int*)malloc(sizeof(int)*(lengthOfColumnPtr));
				//std::cout << "\nval (" << numberOfNonZeroElements << "): ";
				for (long i = 0; i < p * u; i++) {
					valueArray[i] = (rand() % MODULUS_PRIME);
					//std::cout << valueArray[i] << ",";
				}
				//std::cout << "\n\nRow: ";
				int t=0;
				int sum=0;
				columnPtr[0] = 0;
				int lengthOfCPReduced = 0;
				for (int i = 0; i < r; ++i)
				{
					for (auto it = begin(cols[i]); it != end(cols[i]); ++it)
					{
						rowIndex[t++] = (*it)-1;
						//std::cout << rowIndex[t-1] << ",";
					}
					if (cols[i].size())
					{
						columnPtr[lengthOfCPReduced+1]=columnPtr[lengthOfCPReduced]+cols[i].size();
						lengthOfCPReduced++;
					}
					sum+=cols[i].size();
				}

				//std::cout << "\n\nCol (" << cols.size() <<"): ";
				std::chrono::high_resolution_clock::time_point pruMatrixCreationEnd = std::chrono::high_resolution_clock::now();
				std::chrono::duration<int,std::milli> pruMatrixCreationDuration = std::chrono::duration_cast<std::chrono::duration<int,std::milli>>(pruMatrixCreationEnd-pruTotalTimeStart);
//				std::cout<<"\nTrial: "<<trials+1<< "\tMatrix Creation: "<< (double) pruMatrixCreationDuration.count()/1000 <<"s";
				myfile <<"\nTrial: "<<trials+1<< "\tMatrix Creation: "<< (double) pruMatrixCreationDuration.count()/1000 <<"s";
				/*
				 * CUDA started
				 *
				 **/
				std::chrono::high_resolution_clock::time_point pruMatrixCopyStart = std::chrono::high_resolution_clock::now();
				int lengthOfResultVectorReduced = lengthOfCPReduced-1;
				int THREADSPERBLOCKDVSM = lengthOfResultVectorReduced < 1024 ? lengthOfResultVectorReduced : 1024;
				int NDVSM = (lengthOfResultVectorReduced+THREADSPERBLOCKDVSM-1) / THREADSPERBLOCKDVSM;

				//			unsigned long long int   *resultSparseVectorValueArray = (unsigned long long int *)malloc(sizeof(unsigned long long int)*lengthOfResultVectorReduced*numberOfThreads);
				unsigned long long int   *resultSparseVectorValueArrayDD;
				unsigned long long int   *resultSparseVectorValueArray_d;
				unsigned *queryVector;

				int *rowIndex_d, *columnPtr_d;
				unsigned *valueArray_d, *queryVector_d;


				cudaMalloc((void**)&valueArray_d,(numberOfNonZeroElements*sizeof(unsigned)));
				cudaMalloc((void**)&rowIndex_d,(numberOfNonZeroElements*sizeof(int)));
				cudaMalloc((void**)&columnPtr_d,(lengthOfCPReduced)*sizeof(int));
				cudaMalloc((void**)&queryVector_d,numberOfThreads*p*sizeof(unsigned));
				cudaMalloc((void**)&resultSparseVectorValueArray_d,(numberOfThreads*lengthOfResultVectorReduced*sizeof(unsigned long long int)));

				cudaMemcpy( valueArray_d, valueArray, numberOfNonZeroElements*sizeof(unsigned), cudaMemcpyHostToDevice );
				cudaMemcpy( rowIndex_d, rowIndex, numberOfNonZeroElements*sizeof(int), cudaMemcpyHostToDevice );
				cudaMemcpy( columnPtr_d, columnPtr, lengthOfCPReduced*sizeof(int), cudaMemcpyHostToDevice );

				cudaStream_t *streams = new cudaStream_t[numberOfThreads];

				for (int i = 0; i < numberOfThreads; ++i) cudaStreamCreate(&streams[i]);
				cudaMallocHost((void**)&queryVector,numberOfThreads*p*sizeof(unsigned));
				cudaMallocHost((void**)&resultSparseVectorValueArrayDD,numberOfThreads*lengthOfResultVectorReduced*sizeof(unsigned long long int));
				generateQVector(queryVector,p*numberOfThreads);

				int i = 0;
				int numberOfOps=0;
				std::chrono::duration<int,std::nano> timeSpend;
				std::chrono::nanoseconds zeroSec{0};
				std::chrono::nanoseconds nsInOneSec{1000000000};
				timeSpend = zeroSec;
				while(timeSpend < nsInOneSec)
				{
					std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
					cudaMemcpyAsync(queryVector_d+i*p, queryVector+i*p, p*sizeof(unsigned), cudaMemcpyHostToDevice, streams[i]);
					dvsmMultKernel<<<NDVSM,THREADSPERBLOCKDVSM,0,streams[i]>>>(lengthOfResultVectorReduced,columnPtr_d, rowIndex_d, valueArray_d, queryVector_d+i*p, resultSparseVectorValueArray_d+i*p);
					cudaMemcpyAsync(resultSparseVectorValueArrayDD+i*p, resultSparseVectorValueArray_d+i*p, lengthOfResultVectorReduced*sizeof(unsigned long long int), cudaMemcpyDeviceToHost, streams[i]);
					std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
					std::chrono::duration<int,std::nano> time_span = std::chrono::duration_cast<std::chrono::duration<int,std::nano>>(t2 - t1);
					timeSpend = timeSpend + time_span ;

					i = (i + 1) % numberOfThreads;
					numberOfOps++;

				}



//				std::cout<< "\t\tNumberOfOperations: "<< numberOfOps -1 ;
				myfile<< "\t\tNumberOfOperations: "<< numberOfOps -1 ;

				cudaFree( valueArray_d );
				cudaFree( rowIndex_d );
				cudaFree( columnPtr_d );
				cudaFree( queryVector_d );
				cudaFree( resultSparseVectorValueArray_d );

				cudaFreeHost(queryVector);
				cudaFreeHost(resultSparseVectorValueArrayDD);
				free(valueArray);
				free(rowIndex);
				free(columnPtr);
				//			free(resultSparseVectorValueArray);
				for (int i = 0; i < numberOfThreads; ++i) cudaStreamDestroy(streams[i]);
				std::chrono::high_resolution_clock::time_point pruTotalTimeEnd = std::chrono::high_resolution_clock::now();
				std::chrono::duration<int,std::milli> pruTotalTimeDuration = std::chrono::duration_cast<std::chrono::duration<int,std::milli>>(pruTotalTimeEnd-pruTotalTimeStart);

				std::chrono::duration<int,std::milli> pruTotalMatrixCopyDuration = std::chrono::duration_cast<std::chrono::duration<int,std::milli>>(pruTotalTimeEnd-pruMatrixCopyStart);
//				std::cout<< "\tMatrixCopy: "<< ((double) pruTotalMatrixCopyDuration.count()/1000) - 1.0 <<"s";
				myfile<< "\tMatrixCopy: "<< ((double) pruTotalMatrixCopyDuration.count()/1000) - 1.0 <<"s";


//				std::cout<< "\tTotalTime: "<< (double) pruTotalTimeDuration.count()/1000 <<"s";
				myfile<< "\tTotalTime: "<< (double) pruTotalTimeDuration.count()/1000 <<"s";
			}

		}
	}
	myfile.close();

}

void generateQVector(unsigned *queryVector, int p){
	int i;
	for(i=0;i<p;i++){
		queryVector[i] = rand() % MODULUS_PRIME + 1;
	}
}


void printMatrix(int **a,int r, int c) {
	int i=0,j=0;
	for(;i<r;i++){
		for(j=0;j<c;j++){
			printf("%d ",a[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void printVector(unsigned *a,int c) {
	int j=0;
	for(j=0;j<c;j++){
		printf("%d ",a[j]);
	}
	printf("\n");
}

void printVector2(unsigned long long int *a,int c) {
	int j=0;
	for(j=0;j<c;j++){
		printf("%d ",a[j]);
	}
	printf("\n");
}


int checkIfEqual(unsigned long long int  *resultSparseVectorValueArray, unsigned long long int  *resultSparseVectorValueArrayDD, int length) {
	int i;
	for(i=0;i<length;i++){
		if(resultSparseVectorValueArray[i]!=resultSparseVectorValueArrayDD[i]) {
			return i;
		}
	}
	return 1;
}

int checkIfEqual2(unsigned  *resultSparseVectorValueArray, unsigned  *resultSparseVectorValueArrayDD, int length) {
	int i;
	for(i=0;i<length;i++){
		if(resultSparseVectorValueArray[i]!=resultSparseVectorValueArrayDD[i]) {
			return i;
		}
	}
	return 1;
}

__global__ void	dvsmMultKernel(int lengthOfResultVectorReduced, int *columnPtr,int *rowIndex,unsigned *valueArray,
		unsigned *queryVector,unsigned long long int *resultSparseVectorValueArray){
	int col = blockDim.x * blockIdx.x + threadIdx.x;
	if(col<lengthOfResultVectorReduced){

		unsigned long long int temp = 0;
		int j;
		for(j=columnPtr[col];j<columnPtr[col+1];j++) {
			temp += valueArray[j]*queryVector[rowIndex[j]];
		}
		resultSparseVectorValueArray[col]= temp % MODULUS_PRIME; //mul_m(1,temp,MODULUS_PRIME,INVK);//
	}

}




/*printf("GPU (Device) result sparse vector: Length %d\n",lengthOfResultVectorReduced);
				printVector2(resultSparseVectorValueArrayDD+i*p,lengthOfResultVectorReduced%20);
				for(int w=0;w<lengthOfResultVectorReduced;w++) {
					unsigned long long int  temp = 0;
					for(int j=columnPtr[i];j<columnPtr[w+1];j++) {
						temp += valueArray[j]*(queryVector+i*p)[rowIndex[j]];
					}
					(resultSparseVectorValueArray+i*p)[w]=temp % MODULUS_PRIME;
				}

				printf("CPU (Host) result sparse vector: Length %d\n",lengthOfResultVectorReduced);
				printVector2(resultSparseVectorValueArray+i*p,lengthOfResultVectorReduced%20);

				if(checkIfEqual(resultSparseVectorValueArray+i*p, resultSparseVectorValueArrayDD+i*p, lengthOfResultVectorReduced)){
					printf("Vectors are matched.\n");
				} else {
					printf("Vectors are not matched.\n");
				}*/
