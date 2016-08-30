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

void runKernel(unsigned *queryVector,int p, unsigned *queryVector_d, int NDVSM, int THREADSPERBLOCKDVSM, int lengthOfResultVectorReduced,int *columnPtr_d,int *rowIndex_d,unsigned *valueArray_d,unsigned long long int *resultSparseVectorValueArray_d, unsigned long long int *resultSparseVectorValueArrayDD);

static const int numberOfThreads = 500;
std::mutex mtxKernel;

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

int main()
{
	std::ofstream myfile;
	myfile.open ("testDataStatNOPsE.txt");
	srand (time(NULL));
	const long max_u = 16, r = 1L << 18;
	for (long p = 2; p <= r; p <<= 1)
	{
		for (long u = 1; u <= max_u; ++u)
		{
			//top:
			std::cout << "************************************************************************\n";
			std::cout << "p: " << p << "; r: " << r << "; u: " << u << "\n";

			myfile << "************************************************************************\n";
			myfile << "p: " << p << "; r: " << r << "; u: " << u << "\n";

			//			long k = 0;
			std::vector<std::set<long>> cols;
			std::vector<std::set<long>*> cols2;
			cols.resize(r);
			for (auto it = begin(cols); it != end(cols); ++it) cols2.push_back(&(*it));

			for (long i = 1; i <= p; ++i)
			{
				for (long j = 1; j <= u; )
				{
					long c = rand() % cols2.size();
					if (cols2[c]->size() < u && cols2[c]->insert(i).second)
					{
						j++;
					}
					else
					{
						long a = rand() % r;
						if (cols[a].size() > 0 && cols[a].find(i) == end(cols[a]))
						{
							auto elt = begin(cols[a]);
							std::advance(elt, rand() % cols[a].size());
							long tmp = *elt;
							if (cols2[c]->find(tmp) == end(*(cols2[c])))
							{
								cols[a].erase(elt);
								cols[a].insert(i);
								cols2[c]->insert(tmp);
								j++;
							}
						}
					}
					if (cols2[c]->size() == u) cols2.erase(begin(cols2) + c);
				}
			}

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


			/*
			 * CUDA started
			 *
			 **/
			int lengthOfResultVectorReduced = lengthOfCPReduced-1;
			int THREADSPERBLOCKDVSM = lengthOfResultVectorReduced < 1024 ? lengthOfResultVectorReduced : 1024;
			int NDVSM = (lengthOfResultVectorReduced+THREADSPERBLOCKDVSM-1) / THREADSPERBLOCKDVSM;


			unsigned long long int   *resultSparseVectorValueArrayDD = (unsigned long long int  *)malloc(sizeof(unsigned long long int)*lengthOfResultVectorReduced*numberOfThreads);
			unsigned long long int   *resultSparseVectorValueArray_d;

			unsigned *queryVector =  (unsigned*)malloc(sizeof(unsigned)*p*numberOfThreads);

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

			unsigned long numberOfOps;
			generateQVector(queryVector,p*numberOfThreads);
			std::thread thrds[numberOfThreads];
			std::chrono::duration<int,std::nano> timeSpend;
			std::chrono::nanoseconds zeroSec{0};
			timeSpend = zeroSec;
			//std::chrono::nanoseconds nsInOneSec{1000000000};
			std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
			for(numberOfOps=0; numberOfOps < numberOfThreads; numberOfOps++){

				runKernel(queryVector+numberOfOps,p,queryVector_d+numberOfOps,NDVSM,THREADSPERBLOCKDVSM,lengthOfResultVectorReduced,columnPtr_d,rowIndex_d,valueArray_d,resultSparseVectorValueArray_d+numberOfOps,resultSparseVectorValueArrayDD+numberOfOps);

			}
			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<int,std::nano> time_span = std::chrono::duration_cast<std::chrono::duration<int,std::nano>>(t2 - t1);


			timeSpend = timeSpend + time_span ;

			std::cout << "Number of "<< numberOfThreads <<" operations at GPU takes: "<< (double)timeSpend.count()/1000 << "ms\n";
			myfile << "Number of "<< numberOfThreads <<" operations at GPU takes: "<< (double)timeSpend.count()/1000 << "ms\n";

			cudaFree( valueArray_d );
			cudaFree( rowIndex_d );
			cudaFree( columnPtr_d );
			cudaFree( queryVector_d );
			cudaFree( resultSparseVectorValueArray_d );

			free(queryVector);
			free(resultSparseVectorValueArrayDD);
			free(valueArray);
			free(rowIndex);
			free(columnPtr);


			/*
			 * CUDA finished
			 *
			 **/	}
	}
	myfile.close();

}
void runKernel(unsigned *queryVector,int p, unsigned *queryVector_d, int NDVSM, int THREADSPERBLOCKDVSM, int lengthOfResultVectorReduced,int *columnPtr_d,int *rowIndex_d,unsigned *valueArray_d,unsigned long long int *resultSparseVectorValueArray_d, unsigned long long int *resultSparseVectorValueArrayDD) {
	cudaMemcpy( queryVector_d, queryVector, p*sizeof(unsigned), cudaMemcpyHostToDevice );
	//	mtxKernel.lock();
	dvsmMultKernel<<< NDVSM, THREADSPERBLOCKDVSM>>>(lengthOfResultVectorReduced,columnPtr_d, rowIndex_d, valueArray_d, queryVector_d, resultSparseVectorValueArray_d);
	//	mtxKernel.unlock();
	cudaMemcpy( resultSparseVectorValueArrayDD, resultSparseVectorValueArray_d, (lengthOfResultVectorReduced*sizeof(unsigned long long int)), cudaMemcpyDeviceToHost );
}
