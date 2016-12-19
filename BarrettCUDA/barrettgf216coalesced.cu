// This file is part of BarretCUDA v0.1 
// 
// BarretCUDA is a fast(ish) CUDA implementation of sparse matrix
// multiplication modulo a multi-precision prime.
// 
// Copyright (C) 2016, Ryan Henry and Syed Mahbub Hafiz
// 
// 
// BarretCUDA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// 
// BarretCUDA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with BarretCUDA.  If not, see <http://www.gnu.org/licenses/>.
#include <iostream>
#include <cuda.h>
#include <atomic>
#include <chrono>
#include <ratio>
#include <fstream>
#include "gf2earrays.h"
#include "barret.h"

using namespace std;

#define THREADS_PER_BLOCK(n) 	(n >= 512 ? 512 : n)
#define NUM_BLOCKS(n)	((n + THREADS_PER_BLOCK(n) - 1) / THREADS_PER_BLOCK(n))

#define DEBUG_IDX 4

typedef uint16_t uintX;

int u;

__device__ __forceinline__ uint16_t multiply_GF216_Element(const uint16_t x, const uint16_t y,  const uint16_t * d_GF216_log_table, const uint16_t * d_GF216_exp_table)
{
    if(x == 0 || y == 0) return 0;
    return d_GF216_exp_table[d_GF216_log_table[x]+d_GF216_log_table[y]];
}

 void Check_CUDA_Error(const char *message)
{
    cudaError_t error = cudaGetLastError();
    if(error!=cudaSuccess) {
       fprintf(stderr,"ERROR: %s: %s\n", message, cudaGetErrorString(error) );
       exit(-1);
    }                         
}

void SpMV_gf216_cpu(int u_wrap, uint16_t * response, const uint16_t * query, const SparseMatrix<uint> & matrix)
{
	for (int i = 0; i < matrix.ncols; i++)
	{
	     uint col = matrix.l_cols[i];
	     uint16_t temp = 0;
		
	    for (int j = 0; j < u_wrap; j++)
	    {   
		const uint val = matrix.l_vals[i + j * matrix.ncols];
		//if(i==101) printf("CPU: val: %u, col: %u, rows: %u, query: %u, temp: %u\n",(val      )& 0xFFFF,col,matrix.l_rows[col],query[matrix.l_rows[col]],multiply_GF2E<uint16_t>(((val      )& 0xFFFF), query[matrix.l_rows[col]]));
		temp ^= multiply_GF2E<uint16_t>(((val      )& 0xFFFF), query[matrix.l_rows[col++]]);
		//if(i==101) printf("CPU: val: %u, col: %u, rows: %u, query: %u, temp: %u\n",(val >> 16)& 0xFFFF,col,matrix.l_rows[col],query[matrix.l_rows[col]],multiply_GF2E<uint16_t>(((val  >> 16)& 0xFFFF), query[matrix.l_rows[col]]));	
		temp ^= multiply_GF2E<uint16_t>(((val >> 16)& 0xFFFF), query[matrix.l_rows[col++]]);
			
	    }
	    response[i] = temp;
	if(i==101)		printf("CPU response[%d]: %u\n",i,response[i]);
	}
}

__global__ void SpMV_kernel(int u_wrap, uint16_t * response, const uint16_t * query, const uint nvals,
	const uint * vals, const uint ncols, const uint * cols, const uint * rows,  const uint16_t * d_GF216_log_table, const uint16_t * d_GF216_exp_table)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= ncols) return;
    register uint col = cols[i];
    register uint16_t temp = 0;
    	for (int j = 0; j < u_wrap; j++)
    	{
		register const uint val = vals[i + j * ncols];
		//if(i==101) printf("Kernel: val: %u, col: %u, rows: %u, query: %u, temp: %u\n",(val      )& 0xFFFF,col,rows[col],query[rows[col]],multiply_GF216_Element(val      & 0xFFFF, query[rows[col]],d_GF216_log_table,d_GF216_exp_table));
		temp ^= multiply_GF216_Element(((val      )& 0xFFFF), query[rows[col++]],d_GF216_log_table,d_GF216_exp_table);
		//if(i==101) printf("Kernel: val: %u, col: %u, rows: %u, query: %u, temp: %u\n",(val >> 16)& 0xFFFF,col,rows[col],query[rows[col]],multiply_GF216_Element(((val  >> 16)& 0xFFFF), query[rows[col]],d_GF216_log_table,d_GF216_exp_table));	
		temp ^= multiply_GF216_Element(((val >> 16)& 0xFFFF), query[rows[col++]],d_GF216_log_table,d_GF216_exp_table);
	}
    response[i] = temp;
	if(i == 101) printf("GPU response[%d]: %u\n",i,response[i]);
}

void SpMV(int u_wrap, uint16_t * l_response, const uint16_t * l_query,
	uint16_t * d_response, uint16_t * d_query, const cudaStream_t & stream,
	const SparseMatrix<uint> & matrix, const uint16_t * d_GF216_log_table, const uint16_t * d_GF216_exp_table)
{
    cudaMemcpyAsync(d_query, l_query, matrix.nrows*sizeof(uint16_t), cudaMemcpyHostToDevice, stream);
    Check_CUDA_Error("GF216:cudaMemcpyAsync:cudaMemcpyHostToDevice");
    const dim3 Dg(NUM_BLOCKS(matrix.ncols), 1, 1);
    const dim3 Db(THREADS_PER_BLOCK(matrix.ncols), 1, 1);
    const size_t Ns = 0;

    SpMV_kernel<<< Dg, Db, Ns, stream >>> (u_wrap, d_response, d_query, matrix.nvals, matrix.d_vals, matrix.ncols, matrix.d_cols, matrix.d_rows, d_GF216_log_table,d_GF216_exp_table);
    Check_CUDA_Error("GF216:SpMV_kernel");
    cudaMemcpyAsync(l_response, d_response, matrix.ncols*sizeof(uint16_t), cudaMemcpyDeviceToHost, stream);
    Check_CUDA_Error("GF216:cudaMemcpyAsync:cudaMemcpyDeviceToHost");
}
template <typename T>
void initMatrix(const char * valfile, const char * rowfile,
	const char * colfile, NTL::ZZ & modulus,
	struct SparseMatrix<T> & matrix)
{

    std::ifstream rowstream(rowfile, std::ifstream::in);
    if (!rowstream) { cerr << "Error: opening ROWS files\n"; exit(-1); }
    std::ifstream colstream(colfile, std::ifstream::in);
    if (!colstream) { cerr << "Error: opening COLS files\n"; exit(-1); }
    std::ifstream valstream(valfile, std::ifstream::in);
    if (!valstream) { cerr << "Error: opening VALS files\n"; exit(-1); }

    	NTL::ZZ temp_zz;
	valstream >> temp_zz;
	modulus = NTL::trunc_ZZ(temp_zz, 16);

    //valstream >> u;
    rowstream >> matrix.nrows;
    rowstream >> matrix.nvals;
    colstream >> matrix.ncols;
    matrix.l_cols = (uint *)malloc((matrix.ncols+1) * sizeof(uint));
    cudaMalloc((void**)&matrix.d_cols, (matrix.ncols+1) * sizeof(uint));

    int u_wrap = ((u-1)/2+1);
  //std::cout << "nvals: " << matrix.nvals << " Size l_vals: "<< u_wrap*matrix.ncols* sizeof(uint) <<" u_wrap really: "<< u_wrap <<"\n";
    matrix.l_rows = (uint *)malloc(matrix.ncols * u * sizeof(uint));memset(matrix.l_rows+matrix.nvals,0,(matrix.ncols * u-matrix.nvals) * sizeof(uint));
    cudaMalloc((void**)&matrix.d_rows, matrix.ncols * u * sizeof(uint));
    matrix.l_vals = (uint *)malloc(u_wrap*matrix.ncols* sizeof(uint)); memset(matrix.l_vals,0,u_wrap*matrix.ncols* sizeof(uint));
    cudaMalloc((void**)&matrix.d_vals, u_wrap*matrix.ncols* sizeof(uint));



//std::cout << "modulus:\t" << modulus << "\n";
//std::cout << "matrix.nrows:\t" << matrix.nrows << "\n";
//std::cout << "matrix.ncols:\t" << matrix.ncols << "\n";
//std::cout << "matrix.nvals:\t" << matrix.nvals << "\n";

    for (int i = 0; i < matrix.ncols+1; i++)
    {
	colstream >> matrix.l_cols[i];
    }
    colstream.close();


    cudaMemcpy(matrix.d_cols, matrix.l_cols, (matrix.ncols+1) * sizeof(uint),
	cudaMemcpyHostToDevice);



    NTL::ZZ_pPush p(modulus);
    for (int i = 0; i < matrix.nvals; i++)
    {
	rowstream >> matrix.l_rows[i];
    }
    rowstream.close();
    /*for (int i = matrix.nvals; i < matrix.ncols * u; i++ ) 
    {
	matrix.l_rows[i] = 0;
    }*/
    cudaMemcpy(matrix.d_rows, matrix.l_rows, matrix.nvals* sizeof(uint),
	cudaMemcpyHostToDevice);

	for(int i=0,c=0;i<matrix.ncols;i++)
	{	int number_of_nnz_in_a_col = matrix.l_cols[i+1] - matrix.l_cols[i];
		//std::cout<< "col i: " << number_of_nnz_in_a_col;
		for(int j=0; j < number_of_nnz_in_a_col;j++)
		{	c++;
			
			NTL::ZZ tmp_zz;
			valstream >> tmp_zz;
			uint tmp = NTL::trunc_long(tmp_zz, 16);
			matrix.l_vals[i+(j/2)*matrix.ncols] +=  tmp<<(16*(j%2));		
//			std::cout << "\nc: " << c << " at ("<< i+(j/2)*matrix.ncols << ", " << 8*(j%2) << ")";
		}
		//std::cout << "\n";

	}
	//to_uint<T>(NTL::rep(tmp), matrix.l_vals[i]);
    valstream.close();
  //std::cout << "problem9\n\n";
    cudaMemcpy(matrix.d_vals, matrix.l_vals, u_wrap*matrix.ncols* sizeof(uint),
	cudaMemcpyHostToDevice);
//  std::cout << "problem10\n\n";


}

template <typename T>
void freeMatrix(struct SparseMatrix<T> & matrix)
{
    free(matrix.l_vals);
    free(matrix.l_rows);
    free(matrix.l_cols);
    cudaFree(matrix.d_vals);
    cudaFree(matrix.d_cols);
    cudaFree(matrix.d_rows);
}
int main(int argc, char ** argv)
{
    int nstreams = 2;

    u = atoi(argv[4]);

    int u_wrap = ((u-1)/2+1);
    if (argc < 3)
    {
	std::cout << "Usage: " << argv[0] << " VALUES ROWS COLS\n\n";
	return 1;
    }

    time_t t0 = time(0);
    NTL::SetSeed(to_ZZ(t0));
    std::cout << "Seed: " << t0 << "\n";

    struct SparseMatrix<uint> matrix = { 0 };
    NTL::ZZ modulus;

    initMatrix(argv[1], argv[2], argv[3], modulus, matrix);
    NTL::ZZ_p::init(modulus);

    uint16_t * d_GF216_log_table;//65536
    uint16_t * d_GF216_exp_table;//131070

    cudaMalloc((void**)&d_GF216_log_table, sizeof(GF216_log_table));
    cudaMemcpy(d_GF216_log_table, &GF216_log_table, sizeof(GF216_log_table),cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_GF216_exp_table, sizeof(GF216_exp_table));
    cudaMemcpy(d_GF216_exp_table, &GF216_exp_table, sizeof(GF216_exp_table),	cudaMemcpyHostToDevice);

    uintX * l_query, * d_query;
    cudaMallocHost((void**)&l_query, nstreams * matrix.nrows * sizeof(uintX));
    cudaMalloc((void**)&d_query, nstreams * matrix.nrows * sizeof(uintX));

    uint16_t * l_response, * d_response;
    cudaMallocHost((void**)&l_response, nstreams * matrix.ncols * sizeof(uint16_t));
    cudaMalloc((void**)&d_response, nstreams * matrix.ncols * sizeof(uint16_t));
    uint16_t * cpu_l_response = (uint16_t*) malloc(nstreams*matrix.ncols * sizeof(uint16_t));

    cudaStream_t * streams = new cudaStream_t[nstreams];
    for (int i = 0; i < nstreams; ++i) cudaStreamCreate(&streams[i]);

    for (int i = 0; i < nstreams * matrix.nrows; i++)
    {
	l_query[i] = NTL::RandomBits_long(16);
    }

    std::atomic<int> cnt = ATOMIC_VAR_INIT(0);
    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds onesec{1000000000};
    //while (std::chrono::duration_cast<std::chrono::duration<int,std::nano>>(std::chrono::high_resolution_clock::now() - start) < onesec)
    {
	//int i = cnt % nstreams;
	#pragma omp parallel
	for (int i = 0; i < nstreams; i++)
	{
	    uint16_t * __l_response = l_response + i * matrix.ncols;
	    uint16_t * __d_response = d_response + i * matrix.ncols;
	    uintX * __l_query = l_query + i * matrix.nrows;
	    uintX * __d_query = d_query + i * matrix.nrows;		

	    SpMV(u_wrap,__l_response, __l_query, __d_response,
		__d_query, streams[i], matrix, d_GF216_log_table, d_GF216_exp_table);
	    SpMV_gf216_cpu(u_wrap, cpu_l_response, __l_query, matrix);


	    std::atomic_fetch_add(&cnt, 1);
	    //for (int j = 0; j < matrix.nrows; j++)
	    //{	__l_query[j] = 65535;
		//to_uint<uintX>(NTL::rep(NTL::random_ZZ_p()), __l_query[j]);
		//NTL::BytesFromZZ((unsigned char *)&__l_query[j], NTL::RandomPrime_ZZ(8), 1);//conv<uint>();
	    //}
	}
    }

    std::cout << "Count: " << cnt << "\n";

    // cleanup
    for (int i = 0; i < nstreams; ++i) cudaStreamDestroy(streams[i]);
    delete [] streams;
    cudaFreeHost(l_query);
    cudaFree(d_query);
    cudaFreeHost(l_response);
    cudaFree(d_response);
    cudaFree(d_GF216_log_table);
    cudaFree(d_GF216_exp_table);
    freeMatrix<uint>(matrix);
    return 0;
}



