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

#include <atomic>
#include <fstream>

#include "barret.h"
#include "uintX.h"

#define NUM_BLOCKS 		256
#define THREADS_PER_BLOCK(n)	(n / NUM_BLOCKS)

NTL_CLIENT

template <typename T>
__global__ void SpMV_kernel(T * response, const T * query, const uint nvals,
	const T * vals, const uint ncols, const uint * cols, const uint * rows,
	const T * modulus, const T * mu, const T * subtrahends)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    response[i] = { 0 };
    T hi = { 0 };
    uint overflow = { 0 };
    for (int j = cols[i]; j < cols[i+1]; ++j)
    {
	mad(response[i], hi, overflow, vals[j], query[rows[j]]);
    }
    normalize(response[i], hi, subtrahends + 2*overflow);
    uintXp<T> q = get_q(response[i], hi, *mu);
    uintXp<T> r2 = get_r2(q, *modulus);
    sub(response[i], hi, r2);
}

template <typename T>
void SpMV(NTL::vec_ZZ_p & response, T * l_response, const T * l_query,
	T * d_response,	T * d_query, const cudaStream_t & stream,
	const SparseMatrix<T> & matrix, const BarretParams<T> & barret)
{
    cudaMemcpyAsync(d_query, l_query, matrix.nrows * sizeof(T),
	cudaMemcpyHostToDevice, stream);
    
    const dim3 Dg(NUM_BLOCKS, 1, 1);
    const dim3 Db(THREADS_PER_BLOCK(matrix.ncols), 1, 1);
    const size_t Ns = 0;

    SpMV_kernel<T> <<< Dg, Db, Ns, stream >>> (d_response, d_query,
	matrix.nvals, matrix.vals, matrix.ncols, matrix.cols, matrix.rows,
	barret.modulus, barret.mu, barret.subtrahends);

    cudaMemcpyAsync(l_response, d_response, matrix.ncols * sizeof(T),
	cudaMemcpyDeviceToHost, stream);
    
    response.SetLength(matrix.ncols);
    for (int i = 0; i < matrix.ncols; ++i)
    {
	response[i] = to_ZZ_p<T>(l_response[i]);
    }
}

int main(int argc, char ** argv)
{
    int nstreams = 1;

    if (argc < 3)
    {
	cout << "Usage: " << argv[0] << " VALUES ROWS COLS\n\n";
	return 1;
    }

    struct SparseMatrix<uintX> matrix;
    NTL::ZZ modulus;
    initMatrix(argv[1], argv[2], argv[3], modulus, matrix);
    NTL::ZZ_p::init(modulus);

    struct BarretParams<uintX> barret;
    initBarret<uintX>(modulus, barret);

    uintX * l_query, * d_query;
    cudaMallocHost((void**)&l_query, nstreams * matrix.nrows * sizeof(uintX));
    cudaMalloc((void**)&d_query, nstreams * matrix.nrows * sizeof(uintX));
    uintX * l_response, * d_response;
    cudaMallocHost((void**)&l_response, nstreams * matrix.ncols * sizeof(uintX));
    cudaMalloc((void**)&d_response, nstreams * matrix.ncols * sizeof(uintX));
    
    cudaStream_t * streams = new cudaStream_t[nstreams];
    for (int i = 0; i < nstreams; ++i) cudaStreamCreate(&streams[i]);

    NTL::vec_vec_ZZ_p responses(INIT_SIZE, nstreams,
	NTL::vec_ZZ_p(INIT_SIZE, matrix.ncols));
    for (int i = 0; i < nstreams * matrix.nrows; i++)
    {
	to_uint<uintX>(NTL::rep(NTL::random_ZZ_p()), l_query[i]);
    }

    std::atomic<int> cnt = ATOMIC_VAR_INIT(0);
    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds onesec{1000000000};
    while (std::chrono::duration_cast<std::chrono::duration<int,std::nano>>(std::chrono::high_resolution_clock::now() - start) < onesec)
    {
    	int i = cnt % nstreams;
    	uintX * __l_response = l_response + i * matrix.ncols;
    	uintX * __d_response = d_response + i * matrix.ncols;
    	uintX * __l_query    = l_query    + i * matrix.nrows;
    	uintX * __d_query    = d_query    + i * matrix.nrows;

	SpMV<uintX>(responses[i], __l_response, __l_query, __d_response,
	    __d_query, streams[i], matrix, barret);
	std::atomic_fetch_add(&cnt, 1);
	
	for (int j = 0; j < matrix.nrows; j++)
	{
	    to_uint<uintX>(NTL::rep(NTL::random_ZZ_p()), __l_query[j]);
	}
    }

    // cleanup
    delete [] streams;
    responses.kill();

    cudaFreeHost(l_query);
    cudaFree(d_query);
    cudaFreeHost(l_response);
    cudaFree(d_response);

    freeBarret<uintX>(barret);
    freeMatrix<uintX>(matrix);

    return 0;
}


template <typename T>
void initMatrix(const char * valfile, const char * rowfile,
	const char * colfile, NTL::ZZ & modulus,
	struct SparseMatrix<T> & matrix)
{
    std::ifstream valstream(valfile, std::ifstream::in);
    valstream >> modulus;
    valstream >> matrix.nvals;
    T * vals = (T *)malloc(matrix.nvals * sizeof(T));
    cudaMalloc((void**)&matrix.vals, matrix.nvals * sizeof(T));

    std::ifstream rowstream(rowfile, std::ifstream::in);
    uint * rows = (uint *)malloc(matrix.nvals * sizeof(uint));
    cudaMalloc((void**)&matrix.rows, matrix.nvals * sizeof(uint));

    std::ifstream colstream(colfile, std::ifstream::in);
    colstream >> matrix.ncols;
    uint * cols = (uint *)malloc(matrix.ncols * sizeof(uint));
    cudaMalloc((void**)&matrix.cols, matrix.ncols * sizeof(uint));

    NTL::ZZ_pPush p(modulus);
    for (int i = 0; i < matrix.nvals; i++)
    {
	NTL::ZZ_p tmp;
	valstream >> tmp;
	to_uint<T>(NTL::rep(tmp), vals[i]);
	rowstream >> rows[i];
    }
    valstream.close();
    rowstream.close();
    cudaMemcpy(matrix.vals, vals, matrix.nvals * sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(matrix.rows, rows, matrix.nvals * sizeof(uint), cudaMemcpyHostToDevice);

    for (int i = 0; i < matrix.ncols; i++)
    {
	colstream >> cols[i];
    }
    colstream.close();
    cudaMemcpy(matrix.cols, cols, matrix.ncols * sizeof(uint), cudaMemcpyHostToDevice);
    
    free(vals);
    free(rows);
    free(cols);
}

template <typename T>
void freeMatrix(struct SparseMatrix<T> & matrix)
{
    cudaFree(matrix.vals);
    cudaFree(matrix.cols);
    cudaFree(matrix.rows);
}

template <typename T>
void initBarret(const NTL::ZZ & modulus_zz, BarretParams<T> & barret)
{
    const NTL::ZZ mu_zz = NTL::power2_ZZ(2 * BITS_PER_LIMB * LIMBS_IN(T)) / modulus_zz;
    T modulus;
    to_uint<T>(modulus_zz, modulus);
    T mu;
    to_uint<T>(mu_zz, mu); // loses high order bit!

    T * subtrahends = (T *)malloc(barret.u * 2 * sizeof(T));
    for (int i = 0; i < barret.u; ++i)
    {
	NTL::ZZ subtrahend = ((NTL::to_ZZ(i) << (2 * BITS_PER_LIMB * LIMBS_IN(T))) / modulus_zz) * modulus_zz;
	NTL::BytesFromZZ((unsigned char *)&subtrahends[2*i], subtrahend, 2*sizeof(T));
    }

    cudaMalloc((void**)&barret.modulus, sizeof(T));
    cudaMemcpy(barret.modulus, &modulus, sizeof(T), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&barret.mu, sizeof(T));
    cudaMemcpy(barret.mu, &mu, sizeof(T), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&barret.subtrahends, 2 * barret.u * sizeof(T));
    cudaMemcpy(barret.subtrahends, subtrahends, 2 * barret.u * sizeof(T), cudaMemcpyHostToDevice);

    free(subtrahends);
}

template<typename T>
void freeBarret(struct BarretParams<T> & barret)
{
    cudaFree(barret.modulus);
    cudaFree(barret.mu);
    cudaFree(barret.subtrahends);
}