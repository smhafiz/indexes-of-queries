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
//#include <NTL/vec_ZZ_p.h>
#include "uintX.h"

NTL_CLIENT

template< typename T>
struct BarretParams
{
    const int u = 16;
    T * modulus;
    T * mu;
    uint2X<T> * subtrahends;
};

struct SparseMatrix
{
    SparseMatrix(uint nv, uint nc, uint nr) : nvals(nv), nrows(nr), ncols(nc) {}
    const uint nvals;
    uintX * vals;
    const uint ncols;
    uint * cols;
    const uint nrows;
    uint * rows;
};

template <typename T>
void initBarret(const NTL::ZZ & modulus_zz, struct BarretParams<T> & barret)
{
    const NTL::ZZ mu_zz = NTL::power2_ZZ(2 * BITS_PER_LIMB * sizeof(T)) / modulus_zz;
    const T modulus = to_uint<T>(modulus_zz);
    const T mu = to_uint<T>(mu_zz); // loses high order bit!

    uint2X<T> * subtrahends = (uint2X<T> *) malloc(barret.u * sizeof(uint2X<T>));
    for (int i = 0; i < barret.u; ++i)
    {
	NTL::ZZ subtrahend = ((NTL::to_ZZ(i) << (2 * BITS_PER_LIMB * sizeof(T))) / modulus_zz) * modulus_zz;
	subtrahends[i].lo = to_uint<T>(subtrahend);
	subtrahends[i].hi = to_uint<T>(subtrahend >> BITS_PER_LIMB * sizeof(T));
    }

    cudaMalloc((void**) & barret.modulus, sizeof(T));
    cudaMemcpy(barret.modulus, & modulus, sizeof(T), cudaMemcpyHostToDevice);

    cudaMalloc((void**) & barret.mu, sizeof(T));
    cudaMemcpy(barret.mu, & mu, sizeof(T), cudaMemcpyHostToDevice);

    cudaMalloc((void**) & barret.subtrahends, 2 * barret.u * sizeof(T));
    cudaMemcpy(barret.subtrahends, subtrahends, 2 * barret.u * sizeof(T), cudaMemcpyHostToDevice);

    free(subtrahends);
}

template<typename T>
void killBarret(struct BarretParams<T> & barret)
{
    cudaFree(barret.modulus);
    cudaFree(barret.mu);
    cudaFree(barret.subtrahends);
}

template <typename T>
__global__ void smm_kernel(T * result, const T * query, const struct SparseMatrix & matrix, const struct BarretParams<T> & barret)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    uint2X<uintX> a;
    a.lo = result[i];
    a.hi = { 0 };
    uint overflow;
    for (int j = matrix.cols[i]; j < matrix.cols[i+1]; ++j)
    {
	mad(a, overflow, matrix.vals[j], query[matrix.rows[j]]);
    }
    normalize(a, barret.subtrahends[overflow]);
    uintXp<T> q = get_q(a, *barret.mu);
    uintXp<T> r2 = get_r2(q, *barret.modulus);
    sub(a, r2);
    result[i] = a.lo;
}

template <typename T>
void smm(T * l_response, T * l_query, T * d_response, T * d_query, const cudaStream_t & stream, const struct SparseMatrix & matrix, const struct BarretParams<T> & barret)
{
    //NTL::ZZ_pPush p(barret.modulus);

    cudaMemcpyAsync(d_query, l_query, matrix.nrows * sizeof(uintX), cudaMemcpyHostToDevice, stream);
    
    smm_kernel<uintX> <<< 1, 1, 0, stream >>> (d_response, d_query, matrix, barret);

    cudaMemcpyAsync(l_response, d_response, matrix.ncols * sizeof(uintX), cudaMemcpyDeviceToHost, stream);
    
//    l_response.SetLength(matrix.ncols);
//    for (int i = 0; i < matrix.ncols; ++i)
//    {
//    	l_response[i] = to_ZZ_p<T>(d_response[i]);
//    }
}

int main(int argc, char ** argv)
{
    NTL::ZZ modulus = NTL::RandomPrime_ZZ(96);
    NTL::ZZ_p::init(modulus);

    int nthreads = 1;

    struct BarretParams<uintX> barret;
    initBarret<uintX>(modulus, barret);

    uintX * vals = NULL;
    uint * rows = NULL;
    uint * cols = NULL;

    struct SparseMatrix matrix(0, 0, 0);
    cudaMalloc((void**)&matrix.vals, matrix.nvals * sizeof(uintX));
    cudaMemcpy(matrix.vals, vals, matrix.nvals * sizeof(uintX), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&matrix.cols, matrix.ncols * sizeof(uint));
    cudaMemcpy(matrix.cols, cols, matrix.ncols * sizeof(uint), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&matrix.rows, matrix.nrows * sizeof(uint));
    cudaMemcpy(matrix.rows, rows, matrix.nrows * sizeof(uint), cudaMemcpyHostToDevice);

    uintX * l_query, * d_query;
    cudaMallocHost((void**)&l_query, nthreads * matrix.nrows * sizeof(uintX));
    cudaMalloc((void**)&d_query, nthreads * matrix.nrows * sizeof(uintX));
    uintX * l_response, * d_response;
    cudaMallocHost((void**)&l_response, nthreads * matrix.ncols * sizeof(uintX));
    cudaMalloc((void**)&d_response, nthreads * matrix.ncols * sizeof(uintX));
    
    cudaStream_t * streams = new cudaStream_t[nthreads];
    for (int i = 0; i < nthreads; ++i) cudaStreamCreate(&streams[i]);

    std::atomic<int> cnt = ATOMIC_VAR_INIT(0);
    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds onesec{1000000000};
    while(std::chrono::duration_cast<std::chrono::duration<int,std::nano>>(std::chrono::high_resolution_clock::now() - start) < onesec)
    {
	smm<uintX>(l_response, l_query, d_response, d_query, streams[0], matrix, barret);
	std::atomic_fetch_add(&cnt, 1);
    }




    delete [] streams;

    cudaFreeHost(l_query);
    cudaFree(d_query);
    cudaFreeHost(l_response);
    cudaFree(d_response);

    free(vals);
    free(rows);
    free(cols);

    cudaFree(matrix.vals);
    cudaFree(matrix.cols);
    cudaFree(matrix.rows);

    killBarret<uintX>(barret);

    return 0;
}
