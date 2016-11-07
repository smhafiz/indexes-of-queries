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
#include <chrono>
#include <ratio>
#include <fstream>

#include "barret.h"
#include "uintX.h"

#define NUM_BLOCKS(n) 		(n >= 256 ? 256 : n)
#define THREADS_PER_BLOCK(n)	((n + NUM_BLOCKS(n) - 1) / NUM_BLOCKS(n))

NTL_CLIENT

// partial specialization for uintX
template <typename T> struct _SpMV_specializer<T, uintXp<T>>
{
    static __device__ void device_SpMV(uintXp<T> * response, const T * query,
    	const uint nvals, const T * vals, const uint ncols, const uint * cols,
    	const uint * rows)
    {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= ncols) return;

	uintXp<T> res_lo = { 0 };
	T res_hi = { 0 };
	uint overflow = { 0 };

	// do the SpMV
	for (int j = cols[i]; j < cols[i+1]; ++j)
	{
	    mad(res_lo.lo, res_hi, overflow, vals[j], query[rows[j]]);
	}

	T * subtrahends = (T *)d_subtrahends;
	uintXp<T> * mu = (uintXp<T> *)d_mu;
	T * modulus = (T *)d_modulus;

	// do the Barret reduction
	normalize(res_lo.lo, res_hi, subtrahends[overflow], (overflow ? -1: 0));
	uintXp<T> q = get_q(res_lo.lo, res_hi, *mu);
	uintXp<T> r2 = get_r2(q, *modulus);
	res_lo.hi = sub(res_lo.lo, res_hi, r2);

	// write final result to global memory
	response[i] = res_lo;
    }
};

// full specialization for GF28_Element
template <> struct _SpMV_specializer<GF28_Element, GF28_Element>
{
    static __device__ void device_SpMV(GF28_Element * response,
    	const GF28_Element * query, const uint nvals, const GF28_Element * vals,
    	const uint ncols, const uint * cols, const uint * rows)
    {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= ncols) return;

	GF28_Element res = 0;
	for (int j = cols[i]; j < cols[i+1]; ++j)
	{
	    res ^= d_GF28_mult_table[vals[j]][query[rows[j]]];
	}
	response[i] = res;
    }
};

// full specialization for GF28_Element
template <> struct _SpMV_specializer<GF216_Element, GF216_Element>
{
    static __device__ void device_SpMV(GF216_Element * response,
    	const GF216_Element * query, const uint nvals, const GF216_Element * vals,
    	const uint ncols, const uint * cols, const uint * rows)
    {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= ncols) return;

	GF216_Element res = 0;
	for (int j = cols[i]; j < cols[i+1]; ++j)
	{
	    GF216_Element log_x = d_GF216_log_table[vals[j]];
	    GF216_Element log_y = d_GF216_log_table[query[rows[j]]];
	    res ^= d_GF216_exp_table[log_x+log_y];
	}
	response[i] = res;
    }
};


#ifdef DEBUG
    template <typename T>
    void SpMV_ntl(NTL::vec_ZZ_p & response, const T * query,
	const SparseMatrix<T> & matrix)
    {
	for (int i = 0; i < matrix.ncols; i++)
	{
	    response[i] = NTL::to_ZZ_p(0);
	    for (int j = matrix.l_cols[i]; j < matrix.l_cols[i+1]; ++j)
	    {
		response[i] += to_ZZ_p(matrix.l_vals[j]) * to_ZZ_p(query[matrix.l_rows[j]]);
	    }
	}
    }

    template <typename T>
    void SpMV_ntl_barret(NTL::vec_ZZ_p & response, const T * query,
	const SparseMatrix<T> & matrix, struct BarretParams<T> & barret)
    {
	NTL::vec_ZZ response_ZZ(INIT_SIZE, matrix.ncols);
	for (int i = 0; i < matrix.ncols; i++)
	{
	    response_ZZ[i] = NTL::to_ZZ(0);

	    for (int j = matrix.l_cols[i]; j < matrix.l_cols[i+1]; ++j)
	    {
		response_ZZ[i] += to_ZZ(matrix.l_vals[j]) * to_ZZ(query[matrix.l_rows[j]]);
	    }
	    uint overflow = (uint)NTL::trunc_long(response_ZZ[i] >> 2*BITS_IN(LIMBS_PER_UINTX), BITS_IN(sizeof(uint)));
	    response_ZZ[i] -= barret.l_subtrahends[overflow];
	    NTL::ZZ q1 = response_ZZ[i] >> BITS_IN(LIMBS_PER_UINTX-1);
	    NTL::ZZ q2 = q1 * barret.l_mu;
	    NTL::ZZ q3 = q2 >> BITS_IN(LIMBS_PER_UINTX+1);
	    NTL::ZZ r1 = response_ZZ[i] % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
	    NTL::ZZ r2 = q3 * barret.l_modulus % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
	    NTL::ZZ r = (r1 - r2) % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
	    response[i] = NTL::to_ZZ_p(r);
	}
    }
#endif // DEBUG

template <typename T>
void SpMV(NTL::vec_ZZ_p & response, uintXp<T> * l_response, const T * l_query,
	uintXp<T> * d_response, T * d_query, const cudaStream_t & stream,
	const SparseMatrix<T> & matrix)
{
    cudaMemcpyAsync(d_query, l_query, matrix.nrows * sizeof(T),
	cudaMemcpyHostToDevice, stream);

    const dim3 Dg(NUM_BLOCKS(matrix.ncols), 1, 1);
    const dim3 Db(THREADS_PER_BLOCK(matrix.ncols), 1, 1);
    const size_t Ns = 0;

    SpMV_kernel<T, uintXp<T>> <<< Dg, Db, Ns, stream >>> (d_response, d_query,
	matrix.nvals, matrix.d_vals, matrix.ncols, matrix.d_cols, matrix.d_rows);

    cudaMemcpyAsync(l_response, d_response, matrix.ncols * sizeof(uintXp<T>),
	cudaMemcpyDeviceToHost, stream);

    cudaStreamSynchronize(stream);
    response.SetLength(matrix.ncols);
    for (int i = 0; i < matrix.ncols; ++i)
    {
	response[i] = to_ZZ_p<T>(l_response[i].lo)
	    + NTL::to_ZZ_p(NTL::to_ZZ(l_response[i].hi) << BITS_IN(LIMBS_PER_UINTX));
    }
}

template <typename GF2E_Element>
void SpMV(GF2E_Element * l_response, const GF2E_Element * l_query,
	GF2E_Element * d_response, GF2E_Element * d_query, const cudaStream_t & stream,
	const SparseMatrix<GF2E_Element> & matrix)
{
    cudaMemcpyAsync(d_query, l_query, matrix.nrows * sizeof(GF2E_Element),
	cudaMemcpyHostToDevice, stream);

    const dim3 Dg(NUM_BLOCKS(matrix.ncols), 1, 1);
    const dim3 Db(THREADS_PER_BLOCK(matrix.ncols), 1, 1);
    const size_t Ns = 0;

    SpMV_kernel <GF2E_Element, GF2E_Element> <<< Dg, Db, Ns, stream >>> (d_response, d_query,
	matrix.nvals, matrix.d_vals, matrix.ncols, matrix.d_cols, matrix.d_rows);

    cudaMemcpyAsync(l_response, d_response, matrix.ncols * sizeof(GF2E_Element),
	cudaMemcpyDeviceToHost, stream);
}

int main(int argc, char ** argv)
{
    int nstreams = 4;

    if (argc < 3)
    {
	std::cout << "Usage: " << argv[0] << " VALUES ROWS COLS\n\n";
	return 1;
    }

    time_t t0 = time(0);
    NTL::SetSeed(to_ZZ(t0));
    std::cout << "seed: " << t0 << "\n";

    struct SparseMatrix<uintX> matrix = { 0 };
    NTL::ZZ modulus;

    uint max_overflow;
    initMatrix(argv[1], argv[2], argv[3], modulus, matrix, max_overflow);
    NTL::ZZ_p::init(modulus);

    struct BarretParams<uintX> barret;
    initBarret<uintX>(modulus, barret, max_overflow);

    uintX * l_query, * d_query;
    cudaMallocHost((void**)&l_query, nstreams * matrix.nrows * sizeof(uintX));
    cudaMalloc((void**)&d_query, nstreams * matrix.nrows * sizeof(uintX));

    uintXp<uintX> * l_response, * d_response;
    cudaMallocHost((void**)&l_response, nstreams * matrix.ncols * sizeof(uintXp<uintX>));
    cudaMalloc((void**)&d_response, nstreams * matrix.ncols * sizeof(uintXp<uintX>));

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
	#pragma omp parallel
	for (int i = 0; i < nstreams; i++)
	{
	    uintXp<uintX> * __l_response = l_response + i * matrix.ncols;
	    uintXp<uintX> * __d_response = d_response + i * matrix.ncols;
	    uintX * __l_query = l_query + i * matrix.nrows;
	    uintX * __d_query = d_query + i * matrix.nrows;		

	    SpMV<uintX>(responses[i], __l_response, __l_query, __d_response,
		__d_query, streams[i], matrix);
//	    SpMV_ntl(responses[i], __l_query, matrix);
//	    SpMV_ntl_barret(responses[i], __l_query, matrix, barret);

	    std::atomic_fetch_add(&cnt, 1);
//	    for (int j = 0; j < matrix.nrows; j++)
//	    {
//		to_uint<uintX>(NTL::rep(NTL::random_ZZ_p()), __l_query[j]);
//	    }
	}
    }

    std::cout << "Count: " << cnt << "\n";

    // cleanup
    for (int i = 0; i < nstreams; ++i) cudaStreamDestroy(streams[i]);
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
	struct SparseMatrix<T> & matrix, uint & max_overflow)
{
    std::ifstream valstream(valfile, std::ifstream::in);
    if (!valstream) { cerr << "Error: opening VALS files\n"; exit(-1); }
    std::ifstream rowstream(rowfile, std::ifstream::in);
    if (!rowstream) { cerr << "Error: opening VALS files\n"; exit(-1); }
    std::ifstream colstream(colfile, std::ifstream::in);
    if (!colstream) { cerr << "Error: opening VALS files\n"; exit(-1); }

    valstream >> modulus;

    rowstream >> matrix.nrows;
    valstream >> matrix.nvals;
    matrix.l_rows = (uint *)malloc(matrix.nvals * sizeof(uint));
    cudaMalloc((void**)&matrix.d_rows, matrix.nvals * sizeof(uint));
    matrix.l_vals = (T *)malloc(matrix.nvals * sizeof(T));
    cudaMalloc((void**)&matrix.d_vals, matrix.nvals * sizeof(T));

    colstream >> matrix.ncols;
    matrix.l_cols = (uint *)malloc((matrix.ncols+1) * sizeof(uint));
    cudaMalloc((void**)&matrix.d_cols, (matrix.ncols+1) * sizeof(uint));

//std::cout << "modulus:\t" << modulus << "\n";
//std::cout << "matrix.nrows:\t" << matrix.nrows << "\n";
//std::cout << "matrix.ncols:\t" << matrix.ncols << "\n";
//std::cout << "matrix.nvals:\t" << matrix.nvals << "\n";

    NTL::ZZ_pPush p(modulus);
    for (int i = 0; i < matrix.nvals; i++)
    {
	NTL::ZZ_p tmp;
	valstream >> tmp;
	to_uint<T>(NTL::rep(tmp), matrix.l_vals[i]);
	rowstream >> matrix.l_rows[i];
    }
    valstream.close();
    rowstream.close();
    cudaMemcpy(matrix.d_vals, matrix.l_vals, matrix.nvals * sizeof(T),
	cudaMemcpyHostToDevice);
    cudaMemcpy(matrix.d_rows, matrix.l_rows, matrix.nvals * sizeof(uint),
	cudaMemcpyHostToDevice);

    for (int i = 0; i < matrix.ncols+1; i++)
    {
	colstream >> matrix.l_cols[i];
    }
    colstream.close();

    NTL::ZZ max_col = NTL::to_ZZ(0);
    for (int i = 0; i < matrix.ncols; i++)
    {
    	NTL::ZZ this_col = NTL::to_ZZ(0);
	for (int j = matrix.l_cols[i]; j < matrix.l_cols[i+1]; j++)
	{
	    this_col += to_ZZ(matrix.l_vals[j]);
	}
	max_col = (max_col > this_col) ? max_col : this_col;
    }
    max_col *= (modulus-1);
    max_col >>= (2*BITS_PER_LIMB*LIMBS_PER_UINTX);
    max_overflow = (uint)trunc_long(max_col, 32);

    cudaMemcpy(matrix.d_cols, matrix.l_cols, (matrix.ncols+1) * sizeof(uint),
	cudaMemcpyHostToDevice);
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

template <typename T>
void initBarret(const NTL::ZZ & modulus_zz, BarretParams<T> & barret,
	const uint max_overflow)
{
    barret.l_modulus = modulus_zz;
    barret.l_mu = NTL::power2_ZZ(2*BITS_PER_LIMB*LIMBS_PER_UINTX) / modulus_zz;
    T modulus;
    to_uint<T>(modulus_zz, modulus);
    uintXp<T> mu;
    to_uint<uintXp<T>>(barret.l_mu, mu);

    barret.l_subtrahends.SetLength(max_overflow+1);
    T * subtrahends = (T *)malloc((max_overflow+1) * sizeof(T));
    for (int i = 0; i <= max_overflow; ++i)
    {
	barret.l_subtrahends[i] = ((NTL::to_ZZ(i) << (2*BITS_PER_LIMB*LIMBS_PER_UINTX))
	    / modulus_zz) * modulus_zz;
	NTL::BytesFromZZ((unsigned char *)&subtrahends[i], barret.l_subtrahends[i],
	    LIMBS_PER_UINTX * sizeof(uint));
    }

    cudaMalloc((void**)&d_modulus, sizeof(T));
    cudaGetSymbolAddress((void **)&barret.d_modulus, d_modulus);
    cudaMemcpy(barret.d_modulus, &modulus, sizeof(T), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_mu, sizeof(uintXp<T>));
    cudaGetSymbolAddress((void **)&barret.d_mu, d_mu);
    cudaMemcpy(barret.d_mu, &mu, sizeof(uintXp<T>), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_subtrahends, (max_overflow+1) * sizeof(T));
    cudaGetSymbolAddress((void **)&barret.d_subtrahends, d_subtrahends);
    cudaMemcpy(barret.d_subtrahends, subtrahends, (max_overflow+1) * sizeof(T),
	cudaMemcpyHostToDevice);

    free(subtrahends);
}

template<typename T>
void freeBarret(struct BarretParams<T> & barret)
{
    barret.l_subtrahends.kill();
    barret.l_modulus.kill();
    barret.l_mu.kill();
    cudaFree(barret.d_modulus);
    cudaFree(barret.d_mu);
    cudaFree(barret.d_subtrahends);
}
