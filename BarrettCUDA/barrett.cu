// This file is part of BarrettCUDA v0.1.
//
// BarrettCUDA is a fast(ish) implementation of finite field sparse
// matrix-vector multiplication (SpMV) for Nvidia GPU devices, written
// in CUDA C++. BarrettCUDA supports SpMV for matrices expressed in
// the 'compressed column storage' (CCS) sparse matrix representation
// over (i) the field of integers modulo an arbitrary multi-precision
// prime, or (ii) either of the binary fields GF(2^8) or GF(2^16).
//
// Copyright (C) 2016, Ryan Henry and Syed Mahbub Hafiz.
//
// BarrettCUDA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// BarrettCUDA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with BarrettCUDA. If not, see <http://www.gnu.org/licenses/>.

#include <atomic>
#include <chrono>
#include <ratio>
#include <fstream>

#include "barrett.h"
#include "uintX.h"

#define THREADS_PER_BLOCK(n) (n >= 512 ? 512 : n)
#define NUM_BLOCKS(n) ((n + THREADS_PER_BLOCK(n) - 1) / THREADS_PER_BLOCK(n))

#define DEBUG_IDX 107
#define GF216
#define DEBUG_IF_ENABLED false

NTL_CLIENT

// specialization for uintX
template <typename T> struct _SpMV_specializer<T,0>
{
    static __device__ void device_SpMV(T * response, const T * query,
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
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  a_lo:  \t"); _print_limbs<T>(res_lo.lo, LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  a_hi:  \t"); _print_limbs<T>(res_hi, LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  over:  \t%u", overflow);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  sub_lo:\t"); _print_limbs<T>(subtrahends[overflow], LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  sub_hi:\t"); _print_limbs<T>(subtrahends[2*overflow+1], LIMBS_PER_UINTX);}
    // do the Barrett reduction
    normalize(res_lo.lo, res_hi, subtrahends[overflow], (overflow ? -1: 0));
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  a_lo': \t"); _print_limbs<T>(res_lo.lo, LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  a_hi': \t"); _print_limbs<T>(res_hi, LIMBS_PER_UINTX);}
    uintXp<T> q = get_q(res_lo.lo, res_hi, *mu);
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  q_lo:  \t"); _print_limbs<T>(q.lo, LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  q_hi:  \t"); _print_limbs<uint>(q.hi, 1);}
    uintXp<T> r2 = get_r2(q, *modulus);
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  r2:    \t"); _print_limbs<uintXp<T>>(r2, LIMBS_PER_UINTX+1);}
    res_lo.hi = sub(res_lo.lo, res_hi, r2);
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  a_lo'':\t"); _print_limbs<T>(res_lo.lo, LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  a_hi'':\t"); _print_limbs<uint>(res_lo.hi, 1);}
    if (res_lo.hi) sub_modulus(res_lo, *modulus);
    if (res_lo.hi) sub_modulus(res_lo, *modulus);
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  a_lo'':\t"); _print_limbs<T>(res_lo.lo, LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){printf("\n  a_hi'':\t"); _print_limbs<uint>(res_lo.hi, 1);}
    // write final result to global memory
    response[i] = res_lo.lo;
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){ printf("\nCUDA Kernel  a''': \t"); _print_limbs<T>(response[i], LIMBS_PER_UINTX);}
    }
};

// specialization for GF28_Element
template <> struct _SpMV_specializer<GF28_Element,0>
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
if (!DEBUG_IF_ENABLED && i==DEBUG_IDX) printf("device_SpMV:GF28_Element:GPU response[%d]: %u\n",i,response[i]);
    }
};

// specialization for GF216_Element
template <> struct _SpMV_specializer<GF216_Element,0>
{
    static __device__ void device_SpMV(GF216_Element * response,
    const GF216_Element * query, const uint nvals,
    const GF216_Element * vals, const uint ncols, const uint * cols,
    const uint * rows)
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
        response[i] += to_ZZ_p(matrix.l_vals[j])
             * to_ZZ_p(query[matrix.l_rows[j]]);
    }
    if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\nNTL  a''': \t"; print_limbs<T>(response[i], LIMBS_PER_UINTX);}
    }
}


    template <typename T>
    void SpMV_ntl_barrett(NTL::vec_ZZ_p & response, const T * query,
    const SparseMatrix<T> & matrix, struct BarrettParams<T> & barrett)
    {
    NTL::vec_ZZ response_ZZ(INIT_SIZE, matrix.ncols);
    for (int i = 0; i < matrix.ncols; i++)
    {
        response_ZZ[i] = NTL::to_ZZ(0);

        for (int j = matrix.l_cols[i]; j < matrix.l_cols[i+1]; ++j)
        {
        response_ZZ[i] += to_ZZ(matrix.l_vals[j]) * to_ZZ(query[matrix.l_rows[j]]);
        }
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  a_lo:  \t"; print_limbs<T>(response_ZZ[i], LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  a_hi:  \t"; print_limbs<T>(response_ZZ[i] >> BITS_IN(LIMBS_PER_UINTX), LIMBS_PER_UINTX);}
        uint overflow = (uint)NTL::trunc_long(response_ZZ[i] >> 2*BITS_IN(LIMBS_PER_UINTX), BITS_IN(sizeof(uint)));
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  over:  \t" << overflow;}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  sub_lo:\t"; print_limbs<T>(barrett.l_subtrahends[overflow], LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  sub_hi:\t"; print_limbs<T>(barrett.l_subtrahends[2*overflow+1] >> BITS_IN(LIMBS_PER_UINTX), LIMBS_PER_UINTX);}
        response_ZZ[i] -= barrett.l_subtrahends[overflow];
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  a_lo': \t"; print_limbs<T>(response_ZZ[i], LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  a_hi': \t"; print_limbs<T>(response_ZZ[i] >> BITS_IN(LIMBS_PER_UINTX), LIMBS_PER_UINTX);}
        NTL::ZZ q1 = response_ZZ[i] >> BITS_IN(LIMBS_PER_UINTX-1);
        NTL::ZZ q2 = q1 * barrett.l_mu;
        NTL::ZZ q3 = q2 >> BITS_IN(LIMBS_PER_UINTX+1);
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  q_lo:  \t"; print_limbs<T>(q3, LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  q_hi:  \t"; print_limbs<uint>(q3 >> BITS_IN(LIMBS_PER_UINTX), 1);}

        NTL::ZZ r1 = response_ZZ[i] % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
        NTL::ZZ r2 = q3 * barrett.l_modulus % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  r2:    \t"; print_limbs<uintXp<T>>(r2, LIMBS_PER_UINTX+1);}
        NTL::ZZ r = (r1 - r2) % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  a_lo'':\t"; print_limbs<T>(r, LIMBS_PER_UINTX);}
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\n  a_hi'':\t"; print_limbs<uint>(r >> BITS_IN(LIMBS_PER_UINTX), 1);}
        response[i] = NTL::to_ZZ_p(r);
if (DEBUG_IF_ENABLED && i==DEBUG_IDX){std::cout << "\nNTL_BARRET  a''':     \t"; print_limbs<T>(response[i], LIMBS_PER_UINTX);}
    }
    }
#endif // DEBUG

template <typename T>
void SpMV(T * l_response, const T * l_query,
    T * d_response, T * d_query, const cudaStream_t & stream,
    const SparseMatrix<T> & matrix)
{
    gpuErrchk(cudaMemcpyAsync(d_query, l_query, matrix.nrows * sizeof(T),
    cudaMemcpyHostToDevice, stream));

    const dim3 Dg(NUM_BLOCKS(matrix.ncols), 1, 1);
    const dim3 Db(THREADS_PER_BLOCK(matrix.ncols), 1, 1);
    const size_t Ns = 0;

    SpMV_kernel<T> <<< Dg, Db, Ns, stream >>> (d_response, d_query,
    matrix.nvals, matrix.d_vals, matrix.ncols, matrix.d_cols, matrix.d_rows);
    gpuErrchk(cudaPeekAtLastError());

    gpuErrchk(cudaMemcpyAsync(l_response, d_response, matrix.ncols * sizeof(T),
    cudaMemcpyDeviceToHost, stream));
}

int main(int argc, char ** argv)
{
    int nstreams = 4;

    if (argc < 4)
    {
    std::cout << "Usage: " << argv[0] << " VALUES ROWS COLS\n\n";
    return 1;
    }

    time_t t0 = time(0);
    //t0 = 1482014996;
    NTL::SetSeed(to_ZZ(t0));
    std::cout << "seed: " << t0 << "\n";



    #ifdef uintX
    struct SparseMatrix<uintX> matrix = { 0 };
    NTL::ZZ modulus;

    uint max_overflow;
    initMatrix(argv[1], argv[2], argv[3], modulus, matrix, max_overflow);
    NTL::ZZ_p::init(modulus);

    struct BarrettParams<uintX> barrett;
    initBarrett<uintX>(modulus, barrett, max_overflow);

    uintX * l_query, * d_query;
    gpuErrchk(cudaMallocHost((void**)&l_query,
        nstreams * matrix.nrows * sizeof(uintX)));
    gpuErrchk(cudaMalloc((void**)&d_query,
        nstreams * matrix.nrows * sizeof(uintX)));

    uintX * l_response, * d_response;
    gpuErrchk(cudaMallocHost((void**)&l_response,
        nstreams * matrix.ncols * sizeof(uintX)));
    gpuErrchk(cudaMalloc((void**)&d_response,
        nstreams * matrix.ncols * sizeof(uintX)));
    #endif

    #ifdef GF28
    struct SparseMatrix<GF28_Element> matrix = { 0 };
    NTL::ZZ modulus;

    uint max_overflow;
    initMatrix(argv[1], argv[2], argv[3], modulus, matrix, max_overflow);
    NTL::ZZ_p::init(modulus);

    initGF28();

    GF28_Element * l_query, * d_query;
    gpuErrchk(cudaMallocHost((void**)&l_query,
        nstreams * matrix.nrows * sizeof(GF28_Element)));
    gpuErrchk(cudaMalloc((void**)&d_query,
        nstreams * matrix.nrows * sizeof(GF28_Element)));

    GF28_Element * l_response, * d_response;
    gpuErrchk(cudaMallocHost((void**)&l_response,
        nstreams * matrix.ncols * sizeof(GF28_Element)));
    gpuErrchk(cudaMalloc((void**)&d_response,
        nstreams * matrix.ncols * sizeof(GF28_Element)));
    #endif

    #ifdef GF216
    struct SparseMatrix<GF216_Element> matrix = { 0 };
    NTL::ZZ modulus;

    uint max_overflow;
    initMatrix(argv[1], argv[2], argv[3], modulus, matrix, max_overflow);
    NTL::ZZ_p::init(modulus);

    initGF216();

    GF216_Element * l_query, * d_query;
    gpuErrchk(cudaMallocHost((void**)&l_query,
        nstreams * matrix.nrows * sizeof(GF216_Element)));
    gpuErrchk(cudaMalloc((void**)&d_query,
        nstreams * matrix.nrows * sizeof(GF216_Element)));

    GF216_Element * l_response, * d_response;
    gpuErrchk(cudaMallocHost((void**)&l_response,
        nstreams * matrix.ncols * sizeof(GF216_Element)));
    gpuErrchk(cudaMalloc((void**)&d_response,
        nstreams * matrix.ncols * sizeof(GF216_Element)));
    #endif

    cudaStream_t * streams = new cudaStream_t[nstreams];
    for (int i = 0; i < nstreams; ++i) cudaStreamCreate(&streams[i]);

    NTL::vec_vec_ZZ_p responses(INIT_SIZE, nstreams,
    NTL::vec_ZZ_p(INIT_SIZE, matrix.ncols));

    for (int i = 0; i < nstreams * matrix.nrows; i++)
    {
    #ifdef uintX
        to_uint<uintX>(NTL::random_ZZ_p(), l_query[i]);
    #endif
    #ifdef GF28
        l_query[i] = NTL::RandomBits_long(8);
    #endif
    #ifdef GF216
        l_query[i] = NTL::RandomBits_long(16);
    #endif
    }

    std::atomic<int> cnt = ATOMIC_VAR_INIT(0);
    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds onesec{1000000000};

    //while (std::chrono::duration_cast<std::chrono::duration<int,std::nano>>(std::chrono::high_resolution_clock::now() - start) < onesec)
    {
    #pragma omp parallel
    for (int i = 0; i < nstreams; i++)
    {
       
        #ifdef uintX
        uintX * __l_response = l_response + i * matrix.ncols;
        uintX * __d_response = d_response + i * matrix.ncols;
        uintX * __l_query = l_query + i * matrix.nrows;
        uintX * __d_query = d_query + i * matrix.nrows;
        SpMV<uintX>(__l_response, __l_query, __d_response,
        __d_query, streams[i], matrix);
        #ifdef DEBUG
            SpMV_ntl(responses[i], __l_query, matrix);
            SpMV_ntl_barrett(responses[i], __l_query, matrix, barrett);
        #endif
        #endif
       
        #ifdef GF28
        GF28_Element * __l_response = l_response + i * matrix.ncols;
        GF28_Element * __d_response = d_response + i * matrix.ncols;
        GF28_Element * __l_query = l_query + i * matrix.nrows;
        GF28_Element * __d_query = d_query + i * matrix.nrows;
        SpMV<GF28_Element>(__l_response, __l_query, __d_response,
        __d_query, streams[i], matrix);
        #endif


        #ifdef GF216
        GF216_Element * __l_response = l_response + i * matrix.ncols;
        GF216_Element * __d_response = d_response + i * matrix.ncols;
        GF216_Element * __l_query = l_query + i * matrix.nrows;
        GF216_Element * __d_query = d_query + i * matrix.nrows;
        SpMV<GF216_Element>(__l_response, __l_query, __d_response,
        __d_query, streams[i], matrix);
        #endif


        std::atomic_fetch_add(&cnt, 1);
    }
    }
    std::cout << "Completed " << cnt << " SpMVs in 1 second\n";

    // cleanup
    for (int i = 0; i < nstreams; ++i) gpuErrchk(cudaStreamDestroy(streams[i]));
    delete [] streams;
    responses.kill();

    gpuErrchk(cudaFreeHost(l_query));
    gpuErrchk(cudaFree(d_query));
    gpuErrchk(cudaFreeHost(l_response));
    gpuErrchk(cudaFree(d_response));

    #ifdef uintX
    freeBarrett<uintX>(barrett);
    freeMatrix<uintX>(matrix);
    #endif
    #ifdef GF28
    freeGF28();
    freeMatrix<GF28_Element>(matrix);
    #endif
    #ifdef GF216
    freeGF216();
    freeMatrix<GF216_Element>(matrix);
    #endif
    return 0;
}


template <typename T>
void initMatrix(const char * valfile, const char * rowfile,
    const char * colfile, NTL::ZZ & modulus,
    struct SparseMatrix<T> & matrix, uint & max_overflow)
{
    std::ifstream valstream(valfile, std::ifstream::in);
    if (!valstream) { cerr << "Error: CCS representation: can't open VALS file\n"; exit(-1); }
    std::ifstream rowstream(rowfile, std::ifstream::in);
    if (!rowstream) { cerr << "Error: CCS representation: can't open ROWS file\n"; exit(-1); }
    std::ifstream colstream(colfile, std::ifstream::in);
    if (!colstream) { cerr << "Error: CCS representation: can't open COLS file\n"; exit(-1); }

    NTL::ZZ tmp_zz;
    valstream >> tmp_zz;
    if(!valstream.good())
    {   
    cerr << "Error: CCS representation: VALS file: File has no more value.\n"; exit(-1); //good
    }
    if(tmp_zz <= NTL::to_ZZ(0))
    {
    cerr << "Error: CCS representation: VALS file: Modulus should be positive.\n"; exit(-1);//good
    }
    if(NTL::NumBits(tmp_zz) >= (sizeof(T)*8))
    {
    modulus = NTL::trunc_ZZ(tmp_zz,sizeof(T)*8);
    } else {
    cerr << "Error: CCS representation: VALS file: Modulus length is small\n"; exit(-1);//good
    }


    rowstream >> matrix.nrows;
    if(!rowstream.good())
    {
    cerr << "Error: CCS representation: ROWS file: File has no more value.\n"; exit(-1);
    }
    if(matrix.nrows <= NTL::to_ZZ(0))
    {
    cerr << "Error: CCS representation: ROWS file: Number of rows (p) should be positive.\n"; exit(-1);
    }


    rowstream >> matrix.nvals;
    if(!rowstream.good())
    {
    cerr << "Error: CCS representation: ROWS file: File has no more value.\n"; exit(-1);
    }
    if(matrix.nvals < NTL::to_ZZ(0))
    {
    cerr << "Error: CCS representation: ROWS file: Number of non-zero values should be non-negitive.\n"; exit(-1);
    }


   
    matrix.l_rows = (uint *)malloc(matrix.nvals * sizeof(uint));
    gpuErrchk(cudaMalloc((void**)&matrix.d_rows, matrix.nvals * sizeof(uint)));
    matrix.l_vals = (T *)malloc(matrix.nvals * sizeof(T));
    gpuErrchk(cudaMalloc((void**)&matrix.d_vals, matrix.nvals * sizeof(T)));

    colstream >> matrix.ncols;
    if(!colstream.good())
    {
    cerr << "Error: CCS representation: COLS file: File has no value.\n"; exit(-1);//good
    }
    if(matrix.ncols <= NTL::to_ZZ(0))
    {
    cerr << "Error: CCS representation: COLS file: Number of columns (r) should be positive.\n"; exit(-1);//good
    }



    matrix.l_cols = (uint *)malloc((matrix.ncols+1) * sizeof(uint));
    gpuErrchk(cudaMalloc((void**)&matrix.d_cols, (matrix.ncols+1) * sizeof(uint)));

    NTL::ZZ_pPush p(modulus);
    for (int i = 0; i < matrix.nvals; i++)
    {
    NTL::ZZ_p tmp;
    valstream >> tmp;
    if(!valstream.good())
    {
        cerr << "Error: CCS representation:  VALS file: Number of non-zero values is less than the provided length.\n"; exit(-1);//good
    }
   
    #ifdef uintX
        to_uint<T>(NTL::rep(tmp), matrix.l_vals[i]);
    #endif
    #ifdef GF28
        matrix.l_vals[i] = NTL::trunc_long(tmp_zz, 8);
    #endif
    #ifdef GF216
        matrix.l_vals[i] = NTL::trunc_long(tmp_zz, 16);
    #endif
        rowstream >> matrix.l_rows[i];
    if(!rowstream.good())
    {
        cerr << "Error: CCS representation:  ROWS file: each value should have a row number.\n"; exit(-1);//good
    }
    if(matrix.l_rows[i] < NTL::to_ZZ(0) || matrix.l_rows[i] >= matrix.nrows)
    {
        cerr << "Error: CCS representation: ROWS file: Invalid row number.\n"; exit(-1);//good
    }
   
    }
    valstream.close();
    rowstream.close();
    gpuErrchk(cudaMemcpy(matrix.d_vals, matrix.l_vals, matrix.nvals * sizeof(T),
    cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(matrix.d_rows, matrix.l_rows,
    matrix.nvals * sizeof(uint), cudaMemcpyHostToDevice));
   
    long sum_vals = 0;
    for (int i = 0; i < matrix.ncols+1; i++)
    {
    colstream >> matrix.l_cols[i];
    if(!colstream.good())
    {
        cerr << "Error: CCS representation:  COLS file: Number of column pointer is less than the provided length.\n"; exit(-1);//good
    }
    if (i>0) {
        if( matrix.l_cols[i] < matrix.l_cols[i-1] || matrix.l_cols[i] > matrix.nvals || matrix.l_cols[i] < NTL::to_ZZ(0))
        {
            cerr << "Error: CCS representation: COLS file: Invalid column pointer given.\n"; exit(-1);//good
        }
        sum_vals += matrix.l_cols[i] - matrix.l_cols[i-1];
    }
    }
    if(sum_vals < matrix.nvals)
    {
    cerr << "Error: CCS representation:  COLS file: More non-zero values left to put in columns.\n"; exit(-1);//good
    }
    colstream.close();
    #ifdef uintX
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
    #endif
    gpuErrchk(cudaMemcpy(matrix.d_cols, matrix.l_cols,
    (matrix.ncols+1) * sizeof(uint), cudaMemcpyHostToDevice));
}

template <typename T>
void freeMatrix(struct SparseMatrix<T> & matrix)
{
    free(matrix.l_vals);
    free(matrix.l_rows);
    free(matrix.l_cols);
    gpuErrchk(cudaFree(matrix.d_vals));
    gpuErrchk(cudaFree(matrix.d_cols));
    gpuErrchk(cudaFree(matrix.d_rows));
}

template <typename T>
void initBarrett(const NTL::ZZ & modulus_zz, BarrettParams<T> & barrett,
    const uint max_overflow)
{
    barrett.l_modulus = modulus_zz;
    barrett.l_mu = NTL::power2_ZZ(2*BITS_PER_LIMB*LIMBS_PER_UINTX) / modulus_zz;
    T modulus;
    to_uint<T>(modulus_zz, modulus);
    uintXp<T> mu;
    to_uint<uintXp<T>>(barrett.l_mu, mu);

    barrett.l_subtrahends.SetLength(max_overflow+1);
    T * subtrahends = (T *)malloc((max_overflow+1) * sizeof(T));
    for (int i = 0; i <= max_overflow; ++i)
    {
    barrett.l_subtrahends[i] = ((NTL::to_ZZ(i)
        << (2*BITS_PER_LIMB*LIMBS_PER_UINTX)) / modulus_zz) * modulus_zz;
    NTL::BytesFromZZ((unsigned char *)&subtrahends[i],
        barrett.l_subtrahends[i], LIMBS_PER_UINTX * sizeof(uint));
    }

    gpuErrchk(cudaMalloc((void**)&barrett.d_modulus, sizeof(T)));
    gpuErrchk(cudaMemcpy(barrett.d_modulus, &modulus, sizeof(T),
    cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(d_modulus, &barrett.d_modulus, sizeof(T *)));

    gpuErrchk(cudaMalloc((void**)&barrett.d_mu, sizeof(uintXp<T>)));
    gpuErrchk(cudaMemcpy(barrett.d_mu, &mu, sizeof(uintXp<T>),
    cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(d_mu, &barrett.d_mu, sizeof(uintXp<T> *)));

    gpuErrchk(cudaMalloc((void**)&barrett.d_subtrahends,
    (max_overflow+1) * sizeof(T)));
    gpuErrchk(cudaMemcpy(barrett.d_subtrahends, subtrahends,
    (max_overflow+1) * sizeof(T), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(d_subtrahends, &barrett.d_subtrahends,
    sizeof(T *)));

    free(subtrahends);
}

template<typename T>
void freeBarrett(struct BarrettParams<T> & barrett)
{
    barrett.l_subtrahends.kill();
    barrett.l_modulus.kill();
    barrett.l_mu.kill();
    gpuErrchk(cudaFree(barrett.d_modulus));
    gpuErrchk(cudaFree(barrett.d_mu));
    gpuErrchk(cudaFree(barrett.d_subtrahends));
}

void initGF28()
{
    GF28_Element * __temp;
    gpuErrchk(cudaMalloc((void**)&__temp, sizeof(GF28_mult_table)));
    gpuErrchk(cudaMemcpy(__temp, &GF28_mult_table, sizeof(GF28_mult_table),
    cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(d_GF28_mult_table, &__temp,
    sizeof(GF28_Element *)));
}

void freeGF28()
{
    void * mult_table;
    gpuErrchk(cudaGetSymbolAddress((void**)&mult_table, d_GF28_mult_table));
    gpuErrchk(cudaFree(mult_table));
}

void initGF216()
{
    void * exp_table, * log_table;
    gpuErrchk(cudaMalloc((void**)&exp_table, sizeof(GF216_exp_table)));
    gpuErrchk(cudaMalloc((void**)&log_table, sizeof(GF216_log_table)));
    gpuErrchk(cudaMemcpy(exp_table, &GF216_exp_table, sizeof(GF216_exp_table),
    cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(log_table, &GF216_log_table, sizeof(GF216_log_table),
    cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(d_GF216_exp_table, &exp_table,
    sizeof(GF216_Element *)));
    gpuErrchk(cudaMemcpyToSymbol(d_GF216_log_table, &log_table,
    sizeof(GF216_Element *)));
}

void freeGF216()
{
    void * exp_table, * log_table;
    gpuErrchk(cudaGetSymbolAddress((void**)&exp_table, d_GF216_exp_table));
    gpuErrchk(cudaGetSymbolAddress((void**)&log_table, d_GF216_exp_table));
    gpuErrchk(cudaFree(exp_table));
    gpuErrchk(cudaFree(log_table));
}
