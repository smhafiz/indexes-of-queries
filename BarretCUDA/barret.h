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

#ifndef __BARRET_H_
#define __BARRET_H_

#include <NTL/vec_vec_ZZ_p.h>

#include "uint.h"
#include "gf2e.h"

// for uintX
__constant__ void * d_modulus;
__constant__ void * d_mu;
__constant__ void * d_subtrahends;
// for GF28
__constant__ GF28_Element * d_GF28_mult_table[256];
// for GF216
__constant__ GF216_Element * d_GF216_exp_table;
__constant__ GF216_Element * d_GF216_log_table;

template<typename T>
struct BarretParams
{
    T * d_modulus;
    uintXp<T> * d_mu;
    T * d_subtrahends;
    NTL::ZZ l_modulus;
    NTL::ZZ l_mu;
    NTL::vec_ZZ l_subtrahends;
};

template<typename T>
struct SparseMatrix
{
    uint nvals;
    uint nrows;
    uint ncols;
    T * d_vals;
    uint * d_cols;
    uint * d_rows;
    T * l_vals;
    uint * l_cols;
    uint * l_rows;
};

// template wizardry to support partial specialization for uintX
template <typename T, typename U>
struct _SpMV_specializer
{
    static __device__ void device_SpMV(U * response, const T * query, const uint nvals,
	const T * vals, const uint ncols, const uint * cols, const uint * rows);
};

template <typename T, typename U>
__global__ void SpMV_kernel(U * response, const T * query, const uint nvals,
	const T * vals, const uint ncols, const uint * cols, const uint * rows)
{
    _SpMV_specializer<T,U>::device_SpMV(response, query, nvals, vals, ncols, cols, rows);
}

template <typename T>
void initMatrix(const char * valfile, const char * rowfile,
	const char * colfile, NTL::ZZ & modulus,
	struct SparseMatrix<T> & matrix, uint & max_overflow);

template <typename T>
void freeMatrix(struct SparseMatrix<T> & matrix);

template <typename T>
void initBarret(const NTL::ZZ & modulus_zz, struct BarretParams<T> & barret,
	const uint max_overflow);

template<typename T>
void freeBarret(struct BarretParams<T> & barret);

#endif
