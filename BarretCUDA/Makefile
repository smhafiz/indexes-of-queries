## This file is part of BarretCUDA v0.1 
## 
## BarretCUDA is a fast(ish) CUDA implementation of sparse matrix
## multiplication modulo a multi-precision prime.
## 
## Copyright (C) 2016, Ryan Henry and Syed Mahbub Hafiz
## 
## 
## BarretCUDA is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License,
## or (at your option) any later version.
## 
## BarretCUDA is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with BarretCUDA.  If not, see <http://www.gnu.org/licenses/>.

LIMBS=4

NVCC=nvcc
CXXFLAGS=-Wall -g -O2 -pedantic -std=c++11
NVCCFLAGS=-std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__
LDLIBS=-lntl -lgmp

genuint: genuint.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

rungenuint: genuint
	./genuint $(LIMBS) > uintX.h

barret: barret.cu rungenuint uint.h
	$(NVCC) $(NVCCFLAGS) -o $@ barret.cu $(LDLIBS)

clean:
	-rm -f barret genuint

veryclean: clean
	-rm -f uintX.h