## batchcoding Copyright 2016
## Syed Mahbub Hafiz <shafiz@indiana.edu>,
## Ryan Henry <henry@indiana.edu>
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of version 2 of the GNU General Public License as
## published by the Free Software Foundation.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## There is a copy of the GNU General Public License in the COPYING file
## packaged with this plugin; if you cannot find it, write to the Free
## Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301 USA

CXXFLAGS=-Wall -g -O2 -pedantic -std=c++11
NVCCFLAGS=-std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__

queryindex: vectormatrixcompcudactreams.cu
	nvcc $(NVCCFLAGS) -o $@ $^

createindex: createindex.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ -lntl -lgmp -lpthread -L/usr/local/lib -I/usr/local/lib

createindex_gf28: createindex_gf28.cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ gf2e.cc -L/usr/local/lib -I/usr/local/lib

clean:
	-rm -f createindex queryindex
