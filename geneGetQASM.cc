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
#include <vector>

#define REGISTER 2
#define CARRY 1
#define FREE 0

void printPairs(std::vector<std::vector<std::pair<int, int>>> sample, int n);

int main(int argc, char **argv){
	int n = atoi(argv[1]);
	std::vector<std::vector<std::pair<int, int>>>lows(2 * n);
	std::vector<std::vector<std::pair<int, int>>>highs(2 * n);
	int occupiedBy[2*n];
	int firstOperandStart = n+1;
	int secondOperandStart = firstOperandStart + n;
	int sMaxReq = secondOperandStart+n-1;

	for (int i = 0; i < 2 * n; ++i) {
		occupiedBy[i] = FREE;
		for (int j = 0, a = 0; j <= i; j++, a++) {
			if (i - j < n   && j < n) {
				lows[i].push_back(std::make_pair(firstOperandStart+j,secondOperandStart+ i - j));
				highs[i].push_back(std::make_pair(firstOperandStart+j,secondOperandStart+ i - j));
			}
		}
	}
	printPairs(lows,n);
	
	auto q = highs[0].back();
	std::cout<<"\"mul.hi.u32     %0, %"<<q.first<<", %"<<q.second<<";\\n\\t\"" << "\t\t//(O0 = r"<<q.first<<"  * s"<<q.second<<".hi)\n";
	occupiedBy[0] = REGISTER;
	occupiedBy[1] = REGISTER;
	lows[0].pop_back();
	highs[0].pop_back();

	for (int i = 1; i < 2*n ; i++) {
		while (!lows[i].empty()) {
			int l = i;
			bool lo = true;
			bool carry = false;
			while(l< 2*n){
				if(lo && !lows[l].empty()){
					auto p = lows[l].back();
					if(occupiedBy[l]){
						if(carry) {
							if(p.second < sMaxReq) {
								std::cout << "\"madc.lo.cc.u32\t%"  << (l - 1)%(n+1) << ", %" << p.first << ", %" << p.second << ", %" << (l - 1)%(n+1) << ";\\n\\t\"";
								std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " += r" << p.first << " * s" << p.second << ".lo + c  " << ")\n";
								lo = false;
							} else {
								std::cout<< "\"addc.cc.u32\t%"  << (l - 1)%(n+1) << ", %" << p.first << ", %"  << (l - 1)%(n+1) << ";\\n\\t\"";
								std::cout<< "\t\t//(O"  << (l - 1)%(n+1) << " += r" << p.first << " + c" << ")\n";

								int g=0;
								while(!highs[l].empty()) {
									if((highs[l][g].first==p.first) && (highs[l][g].second = p.second)) {
										highs[l].erase(highs[l].begin()+g);
										break;
									}
									g++;
								}
								if(highs[l].empty()) {lo = true;} else {lo = false;}//TODO
							}

						} else {
							if(p.second < sMaxReq) {
								std::cout << "\"mad.lo.cc.u32\t%" << (l - 1)%(n+1) << ", %" << p.first << ", %" << p.second << ", %"  << (l - 1)%(n+1) << ";\\n\\t\"";
								std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " += r" << p.first << " * s" << p.second << ".lo" << ")\n";
								lo = false;
							} else {
								std::cout<< "\"add.cc.u32\t%"  << (l - 1)%(n+1) << ", %" << p.first << ", %"  << (l - 1)%(n+1) << ";\\n\\t\"";
								std::cout<< "\t\t//(O"  << (l - 1)%(n+1) << " += r" << p.first << ")\n";

								int g=0;
								while(!highs[l].empty()) {
									if((highs[l][g].first==p.first) && (highs[l][g].second = p.second)) {
										highs[l].erase(highs[l].begin()+g);
										break;
									}
									g++;
								}
								if(highs[l].empty()) {lo = true;} else {lo = false;}//TODO
							}

						}
						occupiedBy[l] = REGISTER;
						carry = true;
						lows[l].pop_back();
					} else {
						if(p.second < sMaxReq) {
							std::cout << "\"madc.lo.u32\t%"  << (l - 1)%(n+1) << ", %" << p.first << ", %" << p.second << ", " << 0 << ";\\n\\t\"";
							std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " = r" << p.first << " * s" << p.second << ".lo + c  " << ")\n";
							lo = false;
						} else {
							std::cout << "\"addc.u32\t%"  << (l - 1)%(n+1) << ", %" << p.first << ", " << 0 << ";\\n\\t\"";
							std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " = r" << p.first << " + c  " << ")\n";
							int g=0;
							while(!highs[l].empty()) {
								if((highs[l][g].first==p.first) && (highs[l][g].second = p.second)) {
									highs[l].erase(highs[l].begin()+g);
									break;
								}
								g++;
							}
							lo = true;
						}
						occupiedBy[l] = REGISTER;
						carry = false;
						lows[l].pop_back();

						break;
					}

				} else if(!highs[l - 1].empty()) {
					auto p = highs[l - 1].back();
					if(occupiedBy[l]){
						if(carry) {
							std::cout << "\"madc.hi.cc.u32\t%"  << (l - 1)%(n+1)<< ", %" << p.first << ", %" << p.second << ", %" << (l - 1)%(n+1) << ";\\n\\t\"";
							std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " += r" << p.first << " * s" << p.second << ".hi + c  " << ")\n";

						} else {
							std::cout << "\"mad.hi.cc.u32\t%" << (l - 1)%(n+1) << ", %" << p.first << ", %" << p.second << ", %"  << (l - 1)%(n+1) << ";\\n\\t\"";
							std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " += r" << p.first << " * s" << p.second << ".hi" << ")\n";

						}
						occupiedBy[l] = REGISTER;
						carry = true;
						highs[l - 1].pop_back();
						lo = true;

					} else {
						if(carry){
							std::cout << "\"madc.hi.u32\t%"  << (l - 1)%(n+1) << ", %" << p.first << ", %" << p.second << ", " << 0 << ";\\n\\t\"";
							std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " = r" << p.first << " * s" << p.second << ".hi + c  " << ")\n";
						}
						occupiedBy[l] = REGISTER;
						carry = false;
						highs[l - 1].pop_back();
						lo = true;
						break;
					}
				} else {
					if(carry) {
						if(occupiedBy[l]==REGISTER) {
							std::cout << "\"addc.cc.u32\t%"  << (l - 1)%(n+1) << ", 0, %" << (l - 1)%(n+1) << ";\\n\\t\"";
							std::cout << "\t\t//(O"  << (l - 1)%(n+1) << "  + = +c" << ")\n";
							carry = true;
						} else if(occupiedBy[l]==CARRY) {
							std::cout << "\"addc.u32\t%"  << (l - 1)%(n+1) << ", 0, %" << (l - 1)%(n+1) << ";\\n\\t\"";
							std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " + = +c" << ")\n";
							carry = false;
						} else if(occupiedBy[l]==FREE) {
							std::cout << "\"addc.u32\t%"  << (l - 1)%(n+1) << ", 0, " << 0 << ";\\n\\t\"";
							std::cout << "\t\t//(O"  << (l - 1)%(n+1) << " = +c" << ")\n";
							occupiedBy[l] = CARRY;
							carry = false;
						}
					}
				}
				l++;
			}
		std::cout << "\n";
		}
//		std::cout << "************column" << i << "*********************\n";
	}
	printPairs(lows,n);
}

void printPairs(std::vector<std::vector<std::pair<int, int>>> sample, int n){
	std::cout<<"\nPrinting pairs:\n" ;
	for (int i = 0; i < 2 * n; i++) {
		std::cout << i << " contains ";
		for (int j = 0; j < sample[i].size(); j++) {
			std::cout << sample[i][j].first << ", " << sample[i][j].second << "\t\t";
		}
		std::cout << std::endl;
	}
	std::cout<<"\n";
}
