#include <stdlib.h>
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

using namespace std;

#define CARRY 1
#define ADD   2
#define MUL   2
#define MAD   2
#define ADDC  3
#define MADC  3

#define CARRY_IN(i)	(state[i] |= CARRY)
#define ADD_IN(i)	(state[i] |= ADD)
#define ADDC_IN(i)	(state[i] |= ADDC)
#define MUL_IN(i)	(state[i] |= MUL)
#define MAD_IN(i)	(state[i] |= MAD)
#define MADC_IN(i)	(state[i] |= MADC)

#define IS_OCCUPIED(i)	(state[i] ? true : false)
#define DO_CARRY_OUT(i)	(state[i] & ADD ? true : false)
#define CARRY_OUT_FLAG(i)(state[i] & ADD ? ".cc" : "")
#define CARRY_IN_FLAG(b)(b ? "c" : "")
#define LO_OR_HI(p)	(p.first ? ".hi" : ".lo")

#define GET_FIRST_REG(p) (2*n-1+p.second)
#define GET_SECOND_REG(i,p) (3*n-1+i-p.second-p.first)

void printPairs(std::vector<std::vector<std::pair<int, int>>> * sample, int n);

std::vector<std::vector<std::pair<int,int>>> *lows;
int * state;
char * buf = new char[10];

void propagate(int i, int n, bool cin)
{
	bool cc = DO_CARRY_OUT(i);
	if (!cin) cout << "\n";
	if ((*lows)[i].empty())
	{
		sprintf(buf, "%%%u", i);
		printf("\"add%s%s.u32\t%%%u, %s, 0;\\n\\t\"\n", CARRY_IN_FLAG(cin), CARRY_OUT_FLAG(i), i, (IS_OCCUPIED(i) ? buf : "0"));
		ADDC_IN(i);
	}
	else
	{
		auto p = (*lows)[i].back();
		(*lows)[i].pop_back();
		sprintf(buf, "%%%u", i);
		printf("\"mad%s%s%s.u32\t%%%u, %%%u, %%%u, %s;\\n\\t\"\n", CARRY_IN_FLAG(cin), LO_OR_HI(p), CARRY_OUT_FLAG(i), i, GET_FIRST_REG(p), GET_SECOND_REG(i,p), (IS_OCCUPIED(i) ? buf : "0"));
		MADC_IN(i);
	}
	if (cc && i < 2*n-2) propagate(i+1,n,true);
}

int main(int argc, char **argv){
	int n = atoi(argv[1]);
	state = new int[2*n];
	lows = new std::vector<std::vector<std::pair<int,int>>>();
	lows->resize(2*n);

	for (int i = 0; i < 2 * n; ++i) {
		for (int j = 0, a = 0; j <= i; j++, a++) {
			if (i - j < n - 1 && j < n) {
			//if (i - j < n && j < n) {
				(*lows)[i].push_back(std::make_pair(0, j));
				(*lows)[i+1].push_back(std::make_pair(1, j));
			}
		}
	}
	printPairs(lows,n);
	bool pendingHigh = false;

	auto p = (*lows)[0].back();
	(*lows)[0].pop_back();
	lows->erase(lows->begin());
	std::cout<<"\"mul.hi.u32     %0, %"<<GET_FIRST_REG(p)<<", %"<<GET_SECOND_REG(0,p)<<";\\n\\t\"" << "\t\t//(O0 = r0 * s0.hi)\n";
	MUL_IN(0);

	for (int i = 0; i < 2*n-1; i++)
	{
		while (!(*lows)[i].empty()) {
			propagate(i,n,false);
		}
	}



/*
	for (int i = 0; i < 2 * n ; i++) {
		while (!lows[i].empty()) {
			if (occupied[i] == false ) {
				if (lows[i].back().second < (n - 1)) {
					std::cout << "\"mul.lo.u32\t%"  << i -1 << ", %" << 2*n-1+lows[i].back().first << ", %" << 3*n-1+lows[i].back().second << ";\\n\\t\"";
					std::cout << "\t\t//(O"  << i -1 << " = r" << lows[i].back().first << " * s" << lows[i].back().second << ".lo)" << std::endl;
					//carries[i + 1] = false;
				} else {
					std::cout << "\"mov.u32\t%"  << i -1 << " ,  %" << 2*n-1+lows[i].back().first << ";\\n\\t\"";
					std::cout << "\t\t//(O"  << i -1 << " =  r" << lows[i].back().first << ")\n";
					//carries[i + 1] = false;
					int g=0;
					while(!highs[i].empty()) {
						if((highs[i][g].first==lows[i].back().first) && (highs[i][g].second = lows[i].back().second)) {
							highs[i].erase(highs[i].begin()+g);
							break;
						}
						g++;
					}
				}
				occupied[i] = true;

			} else {
				if (lows[i].back().second < (n - 1)) {
					std::cout << "\"mad.lo.cc.u32\t%" << i -1 << ", %" << 2*n-1+lows[i].back().first << ", %" << 3*n-1+lows[i].back().second << ", %"  << i -1 << ";\\n\\t\"";
					std::cout << "\t\t//(O"  << i -1 << " += r" << lows[i].back().first << " * s" << lows[i].back().second << ".lo" << ")\n";
					carries[i + 1] = true;
				} else {

					if(carries[i]) {
						std::cout<< "\"addc.cc.u32\t%"  << i -1 << ", %" << 2*n-1+lows[i].back().first << ", %"  << i -1 << ";\\n\\t\"";
						std::cout<< "\t\t//(O"  << i -1 << " += r" << lows[i].back().first << " + c" << ")\n";
					} else {
						std::cout<< "\"add.cc.u32\t%"  << i -1 << ", %" << 2*n-1+lows[i].back().first << ", %"  << i -1 << ";\\n\\t\"";
						std::cout<< "\t\t//(O"  << i -1 << " += r" << lows[i].back().first << ")\n";
					}
					carries[i + 1] = true;
					int g=0;
					while(!highs[i].empty()) {
						if((highs[i][g].first==lows[i].back().first) && (highs[i][g].second = lows[i].back().second)) {
							highs[i].erase(highs[i].begin()+g);
							break;
						}
						g++;
					}

				}

			}
			lows[i].pop_back();
			if(!highs[i].empty()){
				if (carries[i + 1] == true && occupied[i + 1] == true) {
					if (highs[i].back().second < (n - 1)) {
						std::cout << "\"madc.hi.cc.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ", %" << 3*n-1+highs[i].back().second << ", %" << i-1+1 << ";\\n\\t\"";
						std::cout << "\t\t//(O"  << i -1 + 1 << " += r" << highs[i].back().first << " * s" << highs[i].back().second << ".hi + c" << ")\n";
						carries[i + 2] = true;
					} else {
						std::cout << "\"addc.cc.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ", %" << i-1+1 << ";\\n\\t\"";
						std::cout << "\t\t//(O"  << i -1 + 1 << " += r" << highs[i].back().first << " + c" << ")\n";
						carries[i + 2] = true;//TODO
						int g=0;
						while(!lows[i].empty()) {
							if((lows[i][g].first==highs[i].back().first) && (lows[i][g].second = highs[i].back().second)) {
								lows[i].erase(lows[i].begin()+g);
								break;
							}
							g++;
						}

					}
				} else if (carries[i + 1] == true && occupied[i + 1] == false) {
					if (highs[i].back().second < (n - 1)) {
						std::cout << "\"madc.hi.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ", %" << 3*n-1+highs[i].back().second << ", " << 0 << ";\\n\\t\"";
						std::cout << "\t\t//(O"  << i -1 + 1 << " = r" << highs[i].back().first << " * s" << highs[i].back().second << ".hi + c" << ")\n";
					} else {
						std::cout << "\"1addc.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ", " << 0 << ";\\n\\t\"";
						std::cout << "\t\t//(O"  << i -1 + 1 << " = r" << highs[i].back().first << " + c" << ")\n";
						int g=0;
						while(!lows[i].empty()) {
							if((lows[i][g].first==highs[i].back().first) && (lows[i][g].second = highs[i].back().second)) {
								lows[i].erase(lows[i].begin()+g);
								break;
							}
							g++;
						}
					}
					//carries[i + 2] = false;//TODO
					occupied[i + 1] = true;
				}
				highs[i].pop_back();
				if ((i + 2) < 2 * n && (carries[i + 2] == true)) {
					bool done = false;
					while(!done){
						if(!lows[i + 2].empty() && !pendingHigh) {
							if (occupied[i + 2] == false) {
								if (lows[i + 2].back().second < (n - 1)) {
									std::cout << "\"madc.lo.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+lows[i + 2].back().first << ", %" << 3*n-1+lows[i + 2].back().second << ", " << 0 << ";\\n\\t\"";
									std::cout << "\t\t//(O"  << i -1 + 2 << " = r" << lows[i + 2].back().first << " * s" << lows[i + 2].back().second << ".lo + c  " << ")\n";
								} else {
									std::cout << "\"2addc.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+lows[i + 2].back().first << ", " << 0 << ";\\n\\t\"";
									std::cout << "\t\t//(O"  << i -1 + 2 << " = r" << lows[i + 2].back().first << " + c  " << ")\n";
									int g=0;
									while(!highs[i+2].empty()) {
										if((highs[i+2][g].first==lows[i+2].back().first) && (highs[i+2][g].second = lows[i+2].back().second)) {
											highs[i+2].erase(highs[i+2].begin()+g);
											break;
										}
										g++;
									}
								}
								occupied[i + 2] = true;
							} else {
								if (lows[i + 2].back().second < (n - 1)) {
									std::cout << "\"madc.lo.cc.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+lows[i + 2].back().first << ", %" << 3*n-1+lows[i + 2].back().second << ", %" << i-1+2 << ";\\n\\t\"";
									std::cout << "\t\t//(O"  << i -1 + 2 << " += r" << lows[i + 2].back().first << " * s" << lows[i + 2].back().second << ".lo + c  " << ")\n";
									std::cout << "\"3addc.u32\t%"  << i -1 + 3 << ", 0, 0" << ";\\n\\t\"";
									std::cout << "\t\t//(O"  << i -1 + 3 << " = +c" << ")\n";
								} else {
									std::cout << "\"addc.cc.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+lows[i + 2].back().first << ", %" << i-1+2 << ";\\n\\t\"";
									std::cout << "\t\t//(O"  << i -1 + 2 << " += r" << lows[i + 2].back().first << " + c  " << ")\n";
									std::cout << "\"4addc.u32\t%"  << i -1 + 3 << ", 0, 0" << ";\\n\\t\"";
									std::cout << "\t\t//(O"  << i -1 + 3 << " = +c" << ")\n";
									int g=0;
									while(!highs[i+2].empty()) {
										if((highs[i+2][g].first==lows[i+2].back().first) && (highs[i+2][g].second = lows[i+2].back().second)) {
											highs[i+2].erase(highs[i+2].begin()+g);
											break;
										}
										g++;
									}

								}
								if (i + 3 < 2 * n)
									occupied[i + 3] = true;

							}
							done = true;
							pendingHigh = true;
							lows[i + 2].pop_back();
						} else if (!highs[i + 1].empty() && pendingHigh) {
							bool restart = false;
							int h=lows[i+1].size();
							while(h--) {
								if((lows[i+1][h].first==highs[i+1].back().first) && (lows[i+1][h].second = highs[i+1].back().second)) {
									done=false;
									restart = true;
									pendingHigh = false;
									break;
								}

							}
							if(!restart) {
								if (occupied[i + 2] == false) {

									if (highs[i + 1].back().second < (n - 1)) {
										std::cout << "\"madc.hi.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+highs[i + 1].back().first << ", %" << 3*n-1+highs[i + 1].back().second << ", "<<0 << ";\\n\\t\"";
										std::cout << "\t\t//(O"  << i -1 + 2 << " = r" << highs[i + 1].back().first << " * s" << highs[i + 1].back().second << ".hi + c   " << ")\n";
									} else {
										std::cout << "\"5addc.u32\t%"  << i -1 + 2 << " = r" << 2*n-1+highs[i + 1].back().first << ", " << 0 << ";\\n\\t\"";
										std::cout << "\t\t//(O"  << i -1 + 2 << " = r" << highs[i + 1].back().first << " + c  "<< ")\n";
										int g=0;
										while(!lows[i+1].empty()) {
											if((lows[i+1][g].first==highs[i+1].back().first) && (lows[i+1][g].second = highs[i+1].back().second)) {
												lows[i+1].erase(lows[i+1].begin()+g);
												break;
											}
											g++;
										}
									}
									occupied[i + 2] = true;
								} else {

									if (highs[i + 1].back().second < (n - 1)) {
										std::cout << "\"madc.hi.u32\t%"  << i -1+ 2 << ", %" << 2*n-1+highs[i + 1].back().first << ", %" << 3*n-1+highs[i + 1].back().second << ", %" << i-1+2 << ";\\n\\t\"";
										std::cout << "\t\t//(O"  << i -1 + 2 << " += r" << highs[i + 1].back().first << " * s" << highs[i + 1].back().second << ".hi + c   " << ")\n";
									} else {
										std::cout << "\"6addc.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+highs[i + 1].back().first << ", %" << i-1+2 << ";\\n\\t\"";
										std::cout << "\t\t//(O"  << i -1 + 2 << " += r" << highs[i + 1].back().first << ")\n";
										int g=0;
										while(!lows[i+1].empty()) {
											if((lows[i+1][g].first==highs[i+1].back().first) && (lows[i+1][g].second = highs[i+1].back().second)) {
												lows[i+1].erase(lows[i+1].begin()+g);
												break;
											}
											g++;
										}
									}

									if (i + 3 < 2 * n) {
										if (occupied[i + 3] == true) {
											std::cout << "\"7addc.u32\t%"  << i -1 + 3 << ", %" << i-1+3 << ", 0" << ";\\n\\t\"";
											std::cout << "\t\t//(O"  << i -1 + 3 << " += +c" << ")\n";
										}
										else {
											std::cout << "\"8addc.u32\t%"  << i -1 + 3 << ", 0, 0" << ";\\n\\t\"";
											std::cout << "\t\t//(O"  << i -1 + 3 << " = +c" << ")\n";
											occupied[i + 3] = true;
										}
									}

								}
								done=true;
								pendingHigh = false;
								highs[i + 1].pop_back();
							}
						}
					}
				}
				std::fill(carries.begin(), carries.end(), false);
			} else {}

		}
		//std::cout << "************column" << i+1 << "*********************\n";
	}
	if(n>1) {
		/*O4 += +c
		addc.cc.u32 %4,0,%4;
		addc.u32 %5* /

		//		if(n%3==0){
		if(occupied[2*n-2])	std::cout<<"\"addc.cc.u32\t%"<<2*n-2-1<<", 0, %" << 2*n-2-1<< ";\\n\\t\"";
		else std::cout<<"\"addc.u32\t%"<<2*n-2-1<<", 0, %" << 0 << ";\\n\\t\"";;
		std::cout<<"\t\t//(O"<<2*n-2-1<<(occupied[2*n-2]?" +":" ")<<"= +c"<<")\n";
		//		} else {
		//			if(occupied[2*n-1]) std::cout<<"\"addc.u32\t%"<<2*n-1-1<< ", 0, %" << 2*n-1-1 << ";\\n\\t\"";
		//			else std::cout<<"\"addc.u32\t%"<<2*n-1-1<< ", 0, %" << 0 << ";\\n\\t\"";
		//			std::cout<<"\t\t//(O"<<2*n-1-1<<(occupied[2*n-1]?" +":" ")<<"= +c"<<")\n";
		//		}
	}
	if(n>2) {// && n%3==0
		std::cout<<"\"addc.u32\t%"<<2*n-1-1<<", 0,"<< (occupied[(2*n-1)]?" %":" ") << (occupied[(2*n-1)]?2*n-1-1:0) <<";\\n\\t\"";
		std::cout<<"\t\t//(O"<<2*n-1-1<< (occupied[(2*n-1)]?" +":" ") <<" =+c"<<")\n";
	}
	myfile.close();
*/
}

void printPairs(std::vector<std::vector<std::pair<int, int>>> * sample, int n){
	std::cout<<"\nPrinting pairs:\n" ;
	for (int i = 0; i < 2 * n; i++) {
		std::cout << i << " contains ";
		for (int j = 0; j < (*sample)[i].size(); j++) {
			auto p = (*sample)[i][j];
			std::cout << "(" << GET_FIRST_REG(p) << ", " << GET_SECOND_REG(i,p) << ")." << (p.first ? "hi" : "lo") << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<<"\n";
}

/*else if (carries[i + 1] == false && occupied[i + 1] == true) {//TODO not need
					if (highs[i].back().second < (n - 1)) {
						std::cout << "\"DDDDDDDmad.hi.cc.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ", %" << 3*n-1+highs[i].back().second << ", %" << i-1+1 << ";\\n\\t\"";
						std::cout << "\t\t//(O"  << i -1 + 1 << " += r" << highs[i].back().first << " * s" << highs[i].back().second << ".hi " << ")\n";
						carries[i + 2] = true;
					} else {
						std::cout << "\"DDDDDDDDDDDadd.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ", %" << i-1+1 << ";\\n\\t\"";
						std::cout << "\t\t//(O"  << i -1 + 1 << " += r" << highs[i].back().first << ")\n";
						carries[i + 2] = true;//TODO
						int g=0;
						while(!lows[i].empty()) {
							if((lows[i][g].first==highs[i].back().first) && (lows[i][g].second = highs[i].back().second)) {
								lows[i].erase(lows[i].begin()+g);
								break;
							}
							g++;
						}
					}
				} else if (carries[i + 1] == false && occupied[i + 1] == false) {//TODO not need

					if (highs[i].back().second < (n - 1)) {
						std::cout << "\"EEEEEEEEEEEmul.hi.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ", %" << 3*n-1+highs[i].back().second << ";\\n\\t\"";
						std::cout << "\t\t//(O"  << i -1 + 1 << " = r" << highs[i].back().first << " * s" << highs[i].back().second << ".hi" << ")\n";
					} else {
						std::cout << "\"EEEEEEEEEEmov.hi.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ";\\n\\t\"";
						std::cout << "\t\t//(O"  << i -1 + 1 << " = r" << highs[i].back().first << ")\n";
						int g=0;
						while(!lows[i].empty()) {
							if((lows[i][g].first==highs[i].back().first) && (lows[i][g].second = highs[i].back().second)) {
								lows[i].erase(lows[i].begin()+g);
								break;
							}
							g++;
						}
					}
					carries[i + 2] = false;//TODO
					occupied[i + 1] = true;

				}*/
