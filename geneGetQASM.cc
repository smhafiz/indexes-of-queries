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

void printPairs(std::vector<std::vector<std::pair<int, int>>> sample, int n);

int main(int argc, char **argv){
	int n = atoi(argv[1]);
	std::vector<std::vector<std::pair<int, int>>>lows(2 * n);
	std::vector<std::vector<std::pair<int, int>>>highs(2 * n);
	std::vector<bool> carries;
	std::vector<bool> occupied;
	std::ofstream myfile;
	std::string fileName = "ptxCodeGCGQ" + std::to_string(n) + ".cc";
	myfile.open(fileName);

	for (int i = 0; i < 2 * n; ++i) {
		carries.push_back(false);
		occupied.push_back(false);
		for (int j = 0, a = 0; j <= i; j++, a++) {
			if (i - j < n && j < n) {
				lows[i].push_back(std::make_pair(j, i - j));
				highs[i].push_back(std::make_pair(j, i - j));
			}
		}
	}
	printPairs(lows,n);
	bool pendingHigh = false;

	//discarding the unnecessary operation cum register
	std::cout<<"\"mul.hi.u32     %0, %"<<2*n-1<<", %"<<3*n-1<<";\\n\\t\"" << "\t\t//(O0 = r0 * s0.hi)\n";
	occupied[0]=true;
	occupied[1]=true;
	lows[0].pop_back();
	highs[0].pop_back();
	//printPairs(highs,n);

	for (int i = 0; i < 2 * n ; i++) {
		while (!lows[i].empty()) {
			if (occupied[i] == false ) {
				if (lows[i].back().second < (n - 1)) {
					std::cout << "\"mul.lo.u32\t%"  << i -1 << ", %" << 2*n-1+lows[i].back().first << ", %" << 3*n-1+lows[i].back().second << ";\\n\\t\"";
					std::cout << "\t\t//(O"  << i -1 << " = r" << lows[i].back().first << " * s" << lows[i].back().second << ".lo)" << std::endl;
					carries[i + 1] = false;
				} else {
					std::cout << "\"mov.u32\t%"  << i -1 << " ,  %" << 2*n-1+lows[i].back().first << ";\\n\\t\"";
					std::cout << "\t\t//(O"  << i -1 << " =  r" << lows[i].back().first << ")\n";
					carries[i + 1] = false;
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
						std::cout << "\"addc.u32\t%"  << i -1 + 1 << ", %" << 2*n-1+highs[i].back().first << ", " << 0 << ";\\n\\t\"";
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
					carries[i + 2] = false;//TODO
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
									std::cout << "\"addc.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+lows[i + 2].back().first << ", " << 0 << ";\\n\\t\"";
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
									std::cout << "\"addc.u32\t%"  << i -1 + 3 << ", 0, 0" << ";\\n\\t\"";
									std::cout << "\t\t//(O"  << i -1 + 3 << " = +c" << ")\n";
								} else {
									std::cout << "\"addc.cc.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+lows[i + 2].back().first << ", %" << i-1+2 << ";\\n\\t\"";
									std::cout << "\t\t//(O"  << i -1 + 2 << " += r" << lows[i + 2].back().first << " + c  " << ")\n";
									std::cout << "\"addc.u32\t%"  << i -1 + 3 << ", 0, 0" << ";\\n\\t\"";
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
										std::cout << "\"addc.u32\t%"  << i -1 + 2 << " = r" << 2*n-1+highs[i + 1].back().first << ", " << 0 << ";\\n\\t\"";
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
										std::cout << "\"addc.u32\t%"  << i -1 + 2 << ", %" << 2*n-1+highs[i + 1].back().first << ", %" << i-1+2 << ";\\n\\t\"";
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
											std::cout << "\"addc.u32\t%"  << i -1 + 3 << ", %" << i-1+3 << ", 0" << ";\\n\\t\"";
											std::cout << "\t\t//(O"  << i -1 + 3 << " += +c" << ")\n";
										}
										else {
											std::cout << "\"addc.u32\t%"  << i -1 + 3 << ", 0, 0" << ";\\n\\t\"";
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
		addc.u32 %5*/

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
