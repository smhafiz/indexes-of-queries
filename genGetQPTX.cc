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
	for (int i = 0; i < 2 * n; i++) {
		std::cout << i << " contains ";
		for (int j = 0; j < lows[i].size(); j++) {
			std::cout << lows[i][j].first << ", " << lows[i][j].second << "\t\t";
		}
		std::cout << std::endl;
	}
	bool pendingHigh = false;
	for (int i = 0; i < 2 * n; i++) {
		bool removeHigh = false;
		while (!lows[i].empty()) {
			if (occupied[i] == false) {
				
				if (lows[i].back().second < (n - 1)) {
					std::cout << "\"mul.lo.u32\t%" << i << ", %" << 2*n+lows[i].back().first << ", %" << 3*n+lows[i].back().second << ";\\n\\t\"";
					std::cout << "\t\t\\\\(O" << i << " = r" << lows[i].back().first << " * s" << lows[i].back().second << ".lo)" << std::endl;
					carries[i + 1] = false;
				} else {
					std::cout << "\"mov.lo.u32\t%" << i << " ,  %" << 2*n+lows[i].back().first << ";\\n\\t\"";
					std::cout << "\t\t\\\\(O" << i << " =  r" << lows[i].back().first << ")\n";
					removeHigh = true;
					carries[i + 1] = false;
				}
				occupied[i] = true;

			} else {
				if (lows[i].back().second < (n - 1)) {
					std::cout << "\"mad.lo.cc.u32\t%"<< i << ", %" << 2*n+lows[i].back().first << ", %" << 3*n+lows[i].back().second << ", %" << i << ";\\n\\t\"";
					std::cout << "\t\t\\\\(O" << i << " += r" << lows[i].back().first << " * s" << lows[i].back().second << ".lo" << ")\n";
					carries[i + 1] = true;
				} else {
					
					if(carries[i]) {
						std::cout<< "\"addc.lo.cc.u32\t%" << i << ", %" << 2*n+lows[i].back().first << ", %" << i << ";\\n\\t\"";
						std::cout<< "\t\t\\\\(O" << i << " += r" << lows[i].back().first << " + c" << ")\n";
					} else {
						std::cout<< "\"add.lo.cc.u32\t%" << i << ", %" << 2*n+lows[i].back().first << ", %" << i << ";\\n\\t\"";
						std::cout<< "\t\t\\\\(O" << i << " += r" << lows[i].back().first << ")\n";
					}
					carries[i + 1] = true;
					removeHigh = true;
					
				}

			}
			lows[i].pop_back();
			if (!removeHigh) {
				if (carries[i + 1] == true && occupied[i + 1] == true) {
					if (highs[i].back().second < (n - 1)) {
						std::cout << "\"madc.hi.cc.u32\t%" << i + 1 << ", %" << 2*n+highs[i].back().first << ", %" << 3*n+highs[i].back().second << ", %" << i+1 << ";\\n\\t\"";
						std::cout << "\t\t\\\\(O" << i + 1 << " += r" << highs[i].back().first << " * s" << highs[i].back().second << ".hi + c" << ")\n";
						carries[i + 2] = true;
					} else {
						std::cout << "\"addc.hi.cc.u32\t%" << i + 1 << ", %" << 2*n+highs[i].back().first << ", %" << i+1 << ";\\n\\t\"";
						std::cout << "\t\t\\\\(O" << i + 1 << " += r" << highs[i].back().first << " + c" << ")\n";
						carries[i + 2] = false;
					}
				} else if (carries[i + 1] == true && occupied[i + 1] == false) {
					if (highs[i].back().second < (n - 1)) {
						std::cout << "\"madc.hi.u32\t%" << i + 1 << ", %" << 2*n+highs[i].back().first << ", %" << 3*n+highs[i].back().second << ", " << 0 << ";\\n\\t\"";
						std::cout << "\t\t\\\\(O" << i + 1 << " = r" << highs[i].back().first << " * s" << highs[i].back().second << ".hi + c" << ")\n";
					} else {
						std::cout << "\"addc.hi.u32\t%" << i + 1 << ", %" << 2*n+highs[i].back().first << ", " << 0 << ";\\n\\t\"";
						std::cout << "\t\t\\\\(O" << i + 1 << " = r" << highs[i].back().first << " + c" << ")\n";
					}
					occupied[i + 1] = true;
				} else if (carries[i + 1] == false && occupied[i + 1] == true) {
					if (highs[i].back().second < (n - 1)) {
						std::cout << "\"mad.hi.cc.u32\t%" << i + 1 << ", %" << 2*n+highs[i].back().first << ", %" << 3*n+highs[i].back().second << ", %" << i+1 << ";\\n\\t\"";
						std::cout << "\t\t\\\\(O" << i + 1 << " += r" << highs[i].back().first << " * s" << highs[i].back().second << ".hi " << ")\n";
						carries[i + 2] = true;
					} else {
						std::cout << "\"add.hi.u32\t%" << i + 1 << ", %" << 2*n+highs[i].back().first << ", %" << i+1 << ";\\n\\t\"";
						std::cout << "\t\t\\\\(O" << i + 1 << " += r" << highs[i].back().first << ")\n";
						carries[i + 2] = false;
					}
				} else if (carries[i + 1] == false && occupied[i + 1] == false) {
					
					if (highs[i].back().second < (n - 1)) {
						std::cout << "\"mul.hi.u32\t%" << i + 1 << ", %" << 2*n+highs[i].back().first << ", %" << 3*n+highs[i].back().second << ";\\n\\t\"";
						std::cout << "\t\t\\\\(O" << i + 1 << " = r" << highs[i].back().first << " * s" << highs[i].back().second << ".hi" << ")\n";
					} else {
						std::cout << "\"mov.hi.u32\t%" << i + 1 << ", %" << 2*n+highs[i].back().first << ";\\n\\t\"";
						std::cout << "\t\t\\\\(O" << i + 1 << " = r" << highs[i].back().first << ")\n";
					}
					occupied[i + 1] = true;
				}
			}
			highs[i].pop_back();
			if (removeHigh) {
				removeHigh = false;
				break;
			}

			if ((i + 2) < 2 * n && (carries[i + 2] == true)) {
				if (!lows[i + 2].empty() && !pendingHigh) {
					if (occupied[i + 2] == false) {
						if (lows[i + 2].back().second < (n - 1)) {
							std::cout << "\"madc.lo.u32\t%" << i + 2 << ", %" << 2*n+lows[i + 2].back().first << ", %" << 3*n+lows[i + 2].back().second << ", " << 0 << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 2 << " = r" << lows[i + 2].back().first << " * s" << lows[i + 2].back().second << ".lo + c  " << ")\n";
						} else {
							std::cout << "\"addc.lo.u32\t%" << i + 2 << ", %" << 2*n+lows[i + 2].back().first << ", " << 0 << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 2 << " = r" << lows[i + 2].back().first << " + c  " << ")\n";
							highs[i + 2].pop_back();
						}
						occupied[i + 2] = true;
					} else {
						
						if (lows[i + 2].back().second < (n - 1)) {
							std::cout << "\"madc.lo.cc.u32\t%" << i + 2 << ", %" << 2*n+lows[i + 2].back().first << ", %" << 3*n+lows[i + 2].back().second << ", %" << i+2 << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 2 << " += r" << lows[i + 2].back().first << " * s" << lows[i + 2].back().second << ".lo + c  " << ")\n";
							std::cout << "\"addc.lo.u32\t%" << i + 3 << ", 0, 0" << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 3 << " = +c" << ")\n";
						} else {
							std::cout << "\"addc.lo.cc.u32\t%" << i + 2 << ", %" << 2*n+lows[i + 2].back().first << ", %" << i+2 << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 2 << " += r" << lows[i + 2].back().first << " + c  " << ")\n";
							std::cout << "\"addc.lo.u32\t%" << i + 3 << ", 0, 0" << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 3 << " = +c" << ")\n";

						}
						if (i + 3 < 2 * n)
							occupied[i + 3] = true;

					}
					pendingHigh = true;
					lows[i + 2].pop_back();
				} else if (!highs[i + 1].empty() && pendingHigh) {
					if (occupied[i + 2] == false) {
						
						if (highs[i + 1].back().second < (n - 1)) {
							std::cout << "\"madc.hi.u32\t%" << i + 2 << ", %" << 2*n+highs[i + 1].back().first << ", %" << 3*n+highs[i + 1].back().second << ", "<<0 << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 2 << " = r" << highs[i + 1].back().first << " * s" << highs[i + 1].back().second << ".hi + c   " << ")\n";
						} else {
							std::cout << "\"addc.hi.u32\t%" << i + 2 << " = r" << 2*n+highs[i + 1].back().first << ", " << 0 << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 2 << " = r" << highs[i + 1].back().first << " + c  "<< ")\n";
						}
						occupied[i + 2] = true;
					} else {
						
						if (highs[i + 1].back().second < (n - 1)) {
							std::cout << "\"madc.hi.u32\t%" << i + 2 << ", %" << 2*n+highs[i + 1].back().first << ", %" << 3*n+highs[i + 1].back().second << ", %" << i+2 << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 2 << " += r" << highs[i + 1].back().first << " * s" << highs[i + 1].back().second << ".hi + c   " << ")\n";
						} else {
							std::cout << "\"addc.hi.u32\t%" << i + 2 << ", %" << 2*n+highs[i + 1].back().first << ", %" << i+2 << ";\\n\\t\"";
							std::cout << "\t\t\\\\(O" << i + 2 << " += r" << highs[i + 1].back().first << ")\n";
						}

						if (i + 3 < 2 * n) {
							if (occupied[i + 3] == true) {
								std::cout << "\"addc.hi.u32\t%" << i + 3 << ", %" << i+3 << ", 0" << ";\\n\\t\"";
								std::cout << "\t\t\\\\(O" << i + 3 << " += +c" << ")\n";
							}
							else {
								std::cout << "\"addc.hi.u32\t%" << i + 3 << ", 0, 0" << ";\\n\\t\"";
								std::cout << "\t\t\\\\(O" << i + 3 << " = +c" << ")\n";
								occupied[i + 3] = true;
							}
						}

					}
					pendingHigh = false;
					highs[i + 1].pop_back();
				}
			}
		std::fill(carries.begin(), carries.end(), false);
		}
		std::cout << std::endl << std::endl;
	}
	if(n>1) {
		std::cout<<"O"<<(n%3==0?2*n-2:2*n-1)<<(occupied[(n%3==0?2*n-2:2*n-1)]?" +":" ")<<"= +c\n";
	}
	if(n>2 && n%3==0) { 
		std::cout<<"\"addc.u32\t%"<<2*n-1<<", 0,"<< (occupied[(2*n-1)]?" %":" ") << (occupied[(2*n-1)]?2*n-1:0) <<";\\n\\t\"";
		std::cout<<"\t\t\\\\(O"<<2*n-1<< (occupied[(2*n-1)]?" +":" ") <<" =+c"<<")\n";
	}
	myfile.close();
}
