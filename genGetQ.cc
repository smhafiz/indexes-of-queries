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
			std::cout << lows[i][j].first << ", " << lows[i][j].second
					<< "\t\t";
		}
		std::cout << std::endl;
	}
	bool pendingHigh = false;
	for (int i = 0; i < 2 * n; i++) {
		bool removeHigh = false;
		while (!lows[i].empty()) {
			if (occupied[i] == false) {
				std::cout << "O" << i << " = r" << lows[i].back().first;
				if (lows[i].back().second < (n - 1)) {
					std::cout << " * s" << lows[i].back().second << ".lo"
							<< std::endl;
					carries[i + 1] = false;
				} else {
					removeHigh = true;
					carries[i + 1] = false;
					std::cout << std::endl;
				}
				occupied[i] = true;
			} else {
				std::cout << "O" << i << " += r" << lows[i].back().first;
				if (lows[i].back().second < (n - 1)) {
					std::cout << " * s" << lows[i].back().second << ".lo" << std::endl;
					carries[i + 1] = true;
				} else {
					
					if(carries[i]) {std::cout<< " + c" << std::endl;}
					else {std::cout << std::endl;}
					carries[i + 1] = true;
					removeHigh = true;
					
				}

			}
			lows[i].pop_back();
			if (!removeHigh) {
				if (carries[i + 1] == true && occupied[i + 1] == true) {
					std::cout << "O" << i + 1 << " += r"
							<< highs[i].back().first;
					if (highs[i].back().second < (n - 1)) {
						std::cout << " * s" << highs[i].back().second
								<< ".hi + c" << std::endl;
						carries[i + 2] = true;
					} else {
						std::cout << " + c" << std::endl;
						carries[i + 2] = false;
					}

				} else if (carries[i + 1] == true && occupied[i + 1] == false) {
					std::cout << "O" << i + 1 << " = r"
							<< highs[i].back().first;
					if (highs[i].back().second < (n - 1)) {
						std::cout << " * s" << highs[i].back().second
								<< ".hi + c" << std::endl;
					} else {
						std::cout << " + c" << std::endl;
					}
					occupied[i + 1] = true;
				} else if (carries[i + 1] == false && occupied[i + 1] == true) {
					std::cout << "O" << i + 1 << " += r"
							<< highs[i].back().first;
					if (highs[i].back().second < (n - 1)) {
						std::cout << " * s" << highs[i].back().second << ".hi "
								<< std::endl;
						carries[i + 2] = true;
					} else {
						std::cout << std::endl;
						carries[i + 2] = false;
					}

				} else if (carries[i + 1] == false
						&& occupied[i + 1] == false) {
					std::cout << "O" << i + 1 << " = r"
							<< highs[i].back().first;
					if (highs[i].back().second < (n - 1)) {
						std::cout << " * s" << highs[i].back().second << ".hi"
								<< std::endl;
					} else {
						std::cout << std::endl;
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
						std::cout << "O" << i + 2 << " = r"
								<< lows[i + 2].back().first;
						if (lows[i + 2].back().second < (n - 1)) {
							std::cout << " * s" << lows[i + 2].back().second
									<< ".lo + c  " << std::endl;
						} else {
							highs[i + 2].pop_back();

							std::cout << " + c  " << std::endl;
						}
						occupied[i + 2] = true;
					} else {
						std::cout << "O" << i + 2 << " += r"
								<< lows[i + 2].back().first;
						if (lows[i + 2].back().second < (n - 1)) {
							std::cout << " * s" << lows[i + 2].back().second
									<< ".lo + c  " << std::endl;
							std::cout << "O" << i + 3 << " = +c" << "\n";
						} else {
							std::cout << " + c  " << std::endl;
							std::cout << "O" << i + 3 << " = +c" << "\n";

						}
						if (i + 3 < 2 * n)
							occupied[i + 3] = true;

					}
					pendingHigh = true;
					lows[i + 2].pop_back();
				} else if (!highs[i + 1].empty() && pendingHigh) {
					if (occupied[i + 2] == false) {
						std::cout << "O" << i + 2 << " = r"
								<< highs[i + 1].back().first;
						if (highs[i + 1].back().second < (n - 1))
							std::cout << " * s" << highs[i + 1].back().second
									<< ".hi + c   " << std::endl;
						else
							std::cout << std::endl;

						occupied[i + 2] = true;
					} else {
						std::cout << "O" << i + 2 << " += r"
								<< highs[i + 1].back().first;
						if (highs[i + 1].back().second < (n - 1))
							std::cout << " * s" << highs[i + 1].back().second
									<< ".hi + c   " << std::endl;
						else
							std::cout << std::endl;

						if (i + 3 < 2 * n) {
							if (occupied[i + 3] == true)
								std::cout << "O" << i + 3 << " += +c" << "\n";
							else {
								std::cout << "O" << i + 3 << " = +c" << "\n";
								occupied[i + 3] = true;
							}
						}

					}
					pendingHigh = false;
					highs[i + 1].pop_back();
				}
			}
		}
		std::cout << std::endl;
	}
	if(n>1)std::cout<<"O"<<(n%3==0?2*n-2:2*n-1)<<(occupied[(n%3==0?2*n-2:2*n-1)]?" +":" ")<<"= +c\n";
	if(n>2) std::cout<<"O"<<(n%3==0?2*n-1:2*n)<<(n%3==0?" +":" ")<<" =+c\n";
	myfile.close();
}

		/*std::cout<< "********************************************************************\n";
		 for(int e=0;e<2*n;e++){
		 std::cout<< e << " contains ";
		 for(int f=0;f<lows[e].size();f++){
		 std::cout<< "L" <<lows[e][f].first << ", " << lows[e][f].second << "\t";
		 if(f < highs[e].size()) std::cout<< "H" << highs[e][f].first << ", " << highs[e][f].second << "\t\t";
		 }
		 std::cout<< std::endl;
		 }
		 std::cout<< "********************************************************************\n";*/
