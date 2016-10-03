#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>

#include <NTL/ZZ.h>

int main(int argc, char **argv){
	int bitLength = atoi(argv[1]);
	int words128=0, words96=0, words64=0, words32=0;
	words128 = bitLength/128;
	int remainingBits = bitLength%128;
	if(remainingBits > 3*32 && remainingBits <= 4*32) { 
		words128++;
	} else if(remainingBits > 2*32 && remainingBits <= 3*32) {
		words96 = 1;
	} else if(remainingBits > 1*32 && remainingBits  <= 2*32) {
		words64 = 1;
	} else if(remainingBits > 0*32 && remainingBits  <= 32) {
		words32 = 1;
	}

	std::cout << words128 << ", " << words96 << ", "  << words64 << ", " << words32 << "\n";
	int roundedBits = words128*128 + words96*96 + words64 * 64 + words32 * 32;
	std::ofstream myfile;
	myfile.open ("outputQMU.cc");
	
	std::string line;
  	std::ifstream constantPart("constants.txt");
  	if (constantPart.is_open()) {
    		while ( std::getline(constantPart,line) ) {
      			myfile << line << "\n";
    		}
	} else std::cout << "Unable to open file"; 
	constantPart.close();
/********************************************************/
	myfile <<  "\nstruct uint" << roundedBits << "{\n";
	int wordsCount = 0;
	for(int i=0;i<words128; i++) {
		myfile << "\tuint128\tw" <<  wordsCount++ << ";\n";
	}
	for(int i=0;i<words96; i++) {
		myfile << "\tuint96\tw" <<  wordsCount++ << ";\n";
	}
	for(int i=0;i<words64; i++) {
		myfile << "\tuint64\tw" <<  wordsCount++ << ";\n";
	}
	for(int i=0;i<words32; i++) {
		myfile << "\tuint32\tw" <<  wordsCount++ << ";\n";
	}
	myfile << "};\n\n\n\t";
/********************************************************/
	myfile << "static inline ZZ to_ZZ(const uint" << roundedBits <<  " & a) { return to_ZZ<uint" << roundedBits  <<">(a); }\n";
	myfile << "static inline ZZ_p to_ZZ_p(const uint"<< roundedBits << " & a) { return to_ZZ_p<uint"<< roundedBits << ">(a); }\n";
	myfile << "static inline uint" << roundedBits << " to_uint" << roundedBits << "(const ZZ & a) { return to_uint<uint" << roundedBits <<">(a); }\n";
	int numberOfWords32 = words128*4 + words96*3 + words64*2 + words32*1;
	int operandsWords = (numberOfWords32-1)/2;
	int madResultWords = operandsWords*2+1;
	int subMultipleResultsWords = operandsWords*2;
	int modulusWords = operandsWords;
	int muWords = operandsWords;
	int qWords = operandsWords+1;
	int rWords = operandsWords+1;
	std::cout<<subMultipleResultsWords<<"\t" << muWords << "\t\n\n";
	bool c=false;
	for(int i=0;i<muWords;i++){
		c = true;
		int l=i;
		for(int j=0;j<subMultipleResultsWords;j++){
			myfile << "\"mad" << (c ? "c" : "") << ".lo.cc.u32\t%" <<   l   << ",%" << subMultipleResultsWords+muWords+i << ",%" << subMultipleResultsWords+muWords+subMultipleResultsWords+j << ",%" << l << ";\\n\\t\"\n\t";
			myfile << "\"mad" << (c ? "c" : "") << ".hi.cc.u32\t%" <<   ++l   << ",%" << subMultipleResultsWords+muWords+i << ",%" << subMultipleResultsWords+muWords+subMultipleResultsWords+j << ",%" << l << ";\\n\\t\"\n\t";
			c = false; 
		}
	}


	myfile.close();
	return 0;

}
