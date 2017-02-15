#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <chrono>

#include "gf2e.h"

int main(int argc, char **argv)
{

    auto start = std::chrono::high_resolution_clock::now();
    //std::chrono::nanoseconds onesec{1000000000};

    //for(int d=0;d<argc;d++)cout << argv[d] << endl;
    //int modulus = 65537;
    int u = argc - 3;

    std::ifstream xcoords(argv[2], std::ifstream::in);
    std::ifstream * irow = new std::ifstream[u];
    std::ifstream * icol = new std::ifstream[u];

    // determine how many servers, and what evals they will hold
    int num_servers;
    xcoords >> num_servers;
    GF28_Element * server_xcoord = new GF28_Element[num_servers];
    for (int i = 0; i < num_servers; ++i)
    {
        uint8_t tmp;
	    xcoords >> tmp;
	    server_xcoord[i] = (GF28_Element)tmp;
    }

    // open the files containing the sparse permutation-like matrices
    for (int i = 1; i <= u; ++i)
    {
	char * infile = new char[strlen(argv[i+2]) + 4];
	sprintf(infile, "%s.row", argv[i+2]);
	irow[i - 1].open(infile, std::ifstream::in);
	sprintf(infile, "%s.col", argv[i+2]);
	icol[i - 1].open(infile, std::ifstream::in);
	delete [] infile;
    }

    // determine the dimensions (and ensure that all files agree!)
    int p, r;
    irow[0] >> p;
    icol[0] >> r;
    for (int i = 1; i < u; ++i)
    {
	int _p,_r;
	irow[i] >> _p;
	icol[i] >> _r;
	if (p != _p || r != _r)
	{
		std::cout << "size mismatch!\n";
		exit(1);
	}
    }

    // open the output files
    std::ofstream orow;
    std::ofstream ocol;
    std::ofstream opoly;
    std::ofstream * oval = new std::ofstream[num_servers];
    orow.open("out.row", std::ofstream::out);
    ocol.open("out.col", std::ofstream::out);
    opoly.open("out.polys", std::ofstream::out);
    char * outfile = new char[9];
    for (int i = 0; i < num_servers; ++i)
    {
	sprintf(outfile, "val.%u", server_xcoord[i]);
	oval[i].open(outfile, std::ofstream::out);
	oval[i] << "256 ";
    }
    delete [] outfile;

    orow << p << " ";
    long opos = orow.tellp();
    orow << "           ";
    ocol << r << " 0 ";
    opoly << "256 ";

    GF28_Element * interp_xcoords = new GF28_Element[u];
    for (int i = 0; i < u; ++i) interp_xcoords[i] = (GF28_Element)i;
    GF28_Element * precomp = new GF28_Element[num_servers * u];
    GF28_Element * invals = new GF28_Element[u];
    memset(invals, 0, u * sizeof(GF28_Element));

    for (int i = 0; i < num_servers; i++)
    {
        for (int j = 0; j < u; j++)
        {
            invals[j] = 1;
            if (j>0) invals[j-1] = 0;
            precomp[i*num_servers+j] = interpolate_GF2E<GF28_Element>(interp_xcoords, invals, u, server_xcoord[i]);
        }
        invals[u-1] = 0;
    }

    // read the input matrices and create the index in CCS format, on the fly
    int * prev_col = new int[u];
    memset(prev_col, 0, sizeof(int)*u);

    int num_vals = 0;
    GF28_Element * buffer = new GF28_Element[num_servers * p];

    int next_col = 0;
    bool * nz = new bool[p];
    for (int i = 0; i < p; i++) nz[i] = false;
    for (int j = 0; j < r; ++j)
    {
	for (int i = 0; i < u; ++i)
	{
	    int to_read = prev_col[i];
	    icol[i] >> prev_col[i];
	    to_read = prev_col[i] - to_read;
	    while (to_read--)
	    {
		int which_row;
		irow[i] >> which_row;
		//NTL::add(polybuf[which_row], polybuf[which_row], lagrange[i]);
		for (int k = 0; k < num_servers; k++)
		{
			buffer[k*num_servers+which_row] ^= precomp[k*num_servers+i];
		}
		nz[which_row] = true;
	    }
	}
	for (int i = 0; i < p; ++i)
	{
	    if (nz[i])
	    {
		next_col++;
		orow << i << " ";
		//opoly << polybuf[i] << " ";
		for (int k = 0; k < num_servers; ++k)
		{
		    oval[k] << buffer[k*num_servers+i] << " ";
		    buffer[k*num_servers+i] = 0;
		}
                num_vals++;
	    }
	    nz[i] = false;
	}
	ocol << next_col << " ";
    }
    orow.seekp(opos);
    orow << num_vals;

    // cleanup
    delete [] prev_col;
    for (int i = 0; i < u; ++i)
    {
	irow[i].close();
	icol[i].close();
    }
    for (int i = 0; i < num_servers; ++i)
    {
	oval[i].close();
    }\
    delete [] irow;
    delete [] icol;
    delete [] oval;
    orow.close();
    ocol.close();
    opoly.close();
    delete [] nz;

//auto time_elapsed = ;
    std::cout << "Time requires to batch CCS files: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << " milliseconds\n";

    return 0;
}
