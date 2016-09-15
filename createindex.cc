#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstring>

#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_p.h>

NTL_CLIENT

int main(int argc, char **argv)
{
    int modulus = 65537;
    int u = argc - 2;
    std::ifstream xcoords(argv++[1], std::ifstream::in);
    std::ifstream * irow = new std::ifstream[u];
    std::ifstream * icol = new std::ifstream[u];

    // determine how many servers, and what evals they will hold
    int num_servers;
    xcoords >> num_servers;
    uint16_t * server_xcoord = new uint16_t[num_servers];
    for (int i = 0; i < num_servers; ++i)
    {
	xcoords >> server_xcoord[i];
    }

    // open the files containing the sparse permutation-like matrices
    for (int i = 1; i <= u; ++i)
    {	
	char * infile = new char[strlen(argv[i]) + 4];
	sprintf(infile, "%s.row", argv[i]);
	irow[i - 1].open(infile, std::ifstream::in);
	sprintf(infile, "%s.col", argv[i]);
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
	sprintf(outfile, "val.%d", server_xcoord[i]);
	oval[i].open(outfile, std::ofstream::out);
	oval[i] << modulus << " ";
    }
    delete [] outfile;

    orow << p << " ";
    ocol << r << " 0 ";
    opoly << modulus << " ";

    // init NTL and precompute Lagrange polynomials
    ZZ_p::init(to_ZZ(modulus));

    ZZ_pX upoly(INIT_SIZE, u);
    vec_ZZ_pX lagrange(INIT_SIZE, u, upoly);
    const ZZ_pX X(INIT_MONO, 1);
    ZZ_pX tmp(INIT_SIZE, u);
    ZZ_p w(INIT_ALLOC);

    vec_ZZ_p interp_xcoords(INIT_SIZE, u);
    ZZ_pX zeros(INIT_SIZE, u + 1);
    for (int i = 0; i < u; ++i) interp_xcoords[i] = to_ZZ_p(i);
    NTL::BuildFromRoots(zeros, interp_xcoords);
    interp_xcoords.kill();

    for (int i = 0; i < u; i++)
    {
	NTL::set(w);
	for (int j = 0; j < u; ++j)
	{
	    if (i == j) continue;
	    w *= (i - j);
	}
	NTL::div(tmp, zeros, (X - i));
	NTL::mul(lagrange[i], tmp, NTL::inv(w));
    }

    // read the input matrices and create the index in CCS format, on the fly
    int * prev_col = new int[u];
    memset(prev_col, 0, sizeof(int)*u);

    vec_ZZ_pX buffer(INIT_SIZE, p, upoly);
    int next_col = 0;
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
		NTL::add(buffer[which_row], buffer[which_row], lagrange[i]);
	    }
	}
	for (int i = 0; i < p; ++i)
	{
	    if (!NTL::IsZero(buffer[i]))
	    {
		next_col++;
		orow << i << " ";
		opoly << buffer[i] << " ";
		for (int k = 0; k < num_servers; ++k)
		{
		    oval[k] << eval(buffer[i], to_ZZ_p(server_xcoord[k])) << " ";
		}
	    }
	    NTL::clear(buffer[i]);
	}
	ocol << next_col << " ";
    }

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
    }
    delete [] irow;
    delete [] icol;
    delete [] oval;
    orow.close();
    ocol.close();
    opoly.close();
    lagrange.kill();
    buffer.kill();

    return 0;
}
