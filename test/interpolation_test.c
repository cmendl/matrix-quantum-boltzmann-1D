//	Interpolation test file
//
//	Copyright (c) 2012-2015, Christian B. Mendl
//	All rights reserved.
//	http://christian.mendl.net
//
//	This program is free software; you can redistribute it and/or
//	modify it under the terms of the Simplified BSD License
//	http://www.opensource.org/licenses/bsd-license.php
//
//	References:
//	- Martin L.R. F"urst, Christian B. Mendl, Herbert Spohn
//	  Matrix-valued Boltzmann equation for the Hubbard chain
//	  Physical Review E 86, 031122 (2012)
//	  (arXiv:1207.6926)
//
//	- Martin L.R. F"urst, Christian B. Mendl, Herbert Spohn
//	  Matrix-valued Boltzmann equation for the non-integrable Hubbard chain
//	  Physical Review E 88, 012108 (2013)
//	  (arXiv:1302.2075)
//
//	- Jianfeng Lu, Christian B. Mendl
//	  Numerical scheme for a spatially inhomogeneous matrix-valued quantum Boltzmann equation
//	  Journal of Computational Physics 291, 303-316 (2015)
//	  (arXiv:1408.1782)
//________________________________________________________________________________________________________________________
//

#include "interpolation.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


int main()
{
	const unsigned int ngrid = 64;		// number of uniform grid points
	const unsigned int neval = 137;		// number of evaluation points

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// function values on uniform grid
	mat2x2_t *f_grid = (mat2x2_t *)malloc(ngrid*sizeof(mat2x2_t));

	// read values from disk
	if (ReadData("../test/interpolation_test_f_grid.dat", f_grid, sizeof(mat2x2_t), ngrid) < 0)
	{
		fprintf(stderr, "Error reading 'f_grid' data from disk, exiting...\n");
		return -1;
	}

	// generate interpolation
	interpolating_func_t ip_func;
	GenerateInterpolation(f_grid, ngrid, &ip_func);

	double err = 0;

	// simple periodicity consistency check
	{
		mat2x2_t ip0, ip1;

		ip0 = EvaluateInterpolation(&ip_func, 0.0);
		ip1 = EvaluateInterpolation(&ip_func, 1.0);
		err += FrobeniusNorm(SubtractMatrices(ip1, ip0));

		ip0 = EvaluateInterpolation(&ip_func,  0.3);
		ip1 = EvaluateInterpolation(&ip_func, -2.7);
		err += FrobeniusNorm(SubtractMatrices(ip1, ip0));
	}

	// evaluation points
	double *x = (double *)malloc(neval*sizeof(double));
	if (ReadData("../test/interpolation_test_x_eval.dat", x, sizeof(double), neval) < 0)
	{
		fprintf(stderr, "Error reading evaluation points from disk, exiting...\n");
		return -1;
	}

	// evaluate interpolating function
	mat2x2_t *f_eval = (mat2x2_t *)malloc(neval*sizeof(mat2x2_t));
	unsigned int i;
	for (i = 0; i < neval; i++)
	{
		f_eval[i] = EvaluateInterpolation(&ip_func, x[i]);
	}

	// read reference values from disk
	mat2x2_t *f_eref = (mat2x2_t *)malloc(neval*sizeof(mat2x2_t));
	if (ReadData("../test/interpolation_test_f_eref.dat", f_eref, sizeof(mat2x2_t), neval) < 0)
	{
		fprintf(stderr, "Error reading 'f_eref' reference data from disk, exiting...\n");
		return -1;
	}

	// compare with reference
	for (i = 0; i < neval; i++)
	{
		err += FrobeniusNorm(SubtractMatrices(f_eval[i], f_eref[i]));
	}
	printf("cumulative error: %g\n", err);

	// clean up
	DeleteInterpolation(&ip_func);
	free(x);
	free(f_eref);
	free(f_eval);
	free(f_grid);

	return 0;
}
