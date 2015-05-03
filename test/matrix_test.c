//	Matrix operations test file
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

#include "matrix.h"
#include "collision.h"
#include <math.h>
#include <stdio.h>


static void PrintMatrix(const mat2x2_t a)
{
	printf("| %g     %g%+gi|\n", a.w0 + a.wz, a.wx, -a.wy);
	printf("| %g%+gi   %g |\n",  a.wx,  a.wy, a.w0 - a.wz);
}


//________________________________________________________________________________________________________________________
//


int main()
{
	const mat2x2_t a = {  1,  2,  3,  4 };
	const mat2x2_t b = {  2,  3, -5,  7 };

	printf("A:\n");
	PrintMatrix(a);

	printf("\nB:\n");
	PrintMatrix(b);

	double err = 0;

	// commutator
	{
		// reference for commutator
		const mat2x2_t comm_ab = { 0, -82, 4, 38 };

		printf("\ni*[A,B]:\n");
		mat2x2_t w = Commutator(a, b);
		PrintMatrix(w);
		err += FrobeniusNorm(SubtractMatrices(w, comm_ab));
	}

	// anti-commutator
	{
		// reference for anti-commutator
		const mat2x2_t acom_ab = { 42, 14, 2, 30 };

		printf("\n{A,B}:\n");
		mat2x2_t w = AntiCommutator(a, b);
		PrintMatrix(w);
		err += FrobeniusNorm(SubtractMatrices(w, acom_ab));
	}

	// Frobenius norm
	{
		// reference for Frobenius norm
		const double nfro = 2*sqrt(15.0);

		double w = FrobeniusNorm(a);
		printf("\nFrobenius norm of A: %g\n", w);
		err += fabs(w - nfro);
	}

	printf("\ncumulative error: %g\n", err);

	return 0;
}
