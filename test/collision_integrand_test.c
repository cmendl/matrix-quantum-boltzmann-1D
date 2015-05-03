//	Collision operator integrand test file
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


mat2x2_t CalculateAAT(const mat2x2_t W1, const mat2x2_t W2, const mat2x2_t W3, const mat2x2_t W4);



int main()
{
	const mat2x2_t w1 = {  1,  2,  3,  4 };
	const mat2x2_t w2 = {  2,  3, -5,  7 };
	const mat2x2_t w3 = { -1,  2, 13, 11 };
	const mat2x2_t w4 = {  5,  1, 17,  2 };

	printf("W1:\n");
	PrintMatrix(w1);

	printf("\nW2:\n");
	PrintMatrix(w2);

	printf("\nW3:\n");
	PrintMatrix(w3);

	printf("\nW4:\n");
	PrintMatrix(w4);

	double err = 0;

	// A + A^*
	{
		// reference
		const mat2x2_t AATref = { 898, -494, 4066, -1494 };

		printf("\nA + A^*:\n");
		mat2x2_t w = CalculateAAT(w1, w2, w3, w4);
		PrintMatrix(w);
		err += FrobeniusNorm(SubtractMatrices(w, AATref));
	}

	printf("\ncumulative error: %g\n", err);

	return 0;
}
