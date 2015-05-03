/// \file interpolation.c
/// \brief Matrix-valued interpolation based on divided differences
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
//	  Matrix-valued Boltzmann equation for the nonintegrable Hubbard chain
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
#include <math.h>
#include <stdlib.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Generate interpolating function of order 3 by calculating divided differences
///
/// \param f matrix-valued function evaluated at grid points; array of length 'num'
/// \param num number of grid points
/// \param ip interpolating function
///
void GenerateInterpolation(const mat2x2_t *f, const unsigned int num, interpolating_func_t *ip)
{
	ip->num = num;

	ip->dd = (divided_diff_t *)malloc(num*sizeof(divided_diff_t));

	// first step: copy 'f', with index shifted by 1 
	unsigned int i;
	for (i = 1; i < num; i++)
	{
		ip->dd[i].dA[0] = f[i-1];
	}
	ip->dd[0].dA[0] = f[num-1];

	unsigned int j;
	for (j = 1; j < 4; j++)
	{
		const double s = (double)num/j;

		for (i = 0; i < num-1; i++)
		{
			ip->dd[i].dA[j] = ScalarMultiplyMatrix(s, SubtractMatrices(ip->dd[i+1].dA[j-1], ip->dd[i].dA[j-1]));
		}
		ip->dd[num-1].dA[j] = ScalarMultiplyMatrix(s, SubtractMatrices(ip->dd[0].dA[j-1], ip->dd[num-1].dA[j-1]));
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete interpolation (free memory)
///
void DeleteInterpolation(interpolating_func_t *ip)
{
	free(ip->dd);
	ip->dd  = NULL;
	ip->num = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Interpolate divided differences at point 'x'
///
static inline mat2x2_t InterpolateDividedDifferences(const divided_diff_t *diff, const double h, const double x)
{
	assert(0 <= x && x < h);

	return
		ScalarMultiplyAddMatrix(x+h,
		ScalarMultiplyAddMatrix(x,
		ScalarMultiplyAddMatrix(x-h,
			diff->dA[3],
			diff->dA[2]),
			diff->dA[1]),
			diff->dA[0]);
}


//________________________________________________________________________________________________________________________
///
/// \brief Evaluate interpolation at point 'x', periodic in the interval [0, 1]
///
mat2x2_t EvaluateInterpolation(const interpolating_func_t *ip, const double x)
{
	// modulo 1
	double intpart;
	double y = modf(x, &intpart);
	if (y < 0) {
		y++;
	}
	assert(0 <= y && y < 1);

	// grid spacing
	const double h = 1.0/ip->num;

	unsigned int i = (unsigned int)(ip->num*y);		// casting to int takes floor of value
	assert(0 <= i && i < ip->num);
	return InterpolateDividedDifferences(&ip->dd[i], h, y - i*h);
}
