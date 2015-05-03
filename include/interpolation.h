/// \file interpolation.h
/// \brief Header file for interpolation functionality
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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "matrix.h"


//________________________________________________________________________________________________________________________
///
/// \brief Matrix-valued divided differences for polynomial interpolation of order 3, periodic in the interval [0, 1];
/// interpolating polynomial in Horner form
///
typedef struct
{
	mat2x2_t dA[4];		//!< matrix-valued divided differences
}
divided_diff_t;


//________________________________________________________________________________________________________________________
///
/// \brief Matrix-valued interpolating function of order 3, periodic in the interval [0, 1]
///
typedef struct
{
	divided_diff_t *dd;		//!< array of divided differences
	unsigned int num;		//!< length of 'dd'
}
interpolating_func_t;


void GenerateInterpolation(const mat2x2_t *f, const unsigned int num, interpolating_func_t *ip);

void DeleteInterpolation(interpolating_func_t *ip);


mat2x2_t EvaluateInterpolation(const interpolating_func_t *ip, const double x);



#endif
