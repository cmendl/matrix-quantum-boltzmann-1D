/// \file matrix.c
/// \brief Implementation of some matrix operations
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

#include "matrix.h"
#include "util.h"
#include <math.h>


//________________________________________________________________________________________________________________________
///
/// \brief Frobenius norm of A
///
double FrobeniusNorm(const mat2x2_t a)
{
	// find largest entry
	double dtmax = fabs(a.w0);
	if (fabs(a.wx) > dtmax) { dtmax = fabs(a.wx); };
	if (fabs(a.wy) > dtmax) { dtmax = fabs(a.wy); };
	if (fabs(a.wz) > dtmax) { dtmax = fabs(a.wz); };

	if (dtmax == 0) {
		return 0;
	}

	double r =
		square(a.w0/dtmax) +
		square(a.wx/dtmax) +
		square(a.wy/dtmax) +
		square(a.wz/dtmax);

	return dtmax*sqrt(2*r);
}
