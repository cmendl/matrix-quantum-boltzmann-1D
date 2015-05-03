/// \file util.h
/// \brief Header file for the utility functions
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

#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include <stddef.h>


//________________________________________________________________________________________________________________________
///
/// \brief Square function x -> x^2
///
static inline double square(const double x)
{
	return x*x;
}


//________________________________________________________________________________________________________________________
//


int ReadData(const char *filename, void *data, const size_t size, const size_t n);

int WriteData(const char *filename, const void *data, const size_t size, const size_t n, const bool append);



#endif
