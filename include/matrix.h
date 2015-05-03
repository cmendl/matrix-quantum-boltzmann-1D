/// \file matrix.h
/// \brief Header file defining a Hermitian 2x2 matrix structure, represented in the basis of the identity and the Pauli matrices
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

#ifndef MATRIX_H
#define MATRIX_H


//________________________________________________________________________________________________________________________
///
/// \brief Hermitian 2x2 matrix, represented in the basis of the identity and the Pauli matrices:
/// w_0 id + w_x sigma_x + w_y sigma_y + w_z sigma_z
///
typedef struct
{
	double w0;		//!< coefficient of the identity matrix
	double wx;		//!< coefficient of the Pauli sigma_x matrix, equal to half the Bloch vector x-component
	double wy;		//!< coefficient of the Pauli sigma_y matrix, equal to half the Bloch vector y-component
	double wz;		//!< coefficient of the Pauli sigma_z matrix, equal to half the Bloch vector z-component
}
mat2x2_t;


//________________________________________________________________________________________________________________________
///
/// \brief A + B
///
static inline mat2x2_t AddMatrices(const mat2x2_t a, const mat2x2_t b)
{
	mat2x2_t ret;
	ret.w0 = a.w0 + b.w0;
	ret.wx = a.wx + b.wx;
	ret.wy = a.wy + b.wy;
	ret.wz = a.wz + b.wz;

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief A - B
///
static inline mat2x2_t SubtractMatrices(const mat2x2_t a, const mat2x2_t b)
{
	mat2x2_t ret;
	ret.w0 = a.w0 - b.w0;
	ret.wx = a.wx - b.wx;
	ret.wy = a.wy - b.wy;
	ret.wz = a.wz - b.wz;

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief x*A
///
static inline mat2x2_t ScalarMultiplyMatrix(const double x, const mat2x2_t a)
{
	mat2x2_t ret;
	ret.w0 = x * a.w0;
	ret.wx = x * a.wx;
	ret.wy = x * a.wy;
	ret.wz = x * a.wz;

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief x*A + B
///
static inline mat2x2_t ScalarMultiplyAddMatrix(const double x, const mat2x2_t a, const mat2x2_t b)
{
	mat2x2_t ret;
	ret.w0 = x * a.w0 + b.w0;
	ret.wx = x * a.wx + b.wx;
	ret.wy = x * a.wy + b.wy;
	ret.wz = x * a.wz + b.wz;

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief tr(A)
///
static inline double MatrixTrace(const mat2x2_t a)
{
	return 2 * a.w0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Matrix commutator i [A,B] == i*(A*B - B*A)
///
static inline mat2x2_t Commutator(const mat2x2_t a, const mat2x2_t b)
{
	mat2x2_t ret;
	ret.w0 = 0;
	ret.wx = 2 * (b.wy*a.wz - b.wz*a.wy);
	ret.wy = 2 * (b.wz*a.wx - b.wx*a.wz);
	ret.wz = 2 * (b.wx*a.wy - b.wy*a.wx);

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief Matrix anti-commutator {A,B} == A*B + B*A
///
static inline mat2x2_t AntiCommutator(const mat2x2_t a, const mat2x2_t b)
{
	mat2x2_t ret;
	ret.w0 = 2 * (a.w0*b.w0 + a.wx*b.wx + a.wy*b.wy + a.wz*b.wz);
	ret.wx = 2 * (a.w0*b.wx + a.wx*b.w0);
	ret.wy = 2 * (a.w0*b.wy + a.wy*b.w0);
	ret.wz = 2 * (a.w0*b.wz + a.wz*b.w0);

	return ret;
}


//________________________________________________________________________________________________________________________
//


double FrobeniusNorm(const mat2x2_t a);



#endif
