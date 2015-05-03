/// \file collision.h
/// \brief Header file for the collision operator
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

#ifndef COLLISION_H
#define COLLISION_H

#include "matrix.h"


//________________________________________________________________________________________________________________________
//

#define GRIDSIZE		128					//!< size of uniform 'k'-grid, must be a power of 2
#define GRIDMASK		(GRIDSIZE - 1)		//!< bit mask for fast modulo 'GRIDSIZE' operation


//________________________________________________________________________________________________________________________
///
/// \brief Point on the collision contour, used for the discretization of the "collision manifold" for the conservative collision operator
///
/// Pre-computed once for the whole simulation
///
typedef struct
{
	double k[4];				//!< 'k' values
	unsigned int ik[4];			//!< indices of 'k' values; interpolation is done once to re-use interpolated values
	double idomega;				//!< mollified 1/(derivative of energy difference with respect to s12)
	unsigned int bindex[2];		//!< neighboring bin indices of k[0]
	double bweight;				//!< relative weight x: x*bin0 + (1-x)*bin1 with x in [0,1]
}
contour_point_t;


//________________________________________________________________________________________________________________________
///
/// \brief Collision contour consisting of an array of points on the contour
///
typedef struct
{
	double *klist;				//!< list of non-uniform 'k' values
	contour_point_t *points;	//!< points on the contour
	unsigned int nk;			//!< length of 'klist'
	unsigned int npts;			//!< length of 'points', i.e., number of points
}
collision_contour_t;


//________________________________________________________________________________________________________________________
///
/// \brief Enumerate paths on collision contours specific for the "next-nearest neighbor" and "exponential decay" dispersion relation
///
typedef enum
{
	GTP_NNN_DIAG,			//!< "diagonal"  contour of next-nearest neighbor dispersion relation
	GTP_NNN_ELLIPSOID,		//!< "ellipsoid" contour of next-nearest neighbor dispersion relation
	GTP_EXP					//!< contour of exponential decay dispersion relation
}
contour_path_t;


//________________________________________________________________________________________________________________________
///
/// \brief Enumerate next-nearest neighbor or exponential decay dispersion relation omega(k)
///
typedef enum
{
	DM_NNN,		//!< next-nearest neighbor dispersion relation
	DM_EXP		//!< exponential decay dispersion relation
}
dispersion_mode_t;


//________________________________________________________________________________________________________________________
//


void ComputeCollisionContour(const double eta, const double epsilon, const contour_path_t path, collision_contour_t *contour);

void DeleteCollisionContour(collision_contour_t *contour);


void CalculateCd(const mat2x2_t *W, const collision_contour_t *contours, const dispersion_mode_t mode, const double eta, const double epsilon, mat2x2_t *dW);


//________________________________________________________________________________________________________________________
//


void ComputeInverseOmegaTable(const dispersion_mode_t mode, const double eta, const double epsilon, double *iomega_table);

void CalculateHeff(const mat2x2_t *W, const double *iomega_table, mat2x2_t *Heff);



#endif
