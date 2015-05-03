/// \file boltzmann_ode.c
/// \brief Main file for solving the time evolution of the Boltzmann equation
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

#include "boltzmann_ode.h"
#include "collision.h"
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Solve Boltzmann ODE (with respect to time) using the midpoint rule;
/// using the dissipative collision operator Cd only
///
/// \param W0 initial Wigner state, array of size 'GRIDSIZE'
/// \param mode dispersion mode
/// \param dt time step
/// \param eta parameter of the dispersion relation omega(k)
/// \param epsilon mollification parameter
/// \param numsteps number of time steps
/// \param Wevolv output: must point to an array of size 'GRIDSIZE x numsteps'
///
void CalculateWevolvDissipative(const mat2x2_t *W0, const dispersion_mode_t mode, const double dt, const double eta, const double epsilon, const unsigned int numsteps, mat2x2_t *Wevolv)
{
	collision_contour_t contours[2];
	if (mode == DM_NNN)
	{
		// collision contours for "diagonal" and "ellipsoid" integration paths
		ComputeCollisionContour(eta, epsilon, GTP_NNN_DIAG,      &contours[0]);
		ComputeCollisionContour(eta, epsilon, GTP_NNN_ELLIPSOID, &contours[1]);
	}
	else
	{
		assert(mode == DM_EXP);
		// single collision contour for "diagonal" integration path
		ComputeCollisionContour(eta, epsilon, GTP_EXP, &contours[0]);
	}

	// first entry of 'Wevolv' is W0
	memcpy(Wevolv, W0, GRIDSIZE*sizeof(mat2x2_t));

	// temporary storage (for midpoint rule)
	mat2x2_t *Wtmp = (mat2x2_t *)malloc(GRIDSIZE*sizeof(mat2x2_t));

	mat2x2_t *Cd = (mat2x2_t *)malloc(GRIDSIZE*sizeof(mat2x2_t));

	unsigned int nt;
	for (nt = 1; nt < numsteps; nt++)
	{
		unsigned int i;

		// midpoint rule

		// Wtmp = Wevolv[nt-1] + dt/2 * Cd
		CalculateCd(&Wevolv[(nt-1)*GRIDSIZE], contours, mode, eta, epsilon, Cd);
		for (i = 0; i < GRIDSIZE; i++)
		{
			Wtmp[i] = ScalarMultiplyAddMatrix(dt/2, Cd[i], Wevolv[(nt-1)*GRIDSIZE+i]);
		}

		// Wevolv[nt] = Wevolv[nt-1] + dt * Cd
		CalculateCd(Wtmp, contours, mode, eta, epsilon, Cd);
		for (i = 0; i < GRIDSIZE; i++)
		{
			Wevolv[nt*GRIDSIZE+i] = ScalarMultiplyAddMatrix(dt, Cd[i], Wevolv[(nt-1)*GRIDSIZE+i]);
		}
	}

	// clean up
	if (mode == DM_NNN) {
		DeleteCollisionContour(&contours[1]);
	}
	DeleteCollisionContour(&contours[0]);
	free(Wtmp);
	free(Cd);
}


//________________________________________________________________________________________________________________________
///
/// \brief Solve Boltzmann ODE (with respect to time) using the midpoint rule;
/// using the conservative collision operator Cc only
///
/// \param W0 initial Wigner state, array of size 'GRIDSIZE'
/// \param dt time step
/// \param iomega_table precomputed mollified inverse omega table
/// \param numsteps number of time steps
/// \param Wevolv output: must point to an array of size 'GRIDSIZE x numsteps'
///
void CalculateWevolvConservative(const mat2x2_t *W0, const double dt, const double *iomega_table, const unsigned int numsteps, mat2x2_t *Wevolv)
{
	// first entry of 'Wevolv' is W0
	memcpy(Wevolv, W0, GRIDSIZE*sizeof(mat2x2_t));

	// temporary storage (for midpoint rule)
	mat2x2_t *Wtmp = (mat2x2_t *)malloc(GRIDSIZE*sizeof(mat2x2_t));

	mat2x2_t *Heff = (mat2x2_t *)malloc(GRIDSIZE*sizeof(mat2x2_t));

	unsigned int nt;
	for (nt = 1; nt < numsteps; nt++)
	{
		unsigned int i;

		// midpoint rule

		// Wtmp = Wevolv[nt-1] + dt/2 * Cc
		CalculateHeff(&Wevolv[(nt-1)*GRIDSIZE], iomega_table, Heff);
		for (i = 0; i < GRIDSIZE; i++)
		{
			const mat2x2_t C = Commutator(Wevolv[(nt-1)*GRIDSIZE+i], Heff[i]);	// -i[Heff, W1]
			Wtmp[i] = ScalarMultiplyAddMatrix(dt/2, C, Wevolv[(nt-1)*GRIDSIZE+i]);
		}

		// Wevolv[nt] = Wevolv[nt-1] + dt * Cc
		CalculateHeff(Wtmp, iomega_table, Heff);
		for (i = 0; i < GRIDSIZE; i++)
		{
			const mat2x2_t C = Commutator(Wtmp[i], Heff[i]);		// -i[Heff, W1]
			Wevolv[nt*GRIDSIZE+i] = ScalarMultiplyAddMatrix(dt, C, Wevolv[(nt-1)*GRIDSIZE+i]);
		}
	}

	// clean up
	free(Wtmp);
	free(Heff);
}


//________________________________________________________________________________________________________________________
///
/// \brief Solve Boltzmann ODE (with respect to time) using the midpoint rule;
/// using the full (dissipative and conservative) collision operator
///
/// \param W0 initial Wigner state, array of size 'GRIDSIZE'
/// \param mode dispersion mode
/// \param dt time step
/// \param eta parameter of the dispersion relation omega(k)
/// \param epsilonCd mollification parameter for the dissipative collision operator
/// \param iomega_table precomputed mollified inverse omega table
/// \param numsteps number of time steps
/// \param Wevolv output: must point to an array of size 'GRIDSIZE x numsteps'
///
void CalculateWevolv(const mat2x2_t *W0, const dispersion_mode_t mode, const double dt, const double eta, const double epsilonCd, const double *iomega_table, const unsigned int numsteps, mat2x2_t *Wevolv)
{
	collision_contour_t contours[2];
	if (mode == DM_NNN)
	{
		// collision contours for "diagonal" and "ellipsoid" integration paths
		ComputeCollisionContour(eta, epsilonCd, GTP_NNN_DIAG,      &contours[0]);
		ComputeCollisionContour(eta, epsilonCd, GTP_NNN_ELLIPSOID, &contours[1]);
	}
	else
	{
		assert(mode == DM_EXP);
		// single collision contour for "diagonal" integration path
		ComputeCollisionContour(eta, epsilonCd, GTP_EXP, &contours[0]);
	}

	// first entry of 'Wevolv' is W0
	memcpy(Wevolv, W0, GRIDSIZE*sizeof(mat2x2_t));

	// temporary storage (for midpoint rule)
	mat2x2_t *Wtmp = (mat2x2_t *)malloc(GRIDSIZE*sizeof(mat2x2_t));

	mat2x2_t *Cd   = (mat2x2_t *)malloc(GRIDSIZE*sizeof(mat2x2_t));
	mat2x2_t *Heff = (mat2x2_t *)malloc(GRIDSIZE*sizeof(mat2x2_t));

	unsigned int nt;
	for (nt = 1; nt < numsteps; nt++)
	{
		unsigned int i;

		// midpoint rule

		// Wtmp = Wevolv[nt-1] + dt/2 * C
		CalculateCd  (&Wevolv[(nt-1)*GRIDSIZE], contours, mode, eta, epsilonCd, Cd);
		CalculateHeff(&Wevolv[(nt-1)*GRIDSIZE], iomega_table, Heff);
		for (i = 0; i < GRIDSIZE; i++)
		{
			const mat2x2_t C = AddMatrices(Cd[i], Commutator(Wevolv[(nt-1)*GRIDSIZE+i], Heff[i]));	// Cd - i[Heff, W1]
			Wtmp[i] = ScalarMultiplyAddMatrix(dt/2, C, Wevolv[(nt-1)*GRIDSIZE+i]);
		}

		// Wevolv[nt] = Wevolv[nt-1] + dt * C
		CalculateCd  (Wtmp, contours, mode, eta, epsilonCd, Cd);
		CalculateHeff(Wtmp, iomega_table, Heff);
		for (i = 0; i < GRIDSIZE; i++)
		{
			const mat2x2_t C = AddMatrices(Cd[i], Commutator(Wtmp[i], Heff[i]));		// Cd - i[Heff, W1]
			Wevolv[nt*GRIDSIZE+i] = ScalarMultiplyAddMatrix(dt, C, Wevolv[(nt-1)*GRIDSIZE+i]);
		}
	}

	// clean up
	if (mode == DM_NNN) {
		DeleteCollisionContour(&contours[1]);
	}
	DeleteCollisionContour(&contours[0]);
	free(Wtmp);
	free(Heff);
	free(Cd);
}
