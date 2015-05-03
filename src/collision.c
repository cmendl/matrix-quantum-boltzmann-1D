/// \file collision.c
/// \brief Parametrize the "collision manifold" and calculate the collision operator
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

#include "collision.h"
#include "interpolation.h"
#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


// 2*pi
#define M_2PI		6.283185307179586476925


//________________________________________________________________________________________________________________________
///
/// \brief Dispersion relation of the next-nearest neighbor model
///
static inline double OmegaNNN(const double k, const double eta)
{
	return 1 - cos(M_2PI*k) - eta*cos(2*M_2PI*k);
}

//________________________________________________________________________________________________________________________
///
/// \brief Dispersion relation of the exponential decay model
///
static inline double OmegaExp(const double k, const double eta)
{
	return -0.5*(1 + sinh(eta)/(cosh(eta) - cos(M_2PI*k)));
}


//________________________________________________________________________________________________________________________
///
/// \brief Energy conservation term contribution along "trivial" contour
///
static inline double Etriv(const double d12, const double d34)
{
	return 4*sin(M_PI*(d12 + d34))*sin(M_PI*(d12 - d34));
}

//________________________________________________________________________________________________________________________
///
/// \brief Energy conservation term contribution along "trivial" contour, alternative formulation
///
static inline double EtrivAlt(const double k1, const double k3, const double k4)
{
	return 4*sin(M_PI*(k1 - k3))*sin(M_PI*(k1 - k4));
}


//________________________________________________________________________________________________________________________
///
/// \brief Energy conservation term contribution along non-trivial contour (next-nearest neighbor model)
///
static inline double EcontourNNN(const double s12, const double d12, const double d34, const double eta)
{
	return cos(M_PI*s12) + 2*eta*cos(M_2PI*s12)*(cos(M_2PI*d12) + cos(M_2PI*d34));
}

//________________________________________________________________________________________________________________________
///
/// \brief Energy conservation term contribution along non-trivial contour (exponential decay model)
///
static inline double EcontourExp(const double s12, const double d12, const double d34, const double eta)
{
	const double coshe  = cosh(eta);
	const double cos12  = cos(M_2PI*d12);
	const double cos34  = cos(M_2PI*d34);
	const double cossum = cos(M_PI*s12);

	return cossum*(1 + square(coshe) + cos12*cos34 - square(cossum)) - coshe*(cos12 + cos34);
}


//________________________________________________________________________________________________________________________
///
/// \brief Energy conservation term contribution along non-trivial contour (next-nearest neighbor model), alternative formulation
///
static inline double EcontourNNNAlt(const double k1, const double k3, const double k4, const double eta)
{
	return EcontourNNN(k3+k4, k1-0.5*(k3+k4), 0.5*(k3-k4), eta);
}

//________________________________________________________________________________________________________________________
///
/// \brief Energy conservation term contribution along non-trivial contour (exponential decay model), alternative formulation
///
static inline double EcontourExpAlt(const double k1, const double k3, const double k4, const double eta)
{
	return EcontourExp(k3+k4, k1-0.5*(k3+k4), 0.5*(k3-k4), eta);
}


//________________________________________________________________________________________________________________________
///
/// \brief Energy difference omega_1 + omega_2 - omega_3 - omega_4 (next-nearest neighbor model)
///
static inline double EnergyDifferenceNNN(const double k1, const double k3, const double k4, const double eta)
{
	return EtrivAlt(k1, k3, k4)*EcontourNNNAlt(k1, k3, k4, eta);
}

//________________________________________________________________________________________________________________________
///
/// \brief Energy difference omega_1 + omega_2 - omega_3 - omega_4 (exponential decay model)
///
static inline double EnergyDifferenceExp(const double k1, const double k3, const double k4, const double eta)
{
	const double coshe = cosh(eta);
	const double denom =
		(coshe - cos(M_2PI*k1))*
		(coshe - cos(M_2PI*(k3+k4-k1)))*
		(coshe - cos(M_2PI*k3))*
		(coshe - cos(M_2PI*k4));

	return 0.5*sinh(eta)/denom * EtrivAlt(k1, k3, k4) * EcontourExpAlt(k1, k3, k4, eta);
}


//________________________________________________________________________________________________________________________
///
/// \brief Derivative of the energy difference with respect to s12 (next-nearest neighbor), except "trivial" factor
/// (equal to derivative of 'EcontourNNN')
///
static inline double dNontrivEnergyDifferenceNNN(const double s12, const double d12, const double d34, const double eta)
{
	return -M_PI*sin(M_PI*s12) - 4*M_PI*eta*sin(M_2PI*s12)*(cos(M_2PI*d12) + cos(M_2PI*d34));
}

//________________________________________________________________________________________________________________________
///
/// \brief Derivative of the energy difference with respect to s12 (exponential decay), except "trivial" factor
///
static inline double dNontrivEnergyDifferenceExp(const double s12, const double d12, const double d34, const double eta)
{
	const double coshe = cosh(eta);

	// denominator
	double denom, ddenom;
	{
		int i;

		double k2pi[4];
		k2pi[0] = M_PI*(s12 + 2*d12);
		k2pi[1] = M_PI*(s12 - 2*d12);
		k2pi[2] = M_PI*(s12 + 2*d34);
		k2pi[3] = M_PI*(s12 - 2*d34);

		denom  = 1;
		ddenom = 0;
		for (i = 0; i < 4; i++)
		{
			double cosk = cos(k2pi[i]);
			double sink = sin(k2pi[i]);

			denom  *= (coshe - cosk);
			ddenom += sink/(coshe - cosk);
		}
		ddenom *= -M_PI;
	}

	// "contour" factor
	double contour, dcontour;
	{
		const double cos12  = cos(M_2PI*d12);
		const double cos34  = cos(M_2PI*d34);
		const double cossum = cos(M_PI*s12);

		contour  = cossum*(1 + square(coshe) + cos12*cos34 - square(cossum)) - coshe*(cos12 + cos34);
		dcontour = -M_PI*sin(M_PI*s12)*(1 + square(coshe) + cos12*cos34 - 3*square(cossum));
	}

	return 0.5*sinh(eta)*(dcontour + contour*ddenom)/denom;
}


//________________________________________________________________________________________________________________________
///
/// \brief Derivative d(EnergyDifference)/dk4 along "trivial" solution path gamma_2 (k2 == k3) (next-nearest neighbor model)
///
static inline double dEnergyDifferenceNNNGamma2(const double k1, const double k2, const double eta)
{
	return M_2PI*(sin(M_2PI*k2) - sin(M_2PI*k1) + 2*eta*(sin(4*M_PI*k2)-sin(4*M_PI*k1)));
}

//________________________________________________________________________________________________________________________
///
/// \brief Derivative d(EnergyDifference)/dk4 along "trivial" solution path gamma_2 (k2 == k3) (exponential decay model)
///
static inline double dEnergyDifferenceExpGamma2(const double k1, const double k2, const double eta)
{
	const double coshe = cosh(eta);

	return M_PI*sinh(eta)*(sin(M_2PI*k2)/square(coshe - cos(M_2PI*k2)) - sin(M_2PI*k1)/square(coshe - cos(M_2PI*k1)));
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the sum s12 on the "diagonal" contour, see Eq. (28),
/// for the next-nearest neighbor model
/// with dispersion omega(k) = 1 - cos(2 pi k) - eta*cos(4 pi k)
///
static double CalculateSum12NNNDiag(const double d12, const double d34, const double eta)
{
	// independent of signs of d12 and d34, and d12 <-> d34
	const double t = eta*(cos(M_2PI*d12) + cos(M_2PI*d34));

	if (fabs(t) < 0.01)
	{
		// Taylor expansion around t = 0
		const double tsq = t*t;
		return 0.5 + t*(-2.0 + 44.0*tsq/3.0 - 1132.0*tsq*tsq/5.0 + 31096.0*tsq*tsq*tsq/7.0)/M_PI;
	}
	else
	{
		return acos((sqrt(1 + 32*t*t) - 1)/(8*t))/M_PI;
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Check whether parameters d12, d34, eta are valid (next-nearest neighbor model)
///
static bool VerifySum12NNNEllipsoidParams(const double d12, const double d34, const double eta)
{
	// independent of signs of d12 and d34, and d12 <-> d34
	const double t = eta*(cos(M_2PI*d12) + cos(M_2PI*d34));

	return fabs(t) >= 0.5;
}

//________________________________________________________________________________________________________________________
///
/// \brief Calculate the sum s12 on the ellipsoid contour (next-nearest neighbor model)
///
static double CalculateSum12NNNEllipsoid(const double d12, const double d34, const double eta)
{
	// independent of signs of d12 and d34, and d12 <-> d34
	const double t = eta*(cos(M_2PI*d12) + cos(M_2PI*d34));
	assert(fabs(t) >= 0.5);

	return acos(-(sqrt(1 + 32*t*t) + 1)/(8*t))/M_PI;
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the sum s12 on the non-trivial contour for the exponential decay model
/// with dispersion omega(k) = - sum_{j=0}^infty exp(-alpha j) cos(2 pi j k)
///
static double CalculateSum12Exp(const double d12, const double d34, const double eta)
{
	const double coshe = cosh(eta);
	const double cos12 = cos(M_2PI*d12);
	const double cos34 = cos(M_2PI*d34);

	// independent of signs of d12 and d34, and d12 <-> d34
	const double b = 1 + square(coshe) + cos12*cos34;
	const double c = coshe*(cos12 + cos34);

	// Newtons's method for solving x^3 - b*x + c == 0
	double x = 0;
	double f = c;
	do
	{
		x = x - f/(3*x*x - b);
		f = x*x*x - b*x + c;
	}
	while (fabs(f) > 4*DBL_EPSILON);

	return acos(x)/M_PI;
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the integrand of the dissipative collision operator: A + A^*
///
mat2x2_t CalculateAAT(const mat2x2_t W1, const mat2x2_t W2, const mat2x2_t W3, const mat2x2_t W4)
{
	// w3 . eta . w4
	const double w3_eta_w4 = W3.w0*W4.w0 - (W3.wx*W4.wx + W3.wy*W4.wy + W3.wz*W4.wz);

	// gain term
	mat2x2_t tmp = W2;
	tmp.w0 = 1 - tmp.w0;
	mat2x2_t gain = ScalarMultiplyMatrix(2 * w3_eta_w4, tmp);

	// loss term
	const double t = (W3.w0 + W4.w0) - 1;
	mat2x2_t loss;
	loss.w0 = w3_eta_w4 - t*W2.w0;
	loss.wx = t*W2.wx;
	loss.wy = t*W2.wy;
	loss.wz = t*W2.wz;

	return SubtractMatrices(gain, AntiCommutator(loss, W1));
}


//________________________________________________________________________________________________________________________
///
/// \brief Smooth inverse (for mollification, dissipative collision operator)
///
static inline double SmoothInvCd(const double x, const double epsilon)
{
	return 1/sqrt(square(x) + square(epsilon));
}

//________________________________________________________________________________________________________________________
///
/// \brief Smooth inverse (for mollification, conservative collision operator)
///
static inline double SmoothInvCc(const double x, const double epsilon)
{
	return x/(square(x) + square(epsilon));
}


//________________________________________________________________________________________________________________________
///
/// \brief Auxiliary structure used for sorting
///
typedef struct
{
	double k;		//!< k value
	int index;		//!< index
}
k_index_t;

//________________________________________________________________________________________________________________________
///
/// \brief Comparison function used for sorting
///
static int CompareKIndex(const void *s1, const void *s2)
{
	const k_index_t *a = (k_index_t *)s1;
	const k_index_t *b = (k_index_t *)s2;

	if (a->k < b->k)
	{
		return -1;
	}
	else if (a->k > b->k)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Comparison function used for sorting
///
static int CompareContourPoint(const void *s1, const void *s2)
{
	const contour_point_t *a = (contour_point_t *)s1;
	const contour_point_t *b = (contour_point_t *)s2;

	if (a->k[0] < b->k[0])
	{
		return -1;
	}
	else if (a->k[0] > b->k[0])
	{
		return 1;
	}
	else
	{
		if (a->k[2] < b->k[2])
		{
			return -1;
		}
		else if (a->k[2] > b->k[2])
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Sort 'klist' and update indices in 'points'
///
static void SortContourKlist(collision_contour_t *contour)
{
	unsigned int i;

	k_index_t *ki = (k_index_t *)malloc(contour->nk*sizeof(k_index_t));

	// interleave indices with k list
	for (i = 0; i < contour->nk; i++)
	{
		ki[i].k = contour->klist[i];
		ki[i].index = i;
	}

	qsort(ki, contour->nk, sizeof(k_index_t), CompareKIndex);

	// copy 'k' values back
	for (i = 0; i < contour->nk; i++)
	{
		contour->klist[i] = ki[i].k;
	}

	// construct reverse index map
	unsigned int *revind = (unsigned int *)malloc(contour->nk*sizeof(unsigned int));
	for (i = 0; i < contour->nk; i++)
	{
		revind[ki[i].index] = i;
	}

	// update indices in 'points'
	for (i = 0; i < contour->npts; i++)
	{
		unsigned int j;
		for (j = 0; j < 4; j++)
		{
			contour->points[i].ik[j] = revind[contour->points[i].ik[j]];
		}
	}

	free(revind);
	free(ki);
}


//________________________________________________________________________________________________________________________
//

/// \brief Function type of the dispersion relation omega
typedef double (*OmegaPtr)(const double k, const double eta);

/// \brief Function type for calculating the sum s12 
typedef double (*CalculateSum12Ptr)(const double d12, const double d34, const double eta);

/// \brief Function type for calculating the derivative of the energy difference with respect to s12, except "trivial" factor
typedef double (*dNontrivEnergyDifferencePtr)(const double s12, const double d12, const double d34, const double eta);

//________________________________________________________________________________________________________________________
///
/// \brief Compute collision contour
///
void ComputeCollisionContour(const double eta, const double epsilon, const contour_path_t path, collision_contour_t *contour)
{
	const unsigned int scale = (path == GTP_NNN_ELLIPSOID) ? 2 : 1;

	// construct contour

	static const int permtable[8][4] = { { 0, 1, 2, 3 }, { 1, 0, 2, 3 }, { 0, 1, 3, 2 }, { 1, 0, 3, 2 }, { 2, 3, 0, 1 }, { 3, 2, 0, 1 }, { 2, 3, 1, 0 }, { 3, 2, 1, 0 } };
 	static const int signtable[8][2] = { { 1,    1 },    {-1,    1    }, { 1,   -1    }, {-1,   -1    }, { 1,    1 },    {-1,    1    }, { 1,   -1    }, {-1,   -1    } };

	// assign function pointers
	CalculateSum12Ptr CalculateSum12;
	OmegaPtr Omega;
	dNontrivEnergyDifferencePtr dNontrivEnergyDifference;
	switch (path)
	{
		case GTP_NNN_DIAG:
			CalculateSum12 = CalculateSum12NNNDiag;
			Omega = OmegaNNN;
			dNontrivEnergyDifference = dNontrivEnergyDifferenceNNN;
			break;
		case GTP_NNN_ELLIPSOID:
			CalculateSum12 = CalculateSum12NNNEllipsoid;
			Omega = OmegaNNN;
			dNontrivEnergyDifference = dNontrivEnergyDifferenceNNN;
			break;
		case GTP_EXP:
			CalculateSum12 = CalculateSum12Exp;
			Omega = OmegaExp;
			dNontrivEnergyDifference = dNontrivEnergyDifferenceExp;
			break;
		default:
			assert(false);
			return;
	}

	// evaluate omega function at uniform grid
	double *omegaTable = (double *)malloc(GRIDSIZE*sizeof(double));
	unsigned int i;
	for (i = 0; i < GRIDSIZE; i++)
	{
		omegaTable[i] = Omega((double)i/GRIDSIZE, eta);
	}

	// allocate sufficient memory (some will remain unused for 'ellipsoid')
	contour->klist  =          (double *)malloc(GRIDSIZE*(GRIDSIZE/2+1)/scale*sizeof(double));
	contour->points = (contour_point_t *)malloc(GRIDSIZE*GRIDSIZE/scale*sizeof(contour_point_t));
	// number of points is initially zero
	contour->nk   = 0;
	contour->npts = 0;
	// current point
	contour_point_t *gp = contour->points;

	for (i = 0; i <= GRIDSIZE/(2*scale); i++)		// TODO: upper limit GRIDSIZE/4 for 'ellipsoid' could depend on eta
	{
		const double d12 = (double)i/GRIDSIZE;

		unsigned int j;
		for (j = i; j <= GRIDSIZE/(2*scale); j++)	// TODO: upper limit GRIDSIZE/4 for 'ellipsoid' could depend on eta
		{
			const double d34 = (double)j/GRIDSIZE;

			if (path == GTP_NNN_ELLIPSOID && !VerifySum12NNNEllipsoidParams(d12, d34, eta))
			{
				continue;
			}

			// independent of signs of d12 and d34, and d12 <-> d34
			const double s12 = CalculateSum12(d12, d34, eta);

			unsigned int n;
			for (n = 0; n < scale; n++)			// for s12 and 2-s12 (only for 'ellipsoid')
			{
				const double s12h = 0.5*(n == 0 ? s12 : 2-s12);
				double k[4];
				k[0] = s12h + d12;
				k[1] = s12h - d12;
				k[2] = s12h + d34;
				k[3] = s12h - d34;
				// ensure numerically exact same value (analytically the same due to mod 1)
				if (i == GRIDSIZE/2) { k[1] = k[0]; }
				if (j == GRIDSIZE/2) { k[3] = k[2]; }

				// k index
				int ikL[4];
				ikL[0] = contour->nk;
				if (0 < i && i < GRIDSIZE/2) { contour->nk++; }
				ikL[1] = contour->nk;
				contour->nk++;
				if (i < j)
				{
					ikL[2] = contour->nk;
					if (0 < j && j < GRIDSIZE/2) { contour->nk++; }
					ikL[3] = contour->nk;
					contour->nk++;
				}
				else
				{
					ikL[2] = ikL[0];
					ikL[3] = ikL[1];
				}
				assert(contour->nk <= GRIDSIZE*(GRIDSIZE/2+1)/scale);

				unsigned int m;
				for (m = 0; m < 4; m++)
				{
					// modulo 1
					double intpart;
					k[m] = modf(k[m] + 1, &intpart);
					// store k value
					contour->klist[ikL[m]] = k[m];
				}

				// derivative of the energy difference with respect to s12, except trivial factor;
				// independent of signs of d12 and d34, and d12 <-> d34
				const double dec = dNontrivEnergyDifference(n == 0 ? s12 : 2-s12, d12, d34, eta);

				// for each permutation id, (d12 <-> -d12), (d34 <-> -d34), (d12 <-> -d12, d34 <-> -d34), as well as d12 <-> d34
				unsigned int l;
				for (l = 0; l < (i < j ? 8u : 4u); l++)
				{
					// skip d12 <-> -d12 or d34 <-> -d34 for d12 = 0, 1/2 or d34 = 0, 1/2, respectively
					if ((i == 0 || i == GRIDSIZE/2) && signtable[l][l < 4 ? 0 : 1] == -1) { continue; }
					if ((j == 0 || j == GRIDSIZE/2) && signtable[l][l < 4 ? 1 : 0] == -1) { continue; }

					for (m = 0; m < 4; m++)
					{
						// contour point entry
						gp-> k[m] =   k[permtable[l][m]];
						// corresponding indices in 'klist'
						gp->ik[m] = ikL[permtable[l][m]];
					}

					// omega evaluated at current k
					const double omega = Omega(gp->k[0], eta);

					// bin indices
					gp->bindex[0] = (int)(GRIDSIZE*gp->k[0]);			// casting to int takes the floor value
					gp->bindex[1] = (gp->bindex[0] + 1) & GRIDMASK;		// next index
					assert(0 <= gp->bindex[0] && gp->bindex[0] < GRIDSIZE && (double)gp->bindex[0]/GRIDSIZE <= gp->k[0]);

					// affine weight factor (can be > 1 or < 0 at extrema of omega)
					gp->bweight = (omegaTable[gp->bindex[1]] - omega)/(omegaTable[gp->bindex[1]] - omegaTable[gp->bindex[0]]);
					// gp->bweight == 1 for eta == 0

					// inverse of derivative dOmega
					// indepdendent of d12 <-> d34 due to square in 'SmoothInvCd'
					gp->idomega = M_PI*SmoothInvCd(dec*Etriv(signtable[l][0]*d12, signtable[l][1]*d34), epsilon);

					// next contour point
					gp++;
					contour->npts++;
					assert(contour->npts <= GRIDSIZE*GRIDSIZE/scale);
				}
			}
		}
	}
	assert(path == GTP_NNN_ELLIPSOID || contour->nk == GRIDSIZE*(GRIDSIZE/2+1)/scale);

	// sort contour points with respect to k1
	qsort(contour->points, contour->npts, sizeof(contour_point_t), CompareContourPoint);

	// sort 'klist' and update indices in 'points'
	SortContourKlist(contour);

	// consistency check
	for (i = 0; i < contour->npts; i++)
	{
		unsigned int j;
		for (j = 0; j < 4; j++)
		{
			assert(contour->klist[contour->points[i].ik[j]] == contour->points[i].k[j]);
		}
	}

	free(omegaTable);
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete collision contour (free memory)
///
void DeleteCollisionContour(collision_contour_t *contour)
{
	if (contour->klist)
	{
		free(contour->klist);
		contour->klist = NULL;
	}

	if (contour->points)
	{
		free(contour->points);
		contour->points = NULL;
	}

	contour->npts = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Function type for calculating derivative d(EnergyDifference)/dk4 along "trivial" solution path gamma_2
///
typedef double (*dEnergyDifferenceGamma2Ptr)(double k1, double k2, double eta);


//________________________________________________________________________________________________________________________
///
/// \brief Calculate dissipative collision operator Cd along "trivial" integration path
///
/// \param W Wigner matrices for each 'k', array of length 'GRIDSIZE'
/// \param mode dispersion mode
/// \param eta dispersion relation parameter
/// \param epsilon mollification parameter
/// \param dW output, must point to an array of length 'GRIDSIZE' on input
///
static void CalculateCdTriv(const mat2x2_t *W, const dispersion_mode_t mode, const double eta, const double epsilon, mat2x2_t *dW)
{
	dEnergyDifferenceGamma2Ptr dEnergyDifferenceGamma2 = (mode == DM_NNN ? dEnergyDifferenceNNNGamma2 : dEnergyDifferenceExpGamma2);

	unsigned int i;
	for (i = 0; i < GRIDSIZE; i++)
	{
		const double k1 = (double)i/GRIDSIZE;

		unsigned int j;
		for (j = 0; j < GRIDSIZE; j++)
		{
			const double k3 = (double)j/GRIDSIZE;

			// \gamma_2 integration path (k1 = k4, k2 = k3);
			// contribution by \gamma_1 path is equal to \gamma_2 for the convention in 'CalculateAAT()', hence the factor 2;
			// multiply by mollified delta function term
			mat2x2_t dWcur = ScalarMultiplyMatrix(2 * M_PI/GRIDSIZE * SmoothInvCd(dEnergyDifferenceGamma2(k1, k3, eta), epsilon),
				CalculateAAT(W[i], W[j], W[j], W[i]));

			// add to dW[i]
			dW[i] = AddMatrices(dW[i], dWcur);
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate dissipative collision operator Cd on collision contour
///
/// \param Wfunc interpolating function for Wigner function
/// \param contour collision contour
/// \param dW output, must point to an array of length 'GRIDSIZE' on input
///
static void CalculateCdContour(const interpolating_func_t *Wfunc, const collision_contour_t *contour, mat2x2_t *dW)
{
	unsigned int i;

	// interpolate at 'klist'
	mat2x2_t *Wip = (mat2x2_t *)malloc(contour->nk*sizeof(mat2x2_t));
	for (i = 0; i < contour->nk; i++)
	{
		Wip[i] = EvaluateInterpolation(Wfunc, contour->klist[i]);
	}

	for (i = 0; i < contour->npts; i++)
	{
		const contour_point_t *gp = &contour->points[i];

		// multiply by mollified delta function term
		mat2x2_t dWcur = ScalarMultiplyMatrix(gp->idomega/GRIDSIZE,
			CalculateAAT(Wip[gp->ik[0]], Wip[gp->ik[1]], Wip[gp->ik[2]], Wip[gp->ik[3]]));

		// add affine combination to dW
		dW[gp->bindex[0]] = AddMatrices(dW[gp->bindex[0]], ScalarMultiplyMatrix(    gp->bweight, dWcur));
		dW[gp->bindex[1]] = AddMatrices(dW[gp->bindex[1]], ScalarMultiplyMatrix(1 - gp->bweight, dWcur));
	}

	// clean up
	free(Wip);
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate dissipative collision operator Cd
///
/// \param W Wigner matrices for each 'k', array of length 'GRIDSIZE'
/// \param contours contours of the collision manifold
/// \param mode dispersion mode
/// \param eta dispersion relation parameter
/// \param epsilon mollification parameter
/// \param dW output, must point to an array of length 'GRIDSIZE' on input
///
void CalculateCd(const mat2x2_t *W, const collision_contour_t *contours, const dispersion_mode_t mode, const double eta, const double epsilon, mat2x2_t *dW)
{
	// generate interpolating function
	interpolating_func_t Wfunc;
	GenerateInterpolation(W, GRIDSIZE, &Wfunc);

	// initially set to zero
	memset(dW, 0, GRIDSIZE*sizeof(mat2x2_t));

	// "trivial" integration path
	CalculateCdTriv(W, mode, eta, epsilon, dW);

	if (mode == DM_NNN)
	{
		// "diagonal" and "ellipsoid" integration paths
		CalculateCdContour(&Wfunc, &contours[0], dW);
		CalculateCdContour(&Wfunc, &contours[1], dW);
	}
	else
	{
		assert(mode == DM_EXP);

		// "diagonal" integration path
		CalculateCdContour(&Wfunc, &contours[0], dW);
	}

	// clean up
	DeleteInterpolation(&Wfunc);
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate integrand of effective Hamiltonian required for conservative collision operator
///
static inline mat2x2_t HeffIntegrand(const mat2x2_t W2, const mat2x2_t W3, const mat2x2_t W4)
{
	return ScalarMultiplyMatrix(1 - (W3.w0 + W4.w0), W2);
}


//________________________________________________________________________________________________________________________
///
/// \brief Energy difference function type
///
typedef double (*EnergyDifferencePtr)(const double k1, const double k3, const double k4, const double eta);


//________________________________________________________________________________________________________________________
///
/// \brief Pre-compute mollified inverse omega table
///
/// \param mode dispersion mode (next-nearest neighbor or exponential decay)
/// \param eta dispersion relation parameter
/// \param epsilon mollification parameter
/// \param iomega_table output: inverse omega table; must point to an array of length 'GRIDSIZE^3' at input 
///
void ComputeInverseOmegaTable(const dispersion_mode_t mode, const double eta, const double epsilon, double *iomega_table)
{
	EnergyDifferencePtr EnergyDifference = (mode == DM_NNN ? EnergyDifferenceNNN : EnergyDifferenceExp);

	unsigned int i1;
	for (i1 = 0; i1 < GRIDSIZE; i1++)
	{
		const double k1 = (double)i1/GRIDSIZE;

		unsigned int i3;
		for (i3 = 0; i3 < GRIDSIZE; i3++)
		{
			const double k3 = (double)i3/GRIDSIZE;

			unsigned int i4;
			for (i4 = 0; i4 < GRIDSIZE; i4++)
			{
				const double k4 = (double)i4/GRIDSIZE;

				*iomega_table = SmoothInvCc(EnergyDifference(k1, k3, k4, eta), epsilon)/(GRIDSIZE*GRIDSIZE);
				iomega_table++;
			}
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate effective Hamiltonian
///
/// \param W Wigner matrices for each 'k', array of length 'GRIDSIZE'
/// \param iomega_table precomputed mollified inverse omega table
/// \param Heff output, must point to an array of length 'GRIDSIZE' on input
///
void CalculateHeff(const mat2x2_t *W, const double *iomega_table, mat2x2_t *Heff)
{
	memset(Heff, 0, GRIDSIZE*sizeof(mat2x2_t));

	unsigned int i1;
	for (i1 = 0; i1 < GRIDSIZE; i1++)
	{
		unsigned int i3;
		for (i3 = 0; i3 < GRIDSIZE; i3++)
		{
			unsigned int i4;
			for (i4 = 0; i4 < GRIDSIZE; i4++)
			{
				// calculate integrand of effective Hamiltonian
				unsigned int i2 = (i3 + i4 - i1) & GRIDMASK;
				mat2x2_t dHcur = HeffIntegrand(W[i2], W[i3], W[i4]);

				// multiply by mollified principal value term
				dHcur = ScalarMultiplyMatrix(*iomega_table, dHcur);
				iomega_table++;

				// add to dW[i1]
				Heff[i1] = AddMatrices(Heff[i1], dHcur);
			}
		}
	}
}
