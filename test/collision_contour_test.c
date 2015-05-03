//	Collision contour (manifold) test file
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
#include "util.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


static void PrintContourPoint(const contour_point_t *gp)
{
	printf("k:  %g, %g, %g, %g\n", gp-> k[0], gp-> k[1], gp-> k[2], gp-> k[3]);
	printf("ik: %i, %i, %i, %i\n", gp->ik[0], gp->ik[1], gp->ik[2], gp->ik[3]);
	printf("mollified 1/domega: %g\n", gp->idomega);
	printf("bin indices: (%d,%d)\n", gp->bindex[0], gp->bindex[1]);
}


//________________________________________________________________________________________________________________________
//


static inline double maxf(const double x, const double y)
{
	if (x >= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}


static inline double minf(const double x, const double y)
{
	if (x <= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}


//________________________________________________________________________________________________________________________
//


int main()
{
	const double eta = 1.0/3;
	const double epsilon = 0.5;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	assert(GRIDSIZE == 128);
	printf("size of uniform 'k' grid: %i\n", GRIDSIZE);

	printf("eta (dispersion relation parameter): %g\n", eta);
	printf("epsilon (mollification parameter): %g\n", epsilon);

	const contour_path_t paths[] = { GTP_NNN_DIAG, GTP_NNN_ELLIPSOID, GTP_EXP };
	const char *pathnames[]      = {    "nnn_diag",   "nnn_ellipsoid",   "exp" };

	double err = 0;

	unsigned int p;
	for (p = 0; p < 3; p++)		// for each path...
	{
		printf("\n_____________________________________________________________________________\n");
		printf("\nintegration path: '%s'\n", pathnames[p]);

		// compute collision contour
		collision_contour_t contour;
		ComputeCollisionContour(eta, epsilon, paths[p], &contour);

		if (contour.npts == 0)
		{
			printf("contour empty!\n");
			continue;
		}

		printf("number of points:   '%i'\n", contour.npts);
		printf("number of k-values: '%i'\n", contour.nk);

		printf("\npoints[0]:\n");
		PrintContourPoint(&contour.points[0]);

		printf("\npoints[1]:\n");
		PrintContourPoint(&contour.points[1]);

		printf("\npoints[contour.npts-1]:\n");
		PrintContourPoint(&contour.points[contour.npts-1]);

		printf("\nfirst few entries of klist: %g, %g, %g, ..., %g\n",
			contour.klist[0], contour.klist[1], contour.klist[2], contour.klist[contour.nk-1]);

		// compare 'klist' with reference
		{
			// read reference values of 'klist' from disk
			char filepath[1024];
			sprintf(filepath, "../test/collision_contour_ref/%s_klist_ref.dat", pathnames[p]);
			double *klist_ref = (double *)malloc(contour.nk * sizeof(double));
			if (ReadData(filepath, klist_ref, sizeof(double), contour.nk) < 0) {
				return -1;
			}

			// compare
			unsigned int i;
			for (i = 0; i < contour.nk; i++)
			{
				err = maxf(err, fabs(contour.klist[i] - klist_ref[i]));
			}

			// clean up
			free(klist_ref);
		}

		// compare 'k'-values of contour points with reference
		{
			// read reference 'k'-values of contour points from disk
			char filepath[1024];
			sprintf(filepath, "../test/collision_contour_ref/%s_contour_k4_ref.dat", pathnames[p]);
			double *pts_k4_ref = (double *)malloc(4*contour.npts * sizeof(double));
			if (ReadData(filepath, pts_k4_ref, sizeof(double), 4*contour.npts) < 0) {
				return -1;
			}

			// compare
			unsigned int i;
			for (i = 0; i < contour.npts; i++)
			{
				unsigned int j;
				for (j = 0; j < 4; j++)
				{
					err = maxf(err, fabs(contour.points[i].k[j] - pts_k4_ref[4*i+j]));
				}
			}

			// clean up
			free(pts_k4_ref);
		}

		// compare mollified inverse of derivative of energy difference with reference
		{
			// read reference values from disk
			char filepath[1024];
			sprintf(filepath, "../test/collision_contour_ref/%s_idomega_ref.dat", pathnames[p]);
			double *idomega_ref = (double *)malloc(contour.npts * sizeof(double));
			if (ReadData(filepath, idomega_ref, sizeof(double), contour.npts) < 0) {
				return -1;
			}

			// compare
			unsigned int i;
			for (i = 0; i < contour.npts; i++)
			{
				err = maxf(err, fabs(contour.points[i].idomega - idomega_ref[i]));
			}

			// clean up
			free(idomega_ref);
		}

		// clean up
		DeleteCollisionContour(&contour);
	}

	printf("\nmaximum error: %g\n", err);

	return 0;
}
