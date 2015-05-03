/// \file boltzmannML.c
/// \brief MathLink file for the Mathematica interface
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
#include "util.h"
#include <mathlink.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#define strcasecmp _strcmpi
#endif


//________________________________________________________________________________________________________________________
//


void SimulateBoltzmann(double dt, double eta, double epsilonCd, double epsilonCc, int numsteps, const char *modestr, double *W0list, long W0len, int useCc)
{
	// check input arguments
	if (dt <= 0 || epsilonCd <= 0 || epsilonCc <= 0 || numsteps <= 0)
	{
		MLEvaluate(stdlink, "Message[SimulateBoltzmann::InvalidArg]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	dispersion_mode_t mode;
	if (strcasecmp(modestr, "nnn") == 0)
	{
		mode = DM_NNN;
	}
	else if (strcasecmp(modestr, "exp") == 0)
	{
		mode = DM_EXP;
	}
	else
	{
		MLEvaluate(stdlink, "Message[SimulateBoltzmann::InvalidMode]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	if (W0len != GRIDSIZE*sizeof(mat2x2_t)/sizeof(double))
	{
		MLEvaluate(stdlink, "Message[SimulateBoltzmann::InvalidDimensionW0]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}
	mat2x2_t *W0 = (mat2x2_t *)W0list;

	mat2x2_t *Wevolv = (mat2x2_t *)malloc(GRIDSIZE*numsteps*sizeof(mat2x2_t));
	if (Wevolv == NULL)
	{
		MLEvaluate(stdlink, "Message[SimulateBoltzmann::OutOfMemory]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	if (useCc == 1)
	{
		// pre-compute mollified inverse omega table
		double *iomega_table = (double *)malloc(GRIDSIZE*GRIDSIZE*GRIDSIZE*sizeof(double));
		if (iomega_table == NULL)
		{
			MLEvaluate(stdlink, "Message[SimulateBoltzmann::OutOfMemory]");
			// discard 'ReturnPacket'
			MLNextPacket(stdlink);
			MLNewPacket(stdlink);	// discard
			// final output
			MLPutSymbol(stdlink, "$Failed");
			return;
		}
		ComputeInverseOmegaTable(mode, eta, epsilonCc, iomega_table);

		// perform main calculation
		CalculateWevolv(W0, mode, dt, eta, epsilonCd, iomega_table, numsteps, Wevolv);

		free(iomega_table);
	}
	else if (useCc == 0)
	{
		// only dissipative collision operator

		// perform main calculation
		CalculateWevolvDissipative(W0, mode, dt, eta, epsilonCd, numsteps, Wevolv);
	}
	else if (useCc == -1)
	{
		// only conservative collision operator

		// pre-compute mollified inverse omega table
		double *iomega_table = (double *)malloc(GRIDSIZE*GRIDSIZE*GRIDSIZE*sizeof(double));
		if (iomega_table == NULL)
		{
			MLEvaluate(stdlink, "Message[SimulateBoltzmann::OutOfMemory]");
			// discard 'ReturnPacket'
			MLNextPacket(stdlink);
			MLNewPacket(stdlink);	// discard
			// final output
			MLPutSymbol(stdlink, "$Failed");
			return;
		}
		ComputeInverseOmegaTable(mode, eta, epsilonCc, iomega_table);

		// perform main calculation
		CalculateWevolvConservative(W0, dt, iomega_table, numsteps, Wevolv);

		free(iomega_table);
	}
	else
	{
		MLEvaluate(stdlink, "Message[SimulateBoltzmann::InvalidArg]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	// return results to Mathematica (matrices in Pauli representation)
	int dims[3] = { numsteps, GRIDSIZE, 4 };
	MLPutReal64Array(stdlink, (const double *)Wevolv, dims, NULL, 3);

	MLEndPacket(stdlink);
	MLFlush(stdlink);

	// clean up
	free(Wevolv);
}


//________________________________________________________________________________________________________________________
//


#if MACINTOSH_MATHLINK

int main(int argc, char* argv[])
{
	/* Due to a bug in some standard C libraries that have shipped with
	 * MPW, zero is passed to MLMain below.  (If you build this program
	 * as an MPW tool, you can change the zero to argc.)
	 */
	argc = argc; /* suppress warning */
	return MLMain(0, argv);
}

#elif WINDOWS_MATHLINK

#if __BORLANDC__
#pragma argsused
#endif

int PASCAL WinMain(HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char  buff[512];
	char FAR * buff_start = buff;
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;

	hinstPrevious = hinstPrevious; /* suppress warning */

	if (!MLInitializeIcon(hinstCurrent, nCmdShow)) return 1;
	MLScanString(argv, &argv_end, &lpszCmdLine, &buff_start);
	return MLMain((int)(argv_end - argv), argv);
}

#else

int main(int argc, char* argv[])
{
	return MLMain(argc, argv);
}

#endif
