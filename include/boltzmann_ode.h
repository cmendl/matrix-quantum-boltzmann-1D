/// \file boltzmann_ode.h
/// \brief Header file for solving the time evolution of the Boltzmann equation
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

#ifndef BOLTZMANN_ODE_H
#define BOLTZMANN_ODE_H

#include "collision.h"


void CalculateWevolvDissipative(const mat2x2_t *W0, const dispersion_mode_t mode, const double dt, const double eta, const double epsilon, const unsigned int numsteps, mat2x2_t *Wevolv);


void CalculateWevolvConservative(const mat2x2_t *W0, const double dt, const double *iomega_table, const unsigned int numsteps, mat2x2_t *Wevolv);


void CalculateWevolv(const mat2x2_t *W0, const dispersion_mode_t mode, const double dt, const double eta, const double epsilonCd, const double *iomega_table, const unsigned int numsteps, mat2x2_t *Wevolv);



#endif
