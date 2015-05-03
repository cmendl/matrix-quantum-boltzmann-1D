
:Begin:
:Function:		SimulateBoltzmann
:Pattern:		SimulateBoltzmann[dt_Real, eta_Real, epsilonCd_Real, epsilonCc_Real, numsteps_Integer, modestr_String, W0_List, useCc_Integer]
:Arguments:		{ dt,   eta,  epsilonCd, epsilonCc, numsteps, modestr, W0,       useCc }
:ArgumentTypes:	{ Real, Real, Real,      Real,      Integer,  String,  RealList, Integer }
:ReturnType:	Manual
:End:

:Evaluate: SimulateBoltzmann::usage = "SimulateBoltzmann[dt_Real, eta_Real, epsilonCd_Real, epsilonCc_Real, numsteps_Integer, modestr_String, W0_List, useCc_Integer] simulates a matrix-valued quantum Boltzmann equation based on the Hubbard chain for fermions."
:Evaluate: SimulateBoltzmann::InvalidArg = "invalid argument; syntax: SimulateBoltzmann[dt_Real, eta_Real, epsilonCd_Real, epsilonCc_Real, numsteps_Integer, modestr_String, W0_List, useCc_Integer]"
:Evaluate: SimulateBoltzmann::InvalidDimensionW0 = "W0 must be a RealList of length 4*GRIDSIZE"
:Evaluate: SimulateBoltzmann::InvalidMode = "modestr must be either 'nnn' or 'exp'"
:Evaluate: SimulateBoltzmann::OutOfMemory = "malloc failed, probably out of memory"
