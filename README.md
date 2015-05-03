Matrix-valued Boltzmann equation for the nonintegrable Hubbard chain
====================================================================

C source code, MathLink interface and Mathematica demonstration file for the simulations in "Matrix-valued Boltzmann equation for the nonintegrable Hubbard chain" (see References).

Compiling the source code:
- Windows: Visual Studio project files are provided in the *vcproj* folder (standalone test files) and the *vcproj_mlink* folder (Mathematica MathLink interface)
- Linux and similar: see the makefile in the *test* folder (standalone test files) and the *mlink* folder (Mathematica MathLink interface;
  you may have to adapt paths according to your local installation, and manually create the output directory, e.g., *mlink/Linux-x86-64*)

The Mathematica `boltzmann_demo.cdf` demonstration file can be viewed with
the free [CDF Player](http://www.wolfram.com/cdf-player) or opened and edited with [Mathematica](http://www.wolfram.com/mathematica).


License
-------
Copyright (c) 2012-2015, Christian B. Mendl  
All rights reserved.  
http://christian.mendl.net

This program is free software; you can redistribute it and/or
modify it under the terms of the Simplified BSD License
http://www.opensource.org/licenses/bsd-license.php


References
----------
1. Martin L.R. Fürst, Christian B. Mendl, Herbert Spohn  
   Matrix-valued Boltzmann equation for the Hubbard chain  
   Physical Review E 86, 031122 (2012), [arXiv:1207.6926](http://arxiv.org/abs/1207.6926)
2. Martin L.R. Fürst, Christian B. Mendl, Herbert Spohn  
   Matrix-valued Boltzmann equation for the nonintegrable Hubbard chain  
   Physical Review E 88, 012108 (2013), [arXiv:1302.2075](http://arxiv.org/abs/1302.2075)
3. Jianfeng Lu, Christian B. Mendl  
   Numerical scheme for a spatially inhomogeneous matrix-valued quantum Boltzmann equation  
   Journal of Computational Physics 291, 303-316 (2015), [arXiv:1408.1782](http://arxiv.org/abs/1408.1782)
