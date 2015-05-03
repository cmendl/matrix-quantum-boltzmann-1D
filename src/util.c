/// \file util.c
/// \brief Utility functions
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

#include "util.h"
#include <stdio.h>


//________________________________________________________________________________________________________________________
///
/// \brief Read 'n' items of size 'size' from file 'filename', expecting the file size to be exactly n*size
///
int ReadData(const char *filename, void *data, const size_t size, const size_t n)
{
	FILE *fd = fopen(filename, "rb");
	if (fd == NULL)
	{
		fprintf(stderr, "'fopen()' failed during call of 'ReadData()'.\n");
		return -1;
	}

	// obtain the file size
	fseek(fd, 0 , SEEK_END);
	long filesize = ftell(fd);
	rewind(fd);
	// printf("file size: %d\n", filesize);
	if ((size_t)filesize != n*size)
	{
		fprintf(stderr, "'ReadData()' failed: expected file size does not match.\n");
		return -2;
	}

	// copy the file into the data array
	if (fread(data, size, n, fd) != n)
	{
		fprintf(stderr, "'fread()' failed during call of 'ReadData()'.\n");
		return -3;
	}

	fclose(fd);

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Write 'n' items of size 'size' to file 'filename'
///
int WriteData(const char *filename, const void *data, const size_t size, const size_t n, const bool append)
{
	const char *mode = append ? "ab" : "wb";

	FILE *fd = fopen(filename, mode);
	if (fd == NULL)
	{
		fprintf(stderr, "'fopen()' failed during call of 'WriteData()'.\n");
		return -1;
	}

	// write data array to file
	if (fwrite(data, size, n, fd) != n)
	{
		fprintf(stderr, "'fwrite()' failed during call of 'WriteData()'.\n");
		return -3;
	}

	fclose(fd);

	return 0;
}
