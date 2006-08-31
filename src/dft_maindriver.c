/*
//@HEADER
// ******************************************************************** 
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

/*
 * This is a simple driver routine to call the 
 * main setup and calculation program which 
 * can be built into a library
 * Note that the binary built from this bit of code
 * is a stand alone binary and won't requie the library
 * version
 */
#include <mpi.h>
void dftmain (double *);

int
main (int argc, char *argv[])
{
  int i;
  double dumb;
  MPI_Init (&argc, &argv);
/*  for (i=0; i<1000; i++)*/ 
  dftmain (&dumb);
  MPI_Finalize ();
  return (1);
}

