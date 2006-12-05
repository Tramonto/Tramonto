/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
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

