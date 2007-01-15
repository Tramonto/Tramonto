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

/* File for simple global communication operations, e.g. gmax_double */
#include "mpi.h"
#include "dft_comm.h"

MPI_Comm Comm=MPI_COMM_WORLD;

/***************************************************************************/
double gsum_double(double c){
  double r;
  (void) MPI_Allreduce(&c, &r, 1,  MPI_DOUBLE, MPI_SUM, Comm);
  return r;
}
/***************************************************************************/
double gmax_double(double c){
  double r;
  (void) MPI_Allreduce(&c, &r, 1,  MPI_DOUBLE, MPI_MAX, Comm);
  return r;
}
/***************************************************************************/
double gmin_double(double c){
  double r;
  (void) MPI_Allreduce(&c, &r, 1,  MPI_DOUBLE, MPI_MIN, Comm);
  return r;
}
/***************************************************************************/
int gsum_int(int c){
  int r;
  (void) MPI_Allreduce(&c, &r, 1,  MPI_INT, MPI_SUM, Comm);
  return r;
}
/***************************************************************************/
int gmax_int(int c){
  int r;
  (void) MPI_Allreduce(&c, &r, 1,  MPI_INT, MPI_MAX, Comm);
  return r;
}
/***************************************************************************/
int gmin_int(int c){
  int r;
  (void) MPI_Allreduce(&c, &r, 1,  MPI_INT, MPI_MIN, Comm);
  return r;
}
/***************************************************************************/
void gsum_vec_int(int* c, int *r, int n){
  (void) MPI_Allreduce(c, r, n,  MPI_INT, MPI_SUM, Comm);
}
/***************************************************************************/
