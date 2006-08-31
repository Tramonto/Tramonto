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

/* File for simple global communication operations, e.g. gmax_double */
#include "mpi.h"

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
