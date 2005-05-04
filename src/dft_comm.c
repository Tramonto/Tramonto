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
