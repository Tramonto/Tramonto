
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

