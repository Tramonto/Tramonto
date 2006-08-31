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

/* ---------------------------------------------------------
Calculate the thermodynamic properties of an ideal gas fluid.
------------------------------------------------------------*/
#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

/******************************************************************************/
/* calc_ideal_gas: This routine computes ideal gas properties at a density of interest */
double calc_ideal_gas(double *rho,double *betamu)
{
   double betap=0.0;
   int i;
   
   for (i=0;i<Ncomp;i++){
       betamu[i] = log(rho[i]);
                  /*           - 3.0*log(Sigma_ff[icomp][icomp]) -
                               1.5*log(Mass[icomp]*Temp); */
       betap+=rho[i];
   }
   return(betap);
}
/******************************************************************************/
