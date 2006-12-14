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

/* ---------------------------------------------------------
Calculate the thermodynamic properties of an ideal gas fluid.
------------------------------------------------------------*/
#include "dft_thermo_id.h"
/*#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"*/

/******************************************************************************/
/* pressure_ideal_gas: This routine computes the pressure of an ideal gas at the density of interest */
double pressure_ideal_gas(double *rho)
{
   double betap=0.0;
   int i;
   
   for (i=0;i<Ncomp;i++) betap+=rho[i];
   return(betap);
}
/******************************************************************************/
/* chempot_ideal_gas: This routine computes the chemical potential(s) of an ideal gas (mixture) at a density of interest */
void chempot_ideal_gas(double *rho,double *betamu)
{
   int i;
   
   for (i=0;i<Ncomp;i++){
       betamu[i] = log(rho[i]);
                  /*           - 3.0*log(Sigma_ff[icomp][icomp]) -
                               1.5*log(Mass[icomp]*Temp); */
   }
   return;
}
/******************************************************************************/
