/*
//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/* --------------------------------------------------------------------------------------
dft_thermo_att.c: Calculate the thermodynamic properties for attractions in the bulk fluid.
---------------------------------------------------------------------------------------*/
#include "dft_thermo_att.h"

/*************************************************************
ATT_thermo_precalc: Call any functions to precalculate useful global
                    bulk parameters. */
void ATT_thermo_precalc()
{
   calc_Avdw_att();
   return;
}
/*************************************************************
calc_Avdw_att: In this routine calculate the strict mean field
              attraction contribution to the pressure */
void calc_Avdw_att()
{
  int icomp,jcomp;

  for (icomp=0; icomp<Ncomp; icomp++)
      for (jcomp=0; jcomp<Ncomp; jcomp++) Avdw[icomp][jcomp]=0.0;

  for (icomp=0; icomp<Ncomp; icomp++) {
     for (jcomp=0; jcomp<Ncomp;jcomp++){
       Avdw[icomp][jcomp] = int_stencil_bulk(THETA_PAIRPOT_RCUT,icomp,jcomp,NULL);
     }
  }
  return;
}
/*************************************************************
pressure_att: In this routine calculate the strict mean field
              attraction contribution to the pressure */
double pressure_att(double *rho)
{
  int icomp,jcomp;
  double betap_att;

  betap_att = 0.0; 

  for (icomp=0; icomp<Ncomp; icomp++) {
     for (jcomp=0; jcomp<Ncomp;jcomp++){
       betap_att += 0.5*Avdw[icomp][jcomp]*rho[icomp]*rho[jcomp];
     }
  }
  return(betap_att);
}
/*************************************************************
chempot_att: In this routine calculate the strict mean field
                     attraction contribution to the pressure and
                     chemical potential */
void chempot_att(double *rho)
{
  int icomp,jcomp;

  for (icomp=0; icomp<Ncomp; icomp++) Betamu_att[icomp] = 0.0;

  for (icomp=0; icomp<Ncomp; icomp++) {
     for (jcomp=0; jcomp<Ncomp;jcomp++){
       Betamu_att[icomp] += rho[jcomp]*Avdw[icomp][jcomp];
     }
  }
  return;
}
/*************************************************************************
dp_drho_att: the derivative of the attractive part of the pressure 
             with respect to rho*/
double dp_drho_att(double *rho)
{
 double dp_drho;

 dp_drho = Avdw[0][0]*rho[0];
 return (dp_drho);
}
/*************************************************************************
dmu_drho_att: the derivative of the attractive part of the chemical potential 
               with respect to rho*/
double dmu_drho_att(double *rho)
{
   double dmu_drho;

   dmu_drho = Avdw[0][0];
   return (dmu_drho);
}
/******************************************************************************/
