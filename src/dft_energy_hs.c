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

#include "dft_energy_hs.h"

/****************************************************************************/
double integrand_hs_freen(int iunk, int inode_box,double **x)
{
  int i,idim,iv1,iv2;
  double integrand,rho_bar[4+2*NDIM_MAX];
  double n[4+2*NDIM_MAX];

  for(i=0; i<Phys2Nunk[HSRHOBAR]; i++){
      rho_bar[i] = x[i+Phys2Unk_first[HSRHOBAR]][inode_box];
  }
  solutionVec_to_nOrdering(rho_bar,n);
  integrand = phispt_switch(n);
  return(integrand);
}
/*******************************************************************************/
double integrand_hs_freen_bulk(int iunk, int inode_box,double **x)
{
  double integrand,n[4+2*NDIM_MAX],rho_bar[4+2*NDIM_MAX];
  int iv1,iv2,i;

  for(i=0; i<Nrho_bar_s; i++){
     if (Lsteady_state) rho_bar[i] = Rhobar_b_RTF[i];
     else               rho_bar[i] = Rhobar_b[i];
  }
  for (i=0; i<2*Ndim; i++) rho_bar[Nrho_bar_s+i]=0.0;

  solutionVec_to_nOrdering(rho_bar,n);
  integrand = phispt_switch(n);

  return(integrand);
}
/****************************************************************************/

