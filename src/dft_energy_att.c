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

#include "dft_energy_att.h"

/****************************************************************************/
double integrand_att_freen(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     rho_i = x[iunk][inode_box];

     if (rho_i > 0.) integrand = 0.5*rho_i*int_stencil(x,inode_box,iunk,THETA_PAIRPOT_RCUT);
     return(integrand);
}
/****************************************************************************/
double integrand_att_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;
     int icomp,i;

     i= iunk-Phys2Unk_first[DENSITY];

     if (Grafted_Logical==FALSE || (Type_poly==WJDC3 && Grafted[Icomp_to_polID[i]]==FALSE)){

     if (Lseg_densities){
           icomp = Unk2Comp[i];
           if (Type_interface!=UNIFORM_INTERFACE) rho_bulk = Rho_seg_RTF[i];
           else rho_bulk = rho_bulk = Rho_seg_b[i];
     }
     else{
           icomp = i;
           if (Type_interface!=UNIFORM_INTERFACE) rho_bulk = Rho_b_RTF[i];
           else               rho_bulk = Rho_b[i];
     }

     integrand = 0.5*rho_bulk*Betamu_att[icomp];
     }
     else integrand=0.0;
     return(integrand);
}
/****************************************************************************/
