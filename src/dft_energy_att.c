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

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
double integrand_att_freen(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     rho_i = x[iunk][inode_box];

     int_stencil(x,inode_box,iunk,U_ATTRACT);
     if (rho_i > 0.) integrand = 0.5*rho_i*Temporary_sum;
     return(integrand);
}
/****************************************************************************/
double integrand_att_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;
     int icomp,i;

     i= iunk-Phys2Unk_first[DENSITY];
     if (Type_poly==WTC){
           icomp = Unk2Comp[i];
           if (Lsteady_state) rho_bulk = Rho_seg_RTF[i];
           else rho_bulk = rho_bulk = Rho_seg_b[i];
     }
     else{
           icomp = i;
           if (Lsteady_state) rho_bulk = Rho_b_RTF[i];
           else               rho_bulk = Rho_b[i];
     }

     integrand = 0.5*rho_bulk*Betamu_att[icomp];
     return(integrand);
}
/****************************************************************************/
