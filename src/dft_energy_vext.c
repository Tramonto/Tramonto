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

/*#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"*/
#include "dft_energy_vext.h"

/****************************************************************************/
double integrand_vext_freen(int iunk,int inode_box, double **x)
{     double integrand,rho_i;
     int icomp,iseg,i;

     i = iunk-Phys2Unk_first[DENSITY];
     if (Type_poly==WTC) icomp = Unk2Comp[i];
     else                icomp = i;

     rho_i = x[iunk][inode_box];

     if (rho_i > DENSITY_MIN) integrand = rho_i*Vext[B2L_node[inode_box]][icomp];
     else integrand=0.0;
     
     return(integrand);
}
/****************************************************************************/
double integrand_vext_elec_freen(int iunk,int inode_box, double **x)
{     double integrand,rho_i;
     int icomp,iseg,i;

     i = iunk-Phys2Unk_first[DENSITY];
     if (Type_poly==WTC) icomp = Unk2Comp[i];
     else                icomp = i;

     rho_i = x[iunk][inode_box];

     if (rho_i > DENSITY_MIN) integrand = 0.5*rho_i*Vext_coul[B2L_node[inode_box]][icomp];
     else integrand=0.0;
     
     return(integrand);
}
