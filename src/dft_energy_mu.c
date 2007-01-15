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
#include "dft_energy_mu.h"

/****************************************************************************/
double integrand_mu_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i,mu_i;
     int icomp,iseg,i,unk_mu;

     i = iunk-Phys2Unk_first[DENSITY];

     if (Lsteady_state){
        unk_mu = Phys2Unk_first[DIFFUSION]+i;
        mu_i = x[unk_mu][inode_box];
     }
     else {
        if (Lseg_densities){
               icomp = Unk2Comp[i];
               mu_i = Betamu_seg[i];
               if (Type_coul != NONE) mu_i -= log(Rho_seg_b[i]);
        }
        else{
            mu_i = Betamu[i];
            if (Type_coul != NONE && fabs(Charge_f[i])>1.e-15) mu_i -= log(Rho_b[i]);
        }
     }
     rho_i = x[iunk][inode_box];

/*     printf("CHEMPOT IN SOLN: iunk=%d inode_box=%d  rho_i=%9.6f mu_i=%9.6f\n",iunk,inode_box,rho_i,mu_i);*/
     if (rho_i > DENSITY_MIN) integrand = -rho_i*mu_i;
     else integrand=0.0;

     return(integrand);
}

/****************************************************************************/
double integrand_mu_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_i,mu_i;
     int icomp,iseg,i;

     i = iunk-Phys2Unk_first[DENSITY];

     if (Lsteady_state){
        if (Lseg_densities){
               icomp = Unk2Comp[i];
               /*mu_i = Betamu_seg_RTF[i];*/ /* note we need to fix up WTC for diffusion */
               mu_i = Betamu_RTF[i];
               if (Type_coul != NONE) mu_i -= log(Rho_seg_RTF[i]);
               rho_i = Rho_seg_RTF[i];
        }
        else{   mu_i = Betamu_RTF[i];
               if (Type_coul != NONE) mu_i -= log(Rho_b_RTF[i]);
               rho_i = Rho_b_RTF[i];
        }
     }
     else {
        if (Type_poly==WTC){
               icomp = Unk2Comp[i];
               mu_i = Betamu_seg[i];
               if (Type_coul != NONE) mu_i -= log(Rho_seg_b[i]);
               rho_i = Rho_seg_b[i];
        }
        else{   mu_i = Betamu[i];
                if (Type_coul != NONE && fabs(Charge_f[i])>1.e-15) mu_i -= log(Rho_b[i]);
                rho_i = Rho_b[i];
        }
     }

/*     printf("BULK TERM: iunk=%d inode_box=%d rho_i=%9.6f mu_i=%9.6f\n",iunk,inode_box,rho_i,mu_i);*/
     integrand = -rho_i*mu_i;
     return(integrand);
}

