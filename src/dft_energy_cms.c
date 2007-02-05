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

#include "dft_energy_cms.h"

/****************************************************************************/
double integrand_CMS_freen(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     int icomp,pol_number;

     icomp= iunk-Phys2Unk_first[DENSITY];
     rho_i = x[iunk][inode_box];

                                            /* note second term is old adsorption */ 
     pol_number=0;
     while (Nmer_t[pol_number][icomp]==0) pol_number++;

     if (rho_i > 0.) integrand = 0.5*rho_i*int_stencil(x,inode_box,iunk,THETA_CR_DATA)-rho_i/Nmer[pol_number];

     return(integrand);
}
/****************************************************************************/
double integrand_CMS_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk,sum;
     int icomp,jcomp,pol_number;

     icomp= iunk-Phys2Unk_first[DENSITY];
     rho_bulk = Rho_b[icomp];

     sum=0;
     for (jcomp=0;jcomp<Ncomp;jcomp++){
       sum+=int_stencil_bulk(THETA_CR_DATA,icomp,jcomp,NULL)*Rho_b[jcomp];
     }

     pol_number=0;
     while (Nmer_t[pol_number][icomp]==0) pol_number++;

     integrand = 0.5*rho_bulk*sum-rho_bulk/Nmer[pol_number];
     return(integrand);
}
/****************************************************************************/
