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

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

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
