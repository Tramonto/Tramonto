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
double integrand_ideal_gas_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i;

     rho_i = x[iunk][inode_box];

     if (rho_i > DENSITY_MIN)
          integrand = rho_i*(log(rho_i)-1.0);
     else integrand=0.0;
     return(integrand);
}
/****************************************************************************/
double integrand_ideal_gas_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;

     if (Type_poly==WTC) rho_bulk = Rho_seg_b[iunk-Phys2Unk_first[DENSITY]];
     else                rho_bulk = Rho_b[iunk-Phys2Unk_first[DENSITY]];

     integrand = rho_bulk*(log(rho_bulk)-1.0);
     return(integrand);
}
/****************************************************************************/
