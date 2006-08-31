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

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 *  FILE: dft_fill_mf.c
 *
 *  This file contains the fill of the residual equations and Jacobian
 *  matrix for mean field attractions and electrostatics.
 */

#include "dft_globals_const.h"
#include "rf_allo.h"


double load_mean_field(int sten_type, int iunk, int loc_inode,  
                     int icomp, int izone, int *ijk_box, 
                     double **x, int resid_only_flag)
{
   double resid_sum;
   int inode_box,jzone_flag;
   inode_box=L2B_node[loc_inode];

   jzone_flag=FALSE;

   resid_and_Jac_sten_fill_sum_Ncomp(sten_type,x,iunk,loc_inode,inode_box,izone,
                     ijk_box,resid_only_flag,jzone_flag,
                     NULL, &resid_rho_bar,&jac_rho_bar);
   resid_sum=Temporary_sum;
   return(resid_sum);
}
