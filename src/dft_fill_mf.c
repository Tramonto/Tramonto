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

/*
 *  FILE: dft_fill_mf.c
 *
 *  This file contains the fill of the residual equations and Jacobian
 *  matrix for mean field attractions and electrostatics.
 */

#include "dft_fill_mf.h"


double load_mean_field(int sten_type, int iunk, int loc_inode,  
                     int icomp, int izone, int *ijk_box, 
                     double **x, int resid_only_flag)
{
   double resid_sum=0.0,mat_val;
   int inode_box,jzone_flag;
   inode_box=L2B_node[loc_inode];

   jzone_flag=FALSE;

   if (Type_attr !=NONE && ATTInA22Block==FALSE && Unk2Phys[iunk]==MF_EQ){  /* in this case we fill attractions as independent variables - this is the
                                    first step to setting these matrix coefficients to be constant.  otherwise, the
                                    mean field terms are just summed into the Euler-Lagrange Equation. */
           if (resid_only_flag != INIT_GUESS_FLAG){
              resid_sum = -x[iunk][inode_box];
              mat_val=-1.0;
              if (resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_sum);
              if (resid_only_flag==FALSE) dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
           }
   }

   resid_sum+=resid_and_Jac_sten_fill_sum_Ncomp(sten_type,x,iunk,loc_inode,inode_box,izone,
                     ijk_box,resid_only_flag,jzone_flag,
                     NULL, &resid_rho_bar,&jac_rho_bar);
   return(resid_sum);
}
