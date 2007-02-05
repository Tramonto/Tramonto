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
 *  FILE: dft_switch_FMT.c
 *
 *  This file contains switches to select among varous fundamental measures
 *  theories (FMTs) implemented in the code.  These switches provide an 
 *  interface between core logic for filling the matrix problem and the
 *  specific physics implementations.
 *
 */

#include "dft_switch_FMT.h"
/****************************************************************************/
/*phispt_switch: Logic controlling type of FMT functional used.             */
double phispt_switch(double *n)
{ 
  double FMTphispt;
  if (n[3] < 1.0 && n[2] > 0.0 && n[3]>1.e-10){
     switch(Type_func){
         case FMT1:
             FMTphispt=FMT1_energy_density(n);
             break;
         case FMT2:
             FMTphispt=FMT2_energy_density(n);
             break;
         case FMT3:
             FMTphispt=FMT3_energy_density(n);
             break;
         default:
             printf("problem with type of HS FUNCTIONAL: phispt_switch ");
             exit(-1); break;
     }
     return(FMTphispt);
  }
  else return(0.0);
}
/****************************************************************************/
/*FMT1stDeriv_switch: Logic controlling type of FMT functional used.             */
void FMT1stDeriv_switch(int inode_box, double **x, struct RB_Struct *dphi_drb)
{ 
  switch(Type_func){
     case FMT1:
        for (inode_box=0;inode_box<Nnodes_box; inode_box++)
          calc_FMT_derivatives(&FMT1_1stderiv,inode_box,x,dphi_drb);
        break;
     case FMT2: 
        for (inode_box=0;inode_box<Nnodes_box; inode_box++)
          calc_FMT_derivatives(&FMT2_1stderiv,inode_box,x,dphi_drb);
        break;
     case FMT3: 
        for (inode_box=0;inode_box<Nnodes_box; inode_box++)
          calc_FMT_derivatives(&FMT3_1stderiv,inode_box,x,dphi_drb);
        break;
      default:
         printf("problem with type of HS FUNCTIONAL: FMT1stDeriv_switch");
         exit(-1); break;
  }
  return;
}
/****************************************************************************/
/*FMT1stDerivBulk_switch: Logic controlling type of FMT functional used in bulk thermo calculations.     */
void FMT1stDerivBulk_switch(double *n,double *inv_n3, double *dphi_drb)
{ 
  double DOT_12,DOT_22;
  DOT_12=0.0; DOT_22=0.0;  /* vector terms are zero in bulk */
  switch(Type_func){
     case FMT1:
        FMT1_1stderiv(n,DOT_12,DOT_22,inv_n3,dphi_drb);
        break;
     case FMT2: 
        FMT2_1stderiv(n,DOT_12,DOT_22,inv_n3,dphi_drb);
        break;
     case FMT3: 
        FMT3_1stderiv(n,DOT_12,DOT_22,inv_n3,dphi_drb);
        break;
      default:
         printf("problem with type of HS FUNCTIONAL: FMT1stDeriv_switch");
         exit(-1); break;
  }
  return;
}
/****************************************************************************/
struct RB_Struct FMT2ndDerivDelta_switch(double *n, int *offset, double *sign, int icomp)
{ 
   struct  RB_Struct tmp;

   switch(Type_func){
      case FMT1:  tmp = d2phi_drb2_delta_rb_FMT1(n,offset,sign,icomp); break;
      case FMT2:  tmp = d2phi_drb2_delta_rb_FMT2(n,offset,sign,icomp); break;
      case FMT3:  tmp = d2phi_drb2_delta_rb_FMT3(n,offset,sign,icomp); break;
      default:
         printf("problem with type of HS FUNCTIONAL: FMT2ndDerivDelta_switch");
         exit(-1); break;
   }
   return(tmp);
}
/****************************************************************************/
struct RB_Struct FMT2ndDerivTheta_switch(double *n)
{ 
   struct  RB_Struct tmp;

   switch(Type_func){
      case FMT1:  tmp = d2phi_drb2_theta_rb_FMT1(n); break;
      case FMT2:  tmp = d2phi_drb2_theta_rb_FMT2(n); break;
      case FMT3:  tmp = d2phi_drb2_theta_rb_FMT3(n); break;
      default:
         printf("problem with type of HS FUNCTIONAL: FMT2ndDerivTheta_switch");
         exit(-1); break;
   }
   return(tmp);
}
/****************************************************************************/

