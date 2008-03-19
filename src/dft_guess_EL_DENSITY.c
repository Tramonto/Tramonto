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
 *  FILE: dft_guess_ELDENSITY.c
 *
 *  This file contains various  possible initial guess strategies
 *  for the density variables in a FMT based perturbation DFT.
 *
 */

#include "dft_guess_EL_DENSITY.h"

void setup_density(double **xInBox,int iguess)
{
    switch(iguess){
      case CONST_RHO:
            if (Lseg_densities){
                 setup_const_density(xInBox,Rho_seg_b,Nseg_tot,0);
            }
            else              setup_const_density(xInBox,Rho_b,Ncomp,0);
            break;

      case CONST_RHO_V:
            setup_const_density(xInBox,Rho_coex,1,1);
            break;
      case CONST_RHO_L:
            setup_const_density(xInBox,Rho_coex,1,0);
            break;

      case EXP_RHO:
            if (Lseg_densities){
                 setup_exp_density(xInBox,Rho_seg_b,Nseg_tot,0);
            }
            else setup_exp_density(xInBox,Rho_b,Ncomp,0);
            break;
      case EXP_RHO_V:
            setup_exp_density(xInBox,Rho_coex,1,1);
            break;
      case EXP_RHO_L:
            setup_exp_density(xInBox,Rho_coex,1,0);
            break;

      case STEP_PROFILE:
            setup_stepped_profile(xInBox);
            break;

      case CHOP_RHO_L:
            setup_exp_density(xInBox,Rho_coex,1,1);
            chop_profile(xInBox,iguess);
            break;
      case CHOP_RHO_V:
            setup_exp_density(xInBox,Rho_coex,1,0);
            chop_profile(xInBox,iguess);
            break;
      case LINEAR:
            setup_linear_profile(xInBox);
    }  /* end of iguess switch */
    return;
}
/*********************************************************/
/*setup_const_density: in this routine set up a constant
        density profile wherever Zero_density_TF = FALSE */
void setup_const_density(double **xInBox, double *rho,int nloop,int index)
{
  int i,inode_box,iunk,zeroTF;

  for (inode_box=0; inode_box<Nnodes_box; inode_box++){
     for (i=0; i<nloop; i++){
	 iunk = i+Phys2Unk_first[DENSITY];
         if (Lseg_densities) zeroTF=Zero_density_TF[inode_box][Unk2Comp[i]];
         else                zeroTF=Zero_density_TF[inode_box][i];
         if (!zeroTF){
            if (nloop > 1) xInBox[iunk][inode_box] = rho[i];
            else           xInBox[iunk][inode_box] = rho[index];
         }
         else xInBox[iunk][inode_box] = 0.0;
      }
  }
  return;
}
/*********************************************************/
/*setup_stepped_profile: in this routine set up a stepped
        density profile wherever Zero_density_TF = FALSE */
void setup_stepped_profile(double **xInBox)
{
  int loc_inode,i,j,nloop,inode_box,iunk,icomp,inode;
  double nodepos[3];

  if (Lseg_densities) nloop=Nseg_tot;
  else                nloop=Ncomp;

  for (inode_box=0; inode_box<Nnodes_box; inode_box++){
     inode = B2G_node[inode_box];
     node_to_position(inode,nodepos);

     for (i=0;i<Nsteps;i++){
       if (nodepos[Orientation_step[i]] >= Xstart_step[i] &&
           nodepos[Orientation_step[i]] <= Xend_step[i]){

           for (j=0; j<nloop; j++){
               if (Lseg_densities) icomp=Unk2Comp[j];
               else                icomp=j;
	       iunk = Phys2Unk_first[DENSITY]+j;
               if (!Zero_density_TF[inode_box][icomp]){
                   xInBox[iunk][inode_box]=Rho_step[icomp][i];
               }
               else xInBox[iunk][inode_box]=0.0;
           }
       }
    }
  }

  if (Lsteady_state){
     printf("stepped profile not set up for chemical potentials at this time\n");
     exit(-1);
  }
  return;
}
/*********************************************************/
/*setup_exp_density: in this routine set up a density 
                     profile as rho_b*exp(-Vext/kT)*/
void setup_exp_density(double **xInBox, double *rho,int nloop,int index)
{

  int loc_inode,i,inode_box,iunk,icomp;

  for (inode_box=0; inode_box<Nnodes_box; inode_box++){
     loc_inode=B2L_node[inode_box];
     for (i=0; i<nloop; i++) {

        if (Lseg_densities) icomp=Unk2Comp[i];
        else icomp=i;

	iunk = i+Phys2Unk_first[DENSITY];
        if (!Zero_density_TF[inode_box][icomp]){
           if (loc_inode >=0){
               if (nloop > 1) xInBox[iunk][inode_box] = rho[i] * exp(-Vext[loc_inode][icomp]);
               else           xInBox[iunk][inode_box] = rho[index] * exp(-Vext[loc_inode][icomp]);
           }
           else{ /* since Vext is not indexed on the box at this time, we need to approximate here */
               if (nloop > 1) xInBox[iunk][inode_box] = rho[i];
               else           xInBox[iunk][inode_box] = rho[index];
           }
        }
        else xInBox[iunk][inode_box] = 0.0;

     }
  }
  return;
}
/************************************************************/
/*setup_step_2consts: density profile where the left (LBB)
      and right (RTF) sides of the box have different bulk
      densities.  This can either be for chem.pot. gradients
      or for a liquid-vapor profile.  */

void setup_step_2consts(double **xInBox)
{

  int icomp,inode_box,ijk[3],inode;
  int iunk,nloop,i;
  double x_dist;

  if (Lseg_densities) nloop=Nseg_tot;
  else                nloop=Ncomp;

  for (inode_box=0; inode_box<Nnodes_box; inode_box++){
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[0]*ijk[0] - 0.5*Size_x[0];

     for (i=0; i<nloop; i++) {
        if (Lseg_densities) icomp=Unk2Comp[i];
        else                icomp=i;

        iunk = i+Phys2Unk_first[DENSITY];
        if (!Zero_density_TF[inode_box][icomp]){
           if (x_dist< 0.0)  xInBox[iunk][inode_box] = constant_boundary(iunk,-3);
           else              xInBox[iunk][inode_box] = constant_boundary(iunk,-4);
        }
        else xInBox[iunk][inode_box] = 0.0;

     }
  }
  return;
}
/************************************************************/
/*setup_linear profile: density profile that interpolates 
      linearly between densities at the left (LBB) and 
      right (RTF) sides of the box. This can either be for 
      chem.pot. gradients
      or for a liquid-vapor profile.  */

void setup_linear_profile(double **xInBox)
{

  int icomp,inode_box,ijk[3],inode;
  int iunk,i,nloop;
  double x_dist,x_tot,rho_LBB, rho_RTF;
  
  if (Lseg_densities) nloop=Nseg_tot;
  else                nloop=Ncomp;

  if (Lsteady_state) x_tot=Size_x[Grad_dim]-2.0*X_const_mu;
  else x_tot = Size_x[Grad_dim];

  for (inode_box=0; inode_box<Nnodes_box; inode_box++){
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[0]*ijk[0]-X_const_mu;

     for (i=0; i<nloop; i++) {
        if (Lseg_densities) icomp=Unk2Comp[i];
        else                icomp=i;

        iunk = Phys2Unk_first[DENSITY]+i;
        if (!Zero_density_TF[inode_box][icomp]){
           rho_LBB = constant_boundary(iunk,-3);
           rho_RTF = constant_boundary(iunk,-4);

           if (x_dist <0.) xInBox[iunk][inode_box] = rho_LBB;
           else if (x_dist > x_tot) xInBox[iunk][inode_box] = rho_RTF;
           else  xInBox[iunk][inode_box] = rho_LBB + (rho_RTF-rho_LBB)*x_dist/x_tot;
        }
        else xInBox[iunk][inode_box] = 0.0;

     }
  }
  return;
}
/****************************************************************************/
