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
 *  FILE: dft_guess_ELDENSITY.c
 *
 *  This file contains various  possible initial guess strategies
 *  for the density variables in a FMT based perturbation DFT.
 *
 */

/*#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
#include <string.h>*/

#include "dft_guess_EL_DENSITY.h"

void setup_density(double **xOwned,int iguess)
{
    switch(iguess){
      case CONST_RHO:
            if (Lseg_densities){
                 setup_const_density(xOwned,Rho_seg_b,Nseg_tot,0);
            }
            else              setup_const_density(xOwned,Rho_b,Ncomp,0);
            break;

      case CONST_RHO_V:
            setup_const_density(xOwned,Rho_coex,1,1);
            break;
      case CONST_RHO_L:
            setup_const_density(xOwned,Rho_coex,1,0);
            break;

      case EXP_RHO:
            setup_exp_density(xOwned,Rho_b,Ncomp,0);
            break;
      case EXP_RHO_V:
            setup_exp_density(xOwned,Rho_coex,1,1);
            break;
      case EXP_RHO_L:
            setup_exp_density(xOwned,Rho_coex,1,0);
            break;

      case STEP_PROFILE:
            setup_stepped_profile(xOwned);

      case CHOP_RHO_L:
            setup_exp_density(xOwned,Rho_coex,1,1);
            chop_profile(xOwned,iguess);
            break;
      case CHOP_RHO_V:
            setup_exp_density(xOwned,Rho_coex,1,0);
            chop_profile(xOwned,iguess);
            break;
      case LINEAR:
            setup_linear_profile(xOwned);
    }  /* end of iguess switch */
    return;
}
/*********************************************************/
/*setup_const_density: in this routine set up a constant
        density profile wherever Zero_density_TF = FALSE */
void setup_const_density(double **xOwned, double *rho,int nloop,int index)
{
  int loc_inode,i,inode_box,iunk,zeroTF;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (i=0; i<nloop; i++){
	 iunk = i+Phys2Unk_first[DENSITY];
         if (Lseg_densities) zeroTF=Zero_density_TF[inode_box][Unk2Comp[i]];
         else                zeroTF=Zero_density_TF[inode_box][i];
         if (!zeroTF){
            if (nloop > 1) xOwned[iunk][loc_inode] = rho[i];
            else           xOwned[iunk][loc_inode] = rho[index];
         }
         else xOwned[iunk][loc_inode] = 0.0;

         /* set up initial guess if chemical potential is an unknown */
         if (Lsteady_state) {
	     iunk = i+Phys2Unk_first[DIFFUSION];
             if (!zeroTF){
                if (nloop > 1) xOwned[iunk][loc_inode] = Betamu[i];
                else           xOwned[iunk][loc_inode] = Betamu[index];
             }
             else xOwned[iunk][loc_inode] = -10.0; /* zero density == -infinite chem.pot*/
         }
      }
  }
  return;
}
/*********************************************************/
/*setup_stepped_profile: in this routine set up a stepped
        density profile wherever Zero_density_TF = FALSE */
void setup_stepped_profile(double **xOwned)
{
  int loc_inode,i,inode_box,iunk,icomp,inode;
  double nodepos[3];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode = L2G_node[loc_inode];
     node_to_position(inode,nodepos);

     for (i=0;i<Nsteps;i++){
       if (nodepos[Orientation_step[i]] >= Xstart_step[i] &&
           nodepos[Orientation_step[i]] <= Xend_step[i]){
           for (icomp=0; icomp<Ncomp; icomp++){
	       iunk = Phys2Unk_first[DENSITY]+icomp;
               if (!Zero_density_TF[inode_box][icomp]){
                   xOwned[iunk][loc_inode]=Rho_step[icomp][i];
               }
               else xOwned[iunk][loc_inode]=0.0;
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
void setup_exp_density(double **xOwned, double *rho,int nloop,int index)
{

  int loc_inode,i,inode_box,iunk;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (i=0; i<nloop; i++) {

	iunk = i+Phys2Unk_first[DENSITY];
        if (!Zero_density_TF[inode_box][i]){
            if (nloop > 1) xOwned[iunk][loc_inode] = rho[i] * exp(-Vext[loc_inode][i]);
            else           xOwned[iunk][loc_inode] = rho[index] * exp(-Vext[loc_inode][i]);
        }
        else xOwned[iunk][loc_inode] = 0.0;

        if (xOwned[iunk][loc_inode]>1.0) xOwned[iunk][loc_inode] = 0.99;   /* may need this to make sure that rb3<1 always */


         /* set up initial guess if chemical potential is an unknown */
         if (Lsteady_state) {
	    iunk = i+Phys2Unk_first[DIFFUSION];
             if (!Zero_density_TF[inode_box][i]){
                if (nloop > 1) xOwned[iunk][loc_inode] = Betamu[i];
                else           xOwned[iunk][loc_inode] = Betamu[index];
             }
             else xOwned[iunk][loc_inode] = -10.0; /* zero density == -infinite chem.pot*/
         }
     }
  }
  return;
}
/************************************************************/
/*setup_step_2consts: density profile where the left (LBB)
      and right (RTF) sides of the box have different bulk
      densities.  This can either be for chem.pot. gradients
      or for a liquid-vapor profile.  */

void setup_step_2consts(double **xOwned)
{

  int loc_inode,icomp,inode_box,ijk[3],inode;
  int iunk;
  double x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[0]*ijk[0] - 0.5*Size_x[0];

     for (icomp=0; icomp<Ncomp; icomp++) {

        iunk = icomp+Phys2Unk_first[DENSITY];
        if (!Zero_density_TF[inode_box][icomp]){
           if (x_dist< 0.0)  xOwned[iunk][loc_inode] = constant_boundary(iunk,-3);
           else              xOwned[iunk][loc_inode] = constant_boundary(iunk,-4);
        }
        else xOwned[iunk][loc_inode] = 0.0;

       /* set up initial guess if chemical potential is an unknown */
        if (Lsteady_state) {
            iunk = icomp+Phys2Unk_first[DIFFUSION];
            if (!Zero_density_TF[inode_box][icomp]){
               if (x_dist< 0.0) xOwned[iunk][loc_inode] = constant_boundary(iunk,-3);
               else             xOwned[iunk][loc_inode] = constant_boundary(iunk,-4);
            }
            else xOwned[iunk][loc_inode] = -10.0;  /*zero density == -infinite chem.pot*/
        }
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

void setup_linear_profile(double **xOwned)
{

  int loc_inode,icomp,inode_box,ijk[3],inode;
  int iunk;
  double x_dist,x_tot,rho_LBB, rho_RTF;

  if (Lsteady_state) x_tot=Size_x[Grad_dim]-2.0*X_const_mu;
  else x_tot = Size_x[Grad_dim];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[0]*ijk[0]-X_const_mu;

     for (icomp=0; icomp<Ncomp; icomp++) {

        iunk = Phys2Unk_first[DENSITY]+icomp;
        if (!Zero_density_TF[inode_box][icomp]){
           rho_LBB = constant_boundary(iunk,-3);
           rho_RTF = constant_boundary(iunk,-4);

           if (x_dist <0.) xOwned[iunk][loc_inode] = rho_LBB;
           else if (x_dist > x_tot) xOwned[iunk][loc_inode] = rho_RTF;
           else  xOwned[iunk][loc_inode] = rho_LBB + (rho_RTF-rho_LBB)*x_dist/x_tot;
        }
        else xOwned[iunk][loc_inode] = 0.0;

     }
  }
  return;
}
/****************************************************************************/
