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
 *  FILE: dft_guess_CHEMPOT.c
 *
 *  This file contains a routine to set an initial guess for the
 *  chemical potential variable (diffusion problems).
 *
 */

#include "dft_guess_CHEMPOT.h"
 
/************************************************************/
/* setup_chem_pot: for cases with steady state profiles,
   set up regions of constant (electro)chemical potentials */
void setup_chem_pot(double **xOwned)
{
  int loc_inode,inode,ijk[3],icomp,iunk,i,ipol,nloop;
  double x_dist,x_tot;

  if (Type_poly==NONE)nloop=Ncomp;
  else{
     if (Lseg_densities) nloop=Nseg_tot;
     else                nloop=Npol_comp;
  }

  x_tot = Size_x[Grad_dim]-2.*X_const_mu;
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode     = L2G_node[loc_inode];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[Grad_dim]*ijk[Grad_dim]-X_const_mu;

     for (i=0; i<nloop; i++){
        if (Restart==RESTART_FEWERCOMP && i<nloop-Nmissing_densities) i=nloop-Nmissing_densities;
        iunk = i+Phys2Unk_first[DIFFUSION];
        if (Type_poly==NONE)icomp=i;
        else{
           if (Lseg_densities) icomp=Unk2Comp[i];
           ipol=i;
        }
 
        if (Type_poly==NONE && Zero_density_TF[L2B_node[loc_inode]][icomp]==TRUE) xOwned[iunk][loc_inode] = -VEXT_MAX;
        else{
           if (Type_poly==NONE){
              if (Ipot_ff_c == 1){
                xOwned[iunk][loc_inode] = log(xOwned[Phys2Unk_first[DENSITY]+i][loc_inode])
                           + Charge_f[icomp]*(xOwned[Phys2Unk_first[POISSON]][loc_inode]);

              }
              else{
                  if (x_dist<0.)           xOwned[iunk][loc_inode]=Betamu_LBB[i];
                  else if (x_dist > x_tot) xOwned[iunk][loc_inode]=Betamu_RTF[i];
                  else  xOwned[iunk][loc_inode] = Betamu_LBB[i] + (Betamu_RTF[i]-Betamu_LBB[i])* x_dist/x_tot;
              }
           }
           else{
                  if (x_dist<0.)           xOwned[iunk][loc_inode]=Betamu_chain_LBB[ipol];
                  else if (x_dist > x_tot) xOwned[iunk][loc_inode]=Betamu_chain_RTF[ipol];
                  else  xOwned[iunk][loc_inode] = Betamu_chain_LBB[ipol] + 
                              (Betamu_chain_RTF[ipol]-Betamu_chain_LBB[ipol])* x_dist/x_tot;
           }
       }
     }
  }
  return;
}
/************************************************************/
/*calc_init_chem_pot: calculate the chemical potential initial guess based on a density profile and
  the diffusion residual equation. */
void calc_init_chem_pot(double **xInBox,double **xOwned)
{
  int loc_inode,inode_box,ijk_box[3],iunk;
  double resid_DIFFUSION;

  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, xInBox);  /* get owned values to box values */

  for (iunk=Phys2Unk_first[DIFFUSION]; iunk<Phys2Unk_last[DIFFUSION];iunk++){
     for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
        inode_box=L2B_node[loc_inode];
        node_box_to_ijk_box(inode_box, ijk_box);

        resid_DIFFUSION=load_nonlinear_transport_eqn(iunk, loc_inode, inode_box, ijk_box,xInBox,INIT_GUESS_FLAG);
        xInBox[iunk][inode_box]=resid_DIFFUSION;
        xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
     }
  }
  return;
}
/************************************************************/


