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
 *  FILE: dft_fill_pde_diffusion.c
 *
 *  This file contains the fills of 
 *  diffusion Partial differential equaitons.
 */

#include "dft_fill_pde_diffusion.h"

/****************************************************************************/
double load_nonlinear_transport_eqn(int iunk, int loc_inode, int inode_box,
                                    int *ijk_box, double **x,int resid_only_flag)
{

  int ilist=0, icomp, idim, ijk[3],inode;
  int iln, jln, elem, offset[3];
  int nodes_volm_el, nodes_surf_el, junk2[3];
  int off_ref[2][2] = { {0,1}, {-1,0}};
  int jnode_box,flag,count_flag;
  double wt=1.0, area_0=0.0, area_1=0.0, rho_0, rho_1, mu_0, mu_1, tmp;
  double resid,resid_sum=0.0,mat_val;
  int nodeIndices[8],junk_mu,junk_rho;
  double values[2],values_rho[8],values_mu[8];

  /* pre calc basis function stuff for 3D */

  int igp,numEntries;
  double rho, grad_mu[3], grad_mu_dot_grad_phi;
  static double **phi, ***grad_phi, evol;
  static int first=TRUE;
  if (first) {
    phi = (double **) array_alloc(2, 8, 8, sizeof(double));
    grad_phi = (double ***) array_alloc(3, 8, 8, 3, sizeof(double));
    if (Ndim==3) basis_fn_calc(phi, grad_phi, &evol);
    first = FALSE;
  }

  inode = B2G_node[inode_box];

  /* if we are in a region of constant chemical potential, 
   * solve the equation mu = constant ..... in all other
   * cases loop over the elements near this node and fill the
   * transport equation.
   */

   /* iunk is nodal unknown number, icomp is component number */

   icomp = iunk-Phys2Unk_first[DIFFUSION];

   if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
   else if (Nlists_HW > 2)  ilist = icomp;

   node_to_ijk(node_box_to_node(inode_box),ijk);

   if (Zero_density_TF[inode_box][icomp]) {
       /* set mu to -VEXT_MAX in walls, and skip transport eqn fill */
          resid = x[iunk][inode_box] + VEXT_MAX;
          resid_sum+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (!resid_only_flag){
             mat_val = 1.0;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
          }
   }
   /* check if you are in well-mixed bulk region */
   else if (ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
          resid = x[iunk][inode_box] - Betamu_LBB[icomp];
          resid_sum+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (!resid_only_flag){
             mat_val = 1.0;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                          iunk,inode_box,mat_val);
          }
   }
   else if (Size_x[Grad_dim]-ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
          resid = x[iunk][inode_box] - Betamu_RTF[icomp];
          resid_sum+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (!resid_only_flag){
             mat_val = 1.0;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                         iunk,inode_box,mat_val);
          }
   }
   else {
       count_flag=0;
       for (iln=0; iln< Nnodes_per_el_V; iln++) {

         /* find element with inode_box as local node iln */
         elem = node_to_elem_v2(node_box_to_node(inode_box), iln);

         /* skip fill if element is outside domain or in a wall */
         if (elem >= 0 && (Wall_elems[ilist][el_to_el_box(elem)] == -1 ||
                   Lsemiperm[WallType[Wall_elems[ilist][el_to_el_box(elem)]]][icomp]) ){

           if (Ndim==1) {
             if      (iln == 0) wt =  1.0 / (Esize_x[0] * 6.0);
             else if (iln == 1) wt = -1.0 / (Esize_x[0] * 6.0);

             if (iln == 0) {
               area_0 =  Area_IC[inode_box];
               area_1 =  Area_IC[inode_box+1];
             }
             else {  /*iln = 1*/
               area_0 =  Area_IC[inode_box-1];
               area_1 =  Area_IC[inode_box];
             }
           }

           /* process indices to all local mu and rho unks in this element */
           /* this loop should be good for 2D and 3D as well */

           junk_mu = iunk;
           junk_rho = Phys2Unk_first[DENSITY]+iunk-Phys2Unk_first[DIFFUSION];
           flag=FALSE;
           for (jln=0; jln< Nnodes_per_el_V; jln++) {
             for (idim=0; idim<Ndim; idim++){
                 nodes_surf_el = POW_INT(2,idim);
                 nodes_volm_el = POW_INT(2,idim+1);
                 offset[idim] = 
                      off_ref[((nodes_volm_el+iln)%nodes_volm_el)/nodes_surf_el]
                             [((nodes_volm_el+jln)%nodes_volm_el)/nodes_surf_el];
             }

             jnode_box = offset_to_node_box(ijk_box, offset, junk2);
             if (Zero_density_TF[jnode_box][icomp]) flag=TRUE;

             nodeIndices[jln] = jnode_box;
           }
           if (flag) count_flag++;
           if (count_flag==Nnodes_per_el_V) {
               printf("WARNING : zero density nodes in all elements !!!\n");
               printf("inode_box=%d elem: %d iunk: %d \n",inode_box,elem,iunk);
               exit(-1);
                                             }

           /* only fill if every node is non zero density node - otherwise will be 
              including discontinuities in mu */
           if (flag==FALSE){

              if (Ndim==1) {   

                rho_0 = x[junk_rho][nodeIndices[0]];
                rho_1 = x[junk_rho][nodeIndices[1]];
                mu_0  = x[junk_mu][nodeIndices[0]];
                mu_1  = x[junk_mu][nodeIndices[1]];

                tmp = (2.0*rho_0*area_0 + rho_1*area_0 
                       + rho_0*area_1 + 2.0*rho_1*area_1);

                resid = wt*(mu_0 - mu_1)*tmp;
                resid_sum+=resid;
                dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

                if (!resid_only_flag){ 
                   numEntries=2;
                   values[0]=values[1]=wt*(mu_0 - mu_1)*(2.0*area_0 + area_1);
                   dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                     junk_rho, nodeIndices,values,numEntries);
                   values[0] = wt*tmp;  values[1]=-values[0]; 
                   dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                     junk_mu, nodeIndices,values,numEntries);
                }

                /* add in a bulk convection term */
                if (iln == 0) tmp = (2.0*area_0 + area_1)/6.0;
                else          tmp = (area_0 + 2.0*area_1)/6.0;
                resid = Velocity * (rho_1 - rho_0) * tmp;
                resid_sum+=resid;
                dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                if (!resid_only_flag){
                   values[0] = -Velocity*tmp;  values[1]=-values[0];
                   dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                  junk_rho, nodeIndices,values,numEntries);
                }
              }
              else if (Ndim==3){     /* Ndim==3 case */

                /* loop over all quad points */

                for (igp=0; igp<8; igp++) { 
  
                  /* calculate quantities at quad point */
                  rho = 0.0; 
                  grad_mu[0] = grad_mu[1] = grad_mu[2] = 0.0;
                  for (jln=0; jln<8; jln++) {

                     rho += phi[jln][igp] * x[junk_rho][nodeIndices[jln]];
                     grad_mu[0] += grad_phi[jln][igp][0] * x[junk_mu][nodeIndices[jln]];
                     grad_mu[1] += grad_phi[jln][igp][1] * x[junk_mu][nodeIndices[jln]];
                     grad_mu[2] += grad_phi[jln][igp][2] * x[junk_mu][nodeIndices[jln]];
                  }
                  grad_mu_dot_grad_phi = grad_mu[0]*grad_phi[iln][igp][0]
                                       + grad_mu[1]*grad_phi[iln][igp][1]
                                       + grad_mu[2]*grad_phi[iln][igp][2];
                                      
                  resid = evol * rho * grad_mu_dot_grad_phi;
                  resid_sum+=resid;
                  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                  if (!resid_only_flag){
                     for (jln=0; jln<8; jln++) {
                       values_rho[jln]= evol * phi[jln][igp] * grad_mu_dot_grad_phi;
                       values_mu[jln] = evol * rho * (
                                    grad_phi[jln][igp][0] * grad_phi[iln][igp][0]
                                  + grad_phi[jln][igp][1] * grad_phi[iln][igp][1]
                                  + grad_phi[jln][igp][2] * grad_phi[iln][igp][2]);
                     }
                     numEntries=8;
                     dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                     junk_rho, nodeIndices,values_rho,numEntries);
                     dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                     junk_mu, nodeIndices,values_mu,numEntries);
                  }

                }
              }
              else {
                printf("load_transport not written for 2D yet\n");
                exit(-1);
              }  
           }
         }
       }
   }
   return(resid_sum);
}
/****************************************************************************/
double load_linear_transport_eqn(int iunk,int loc_inode,int inode_box, 
                                 int *ijk_box,double **x,int resid_only_flag)
{

  int iwall, isten, idim,ijk[3],icomp;
  int iln, jln, elem, offset[3], el_box;
  int nodes_volm_el, nodes_surf_el, junk2[3];
  /* static variables keep their value for every time the function is called*/
  static double *wt_lp_1el, *wt_s_1el;
  static int   **elem_permute, off_ref[2][2] = { {0,1}, {-1,0}};
  int  jnode_box;
  double resid, resid_sum=0.0, mat_val;

  /* First time through, load weights appropriate for this Ndim */
  
    set_fem_1el_weights(&wt_lp_1el, &wt_s_1el, &elem_permute);

/* if we are at a boundary node and the boundary condition is constant
 * potential, then solve the equation psi = constant ..... in all other
 * cases loop over the elements near this node and fill the laplace and
 * source terms.  Note that the surface charge boundary condition is
 * handled in the main fill routine. */

   /* iunk is nodal unknown number, icomp is component number */

   icomp = iunk-Phys2Unk_first[DIFFUSION];

   node_to_ijk(node_box_to_node(inode_box),ijk);

   iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
   if (iwall != -1) {
          resid = x[iunk][inode_box] + VEXT_MAX;
          resid_sum+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (!resid_only_flag){
             mat_val = 1.0;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
          }
   }
    if (ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
          resid = x[iunk][inode_box] - Betamu_LBB[icomp];
          resid_sum+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (!resid_only_flag){
             mat_val = 1.0;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                         iunk,inode_box,mat_val);
          }
   }
   else if (Size_x[Grad_dim]-ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
          resid = x[iunk][inode_box] - Betamu_RTF[icomp];
          resid_sum+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (!resid_only_flag){
             mat_val = 1.0;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                         iunk,inode_box,mat_val);
          }
   }
   else {
       for (iln=0; iln< Nnodes_per_el_V; iln++) {

         elem = node_to_elem_v2(node_box_to_node(inode_box), iln);
         el_box = el_to_el_box(elem);

         if (elem >= 0 ) {

           for (jln=0; jln< Nnodes_per_el_V; jln++) {

            /* 
             * elem_permute takes the current local node and another local
             * node in the element, and aligns it with a reference
             * element where the weights were calculated
             */

             isten = elem_permute[iln][jln];

             for (idim=0; idim<Ndim; idim++){
                 nodes_surf_el = POW_INT(2,idim);
                 nodes_volm_el = POW_INT(2,idim+1);
                 offset[idim] = 
                      off_ref[((nodes_volm_el+iln)%nodes_volm_el)/nodes_surf_el]
                             [((nodes_volm_el+jln)%nodes_volm_el)/nodes_surf_el];
             }

             jnode_box = offset_to_node_box(ijk_box, offset, junk2);

             /*
              * add in Laplace term
              */
             if (jnode_box >= 0){  /* new flag for boundaries */
                   resid = wt_lp_1el[isten]*x[iunk][jnode_box];
                   resid_sum+=resid;
                   dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                   if (!resid_only_flag){
                      mat_val = wt_lp_1el[isten];
                      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                iunk,jnode_box,mat_val);
                   }
             }

           }     /* end of loop over all local nodes in this fluid element */
         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */
   }
   return(resid_sum);
}
/****************************************************************************/
