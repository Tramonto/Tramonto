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
 *  FILE: dft_fill_pde_elec.c
 *
 *  This file contains the matrix fills for 
 *  for electrostatics partial differential equations.
 */

#include "dft_fill_pde_elec.h"

/*****************************************************************************/
double load_poisson_control(int iunk, int loc_inode, int inode_box, int *ijk_box, double **x,int resid_only_flag)
{
  int idim,l_elec_RTF=FALSE,l_elec_LBB=FALSE;
  double nodepos[3],resid_poisson,resid,mat_value;

  if (Type_interface != UNIFORM_INTERFACE){
     node_to_position(node_box_to_node(inode_box),nodepos);
     idim = Grad_dim;
     if ( nodepos[idim] +0.5*Size_x[idim] - X_const_mu <= 0.00000001) l_elec_LBB=TRUE;
     else if (nodepos[idim] - 0.5*Size_x[idim] + X_const_mu >=-0.00000001) l_elec_RTF=TRUE;
  }

  if (l_elec_LBB){
     resid = -Elec_pot_LBB;
     if(resid_only_flag !=INIT_GUESS_FLAG){
       resid += x[iunk][inode_box];
       if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
       if (resid_only_flag==FALSE){
         mat_value = 1.0;
         if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_value;
         dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);
       }
     }
     resid_poisson = resid;
  }
  else if (l_elec_RTF){
     resid = -Elec_pot_RTF;
     if (resid_only_flag !=INIT_GUESS_FLAG){
        resid += x[iunk][inode_box];
        if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
        if (resid_only_flag==FALSE){
          mat_value = 1.0;
          if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_value;
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);
        }
     }
     resid_poisson = resid;
  }
  else if (!l_elec_LBB && !l_elec_RTF) {
     if(Type_coul==POLARIZE) resid_poisson=load_polarize_poissons_eqn(iunk,loc_inode, inode_box,ijk_box,x,resid_only_flag);
     else                    resid_poisson=load_poissons_eqn(iunk,loc_inode,inode_box,ijk_box,x,resid_only_flag);
  }
  if (resid_only_flag==INIT_GUESS_FLAG) return(-resid_poisson);
  else                                  return resid_poisson;
}
/****************************************************************************/
double load_polarize_poissons_eqn(int iunk, int loc_inode, int inode_box, int *ijk_box, double **x,int resid_only_flag)
{

  int iwall,  isten, icomp,ilist,idim,junk,junkP,jtmp;
  int iln, jln, elem, offset[3], el_box;
  int nodes_volm_el, nodes_surf_el, junk2[3];
  int in_wall,numEntries, nodeIndices[2];
  double values[2];
  double pol_wall,wt,charge_i;
  double resid,resid_sum=0.0,mat_val,prefactor_sum=0.0;

  /* static variables keep their value for every time the function is called*/
  static double *wt_lp_1el, *wt_s_1el;
  static int   **elem_permute, off_ref[2][2] = { {0,1}, {-1,0}};
  double rho_0, rho_1, psi_0, psi_1, tmp;
  int jnode_box,junk_rho;
  int reflect_flag[3];

  /* First time through, load weights appropriate for this Ndim */

  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;
  set_fem_1el_weights(&wt_lp_1el, &wt_s_1el, &elem_permute);

/* if we are at a boundary node and the boundary condition is constant
 * potential, then solve the equation psi = constant ..... in all other
 * cases loop over the elements near this node and fill the laplace and
 * source terms.  Note that the surface charge boundary condition is
 * handled in the main fill routine. */

   iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
   if (iwall != -1 && Type_bc_elec[WallType[iwall]] == CONST_POTENTIAL) {
       resid = - Elec_param_w[iwall];
       if (resid_only_flag != INIT_GUESS_FLAG){
          resid += x[iunk][inode_box];
          if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (resid_only_flag==FALSE){
            mat_val = 1.0;
            if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
            dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk,inode_box,mat_val);
          }
       }
       resid_sum+=resid;
   }
   else {
       prefactor_sum=0.0;
       for (iln=0; iln< Nnodes_per_el_V; iln++) {

         elem = node_to_elem(node_box_to_node(inode_box), iln,reflect_flag);
         el_box = el_to_el_box(elem);

         if (elem >= 0 ) {

           if (Vol_charge_flag) {
               resid = -Charge_vol_els[el_box]*4.0*PI/(Temp_elec * Nnodes_per_el_V);
               resid_sum+=resid;
               if (resid_only_flag != INIT_GUESS_FLAG && 
                   resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
           }

           junkP = Phys2Unk_first[POISSON];
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
                in_wall=TRUE;
                for (icomp=0; icomp<Ncomp; icomp++){
                   if (!Zero_density_TF[jnode_box][icomp])in_wall=FALSE;
                }
                if(resid_only_flag != INIT_GUESS_FLAG || junkP!= iunk || jnode_box != inode_box){
                   if (in_wall){
                      /* pol_wall is a constant wall polarization.  set equal
                         to zero, but leave code as place holder. */
                      /*pol_wall=(KAPPA_H2O-1.0);*/
                      pol_wall=0.0;
                      resid = (1.0+pol_wall)*wt_lp_1el[isten]*x[junkP][jnode_box];
                      mat_val = (1.0+pol_wall)*wt_lp_1el[isten];
                   }
                   else{
                      resid = wt_lp_1el[isten]*x[junkP][jnode_box];
                      mat_val = wt_lp_1el[isten];
                   }
                   if (resid_only_flag != INIT_GUESS_FLAG && 
                       resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                   if (resid_only_flag == FALSE){
                       if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[junkP]*Nnodes]+=mat_val;
                       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junkP,jnode_box,mat_val);
                   }
                   resid_sum+=resid;
                   }
                else {
                   if (in_wall) prefactor_sum +=(1.0+pol_wall)*wt_lp_1el[isten];
                   else prefactor_sum+=wt_lp_1el[isten];
                }
             }

             /* 
              * add in source term (electroneutrality sum) 
              */
             for (junk=Phys2Unk_first[DENSITY]; junk<Phys2Unk_last[DENSITY]; junk++){
               if (Lseg_densities) icomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
               else                icomp=junk-Phys2Unk_first[DENSITY];
 
               if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
               else ilist = icomp;
               if (elem !=-2 && (Wall_elems[ilist][el_box] == -1 ||
                   Lsemiperm[WallType[Wall_elems[ilist][el_box]]][icomp]) ){

                  if (Charge_f[icomp] != 0.0) {
                    resid = -wt_s_1el[isten]*KAPPA_H2O*Charge_f[icomp]*x[junk][jnode_box]*4.0*PI/Temp_elec;
                    if (resid_only_flag != INIT_GUESS_FLAG && 
                        resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                    if (resid_only_flag==FALSE){
                       mat_val = -wt_s_1el[isten]*KAPPA_H2O*Charge_f[icomp]*4.0*PI/Temp_elec;
                       if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[junk]*Nnodes]+=mat_val;
                       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                     junk,jnode_box,mat_val);
                    }
                    resid_sum+=resid;
                  } /* End if (Charge_f[icomp] != 0.0) */
               }
             }   /* end of icomp (source term) loop */

           }     /* end of loop over all local nodes in this fluid element */

           /*
            * Now handle the nonlinear polarization contribution
            */

           /*
            * Set weights for polarization term and
            * set indices to all local rhos and psis (elec. pot.) in this element
            */
           if (Ndim==1) { 
              if      (iln == 0) wt =  1.0 / (Esize_x[0] * 2.0);
              else if (iln == 1) wt = -1.0 / (Esize_x[0] * 2.0);
           }

          for (junk_rho=Phys2Unk_first[DENSITY]; junk_rho<Phys2Unk_last[DENSITY]; junk_rho++){
               if (Lseg_densities) icomp=Unk2Comp[junk_rho-Phys2Unk_first[DENSITY]];
               else                icomp=junk_rho-Phys2Unk_first[DENSITY];

               if (!Zero_density_TF[inode_box][icomp] && Lpolarize[icomp]){

               for (jln=0; jln< Nnodes_per_el_V; jln++) {
                  for (idim=0; idim<Ndim; idim++){
                      nodes_surf_el = POW_INT(2,idim);
                      nodes_volm_el = POW_INT(2,idim+1);
                      offset[idim] =
                           off_ref[((nodes_volm_el+iln)%nodes_volm_el)/nodes_surf_el]
                                  [((nodes_volm_el+jln)%nodes_volm_el)/nodes_surf_el];
                  }
                  nodeIndices[jln] = offset_to_node_box(ijk_box, offset, junk2);
               }

               in_wall=TRUE;
               for (icomp=0; icomp<Ncomp; icomp++){
                   if (!Zero_density_TF[jnode_box][icomp])in_wall=FALSE;
               }

               if (!in_wall){

                  if (Ndim==1) {

                    rho_0 = x[junk_rho][nodeIndices[0]];
                    rho_1 = x[junk_rho][nodeIndices[1]];
                    psi_0  = x[junkP][nodeIndices[0]];
                    psi_1  = x[junkP][nodeIndices[1]];

                    tmp = (rho_0 + rho_1);

                    if(resid_only_flag != INIT_GUESS_FLAG || junkP!= iunk || nodeIndices[0] != inode_box){
                       resid = wt*Pol[icomp]*(psi_0)*tmp;
                       resid_sum+=resid;
                       if (resid_only_flag != INIT_GUESS_FLAG && 
                          resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                       
                       if (resid_only_flag==FALSE){
                          numEntries=2;
                          values[0]=values[1]=wt*Pol[icomp]*(psi_0);
                          if (Iwrite_files==FILES_DEBUG_MATRIX){
                             for (jtmp=0; jtmp<numEntries;jtmp++) 
                                   Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[nodeIndices[jtmp]]+Solver_Unk[junk_rho]*Nnodes]+=values[jtmp];
                          }
                          dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                         junk_rho, nodeIndices,values,numEntries);
                          mat_val=wt*Pol[icomp]*tmp; 
                          if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[nodeIndices[0]]+Solver_Unk[junkP]*Nnodes]+=mat_val;
                          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                              junkP, nodeIndices[0],mat_val);
                       }
                    }
                    else{
                       prefactor_sum+=wt*Pol[icomp]*tmp;
                    }

                    if(resid_only_flag != INIT_GUESS_FLAG || junkP!= iunk || nodeIndices[1] != inode_box){
                       resid = wt*Pol[icomp]*(- psi_1)*tmp;
                       resid_sum+=resid;
                       if (resid_only_flag != INIT_GUESS_FLAG && 
                           resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                       if (resid_only_flag==FALSE){
                          numEntries=2;
                          values[0]=values[1]=wt*Pol[icomp]*(-psi_1);
                          if (Iwrite_files==FILES_DEBUG_MATRIX){
                            for (jtmp=0; jtmp<numEntries;jtmp++) 
                                Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[nodeIndices[jtmp]]+Solver_Unk[junk_rho]*Nnodes]+=values[jtmp];
                          }
                          dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                         junk_rho, nodeIndices,values,numEntries);
                          mat_val=-wt*Pol[icomp]*tmp; 
                          if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[nodeIndices[1]]+Solver_Unk[junkP]*Nnodes]+=mat_val;
                          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                         junkP, nodeIndices[1],mat_val);
                       }
                    }
                    else{
                       prefactor_sum+=(-wt*Pol[icomp]*tmp);
                    }
                  }
               }
             } /*end of test for polarizeable fluid species */
           }   /*end of loop over components for adding in polarization term */


         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */
   }

    /* finally address nodes where there is a constant surface charge */
    if (Nodes_2_boundary_wall[Nlists_HW-1][inode_box] != -1){

       iwall     = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
       if (Type_bc_elec[WallType[iwall]] == CONST_CHARGE ){

         charge_i = 0.0;
         for (idim=0; idim<Ndim; idim++) charge_i += Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];
         resid = -4.0*PI*charge_i/Temp_elec;
         if (resid_only_flag != INIT_GUESS_FLAG && 
             resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
         resid_sum +=resid;

       }   /* check for charge b.c. on this wall */
    }      /* end of check for if this is a boundary node */

    if (resid_only_flag==INIT_GUESS_FLAG) resid_sum /= prefactor_sum;
    return(resid_sum);
}
/****************************************************************************/
double load_poissons_eqn(int iunk, int loc_inode, int inode_box, int *ijk_box, double **x,int resid_only_flag)
{

  int iwall,  isten, icomp,ilist,idim;
  int iln, jln, elem, offset[3], el_box,inormal;
  int nodes_volm_el, nodes_surf_el, junk2[3];
  /* static variables keep their value for every time the function is called*/
  static double *wt_lp_1el, *wt_s_1el, **wt_surfsrc_1el;
  static int   **elem_permute;
  int off_ref[2][2]; /*= { {0,1}, {-1,0}};*/
  int jnode_box,junk,junkP,junk2react;
  int reflect_flag[3],ijk[3],flag_bulk_boundary;
  double resid,mat_val,resid_sum=0.0,elec_param_bulk_boundary,prefactor_sum,charge_i;
  double fac_react,fac_react2,alpha_2,pconst,mat_val2;
  int nshift;
  off_ref[0][0]=0;
  off_ref[0][1]=1;
  off_ref[1][0]=-1;
  off_ref[1][1]=0;

  /* First time through, load weights appropriate for this Ndim */

  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;
  ijk_box_to_ijk(ijk_box,ijk);

  
  set_fem_1el_weights(&wt_lp_1el, &wt_s_1el, &elem_permute);
  set_fem_surfsrc_weights(&wt_surfsrc_1el);

/* if we are at a boundary node and the boundary condition is constant
 * potential, then solve the equation psi = constant ..... in all other
 * cases loop over the elements near this node and fill the laplace and
 * source terms.  Note that the surface charge boundary condition is
 * handled in the main fill routine. */

/*  flag_bulk_boundary=FALSE;
  for (idim=0;idim<Ndim;idim++){
      if ( (ijk[idim]==0 && Type_bc[idim][0] == IN_BULK) ||
           (ijk[idim]==Nodes_x[idim]-1 && Type_bc[idim][1] == IN_BULK) ){
            flag_bulk_boundary=TRUE;
            elec_param_bulk_boundary=0.0;
            resid = x[iunk][inode_box] - elec_param_bulk_boundary;
            resid_sum+=resid;
            if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
            if (!resid_only_flag){
               mat_val=1.0;
               if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
               dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk,inode_box,mat_val);
            }
      }
   }
 
   if (flag_bulk_boundary==FALSE){ */

   iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
   if (iwall != -1 && Type_bc_elec[WallType[iwall]] == CONST_POTENTIAL) {
          resid = - Elec_param_w[iwall];
       if (resid_only_flag != INIT_GUESS_FLAG){ 
          resid += x[iunk][inode_box];
          if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          if (resid_only_flag==FALSE){
             mat_val=1.0;
             if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk,inode_box,mat_val);
          }
       }
       resid_sum+=resid;
   }
   else {
       prefactor_sum=0.0;
       for (iln=0; iln< Nnodes_per_el_V; iln++) {

         /*elem = node_to_elem_v2(node_box_to_node(inode_box), iln);*/
         elem = node_to_elem(node_box_to_node(inode_box), iln,reflect_flag);
         el_box = el_to_el_box(elem);

         if (elem >= 0 ) {
             if (Vol_charge_flag){ 
               resid = -Charge_vol_els[el_box]*4.0*PI
                              / (Temp_elec * Nnodes_per_el_V);
               resid_sum += resid;
               if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
           }

           junkP=Phys2Unk_first[POISSON];
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
                 if (resid_only_flag != INIT_GUESS_FLAG || junkP!= iunk || jnode_box != inode_box){
                   resid = Dielec[el_box]*wt_lp_1el[isten]*x[junkP][jnode_box];
                   if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                   if (resid_only_flag==FALSE){
                      mat_val = Dielec[el_box]*wt_lp_1el[isten];
                      if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[junkP]*Nnodes]+=mat_val;
                      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, junkP,jnode_box,mat_val);
                   }
                   resid_sum+=resid;
                 }
                 else{
                    prefactor_sum+=Dielec[el_box]*wt_lp_1el[isten];
                 }
             }

             /* 
              * add in source term (electroneutrality sum) 
              */
             for (junk=Phys2Unk_first[DENSITY]; junk<Phys2Unk_last[DENSITY]; junk++){
               if (Lseg_densities) icomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
               else                icomp=junk-Phys2Unk_first[DENSITY];

               if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
               else ilist = icomp;

               if (elem !=-2 && (Wall_elems[ilist][el_box] == -1 ||
                   Lsemiperm[WallType[Wall_elems[ilist][el_box]]][icomp]) ){

                  if (Charge_f[icomp] != 0.0) {
                       resid = -wt_s_1el[isten]*Charge_f[icomp]*x[junk][jnode_box]*4.0*PI/Temp_elec;
                       if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                       if (resid_only_flag==FALSE){
                          mat_val = -wt_s_1el[isten]*Charge_f[icomp]* 4.0*PI/Temp_elec;
                          if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[junk]*Nnodes]+=mat_val;
                          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, junk,jnode_box,mat_val);
                       }
                       resid_sum+=resid;
                    
                  } /* End if (Charge_f[icomp] != 0.0) */

               }
             }   /* end of icomp (source term) loop */


/* Deal with charge regulating surfaces **/
/* NEW CODE - IMPLEMENTATION STARTED FOR SPATIALLY VARYING CONSTANT CHARGE AND
              SPATIALLY VARYING CHARGE REGULATING SURFACES.  CODE FOR UNIFORM SURFACE CHARGES
              IS FOUND BELOW THESE LOOPS, AND WILL BE USED TO TEST THIS MORE GENERAL CODE AUGUST 2021, LJDF*/
             if (Nodes_2_boundary_wall[Nlists_HW-1][jnode_box] != -1){
                 iwall     = Nodes_2_boundary_wall[Nlists_HW-1][jnode_box];
                 if (Type_bc_elec[WallType[iwall]] == CONST_CHARGE ){
/*****CONST CHARGE --- GOOD TEST CASE*****
         charge_i = 0.0;
         for (idim=0; idim<Ndim; idim++) charge_i += Charge_w_sum_els[B2L_node[jnode_box]][idim]*Area_surf_el[idim];
         inormal=abs(Surf_normal[0][B2L_node[jnode_box]][0])-1;
         resid = -(wt_surfsrc_1el[isten][inormal])*4.0*PI*charge_i/Temp_elec;
         if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
         resid_sum +=resid;
END CONST CHARGE*/
/******REACT SURF TYPE 1********
         charge_i = 0.0;
                     for (idim=0; idim<Ndim; idim++) charge_i += Charge_w_sum_els[B2L_node[jnode_box]][idim];
                     inormal=abs(Surf_normal[0][B2L_node[jnode_box]][0])-1;
                     fac_react=Rho_b[0]*exp(-charge_i+x[iunk][jnode_box]);
                     resid = -(wt_surfsrc_1el[isten][inormal])*(4.0*PI/Temp_elec)*(1./(1.+fac_react));
                     if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                     resid_sum +=resid;

                     if (resid_only_flag==FALSE){
                            mat_val=-(wt_surfsrc_1el[isten][inormal])*(4.0*PI/Temp_elec)*(-fac_react/POW_DOUBLE_INT((1+fac_react),2));
                            if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
                            dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk,jnode_box,mat_val);
                     }
******END REACT SURF TYPE 1********/
/******REACT SURF TYPE 2********
             charge_i = 0.0;
             for (idim=0; idim<Ndim; idim++) charge_i += Charge_w_sum_els[loc_inode][idim];
             inormal=abs(Surf_normal[0][B2L_node[jnode_box]][0])-1;
             fac_react=exp(charge_i+x[iunk][inode_box]);
             resid = -(wt_surfsrc_1el[isten][inormal])*(4.0*PI/Temp_elec)*Rho_b[0]/(Rho_b[0]+fac_react);
             if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             resid_sum +=resid;

          if (resid_only_flag==FALSE){
               mat_val=-(wt_surfsrc_1el[isten][inormal])*(4.0*PI/Temp_elec)*(-Rho_b[0]*fac_react/POW_DOUBLE_INT((Rho_b[0]+fac_react),2));
             if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk,inode_box,mat_val);
          }

******END REACT SURF TYPE 2********/
/******REACT SURF TYPE 3*********
             if (WallType[iwall]==0){ 
                  pconst=0.7;
                  alpha_2=-24.0;
             }
             else {
                  pconst=0.5;
                  alpha_2=-9.0;
             }
             charge_i = 0.0;
             for (idim=0; idim<Ndim; idim++) charge_i += Charge_w_sum_els[loc_inode][idim];
             inormal=abs(Surf_normal[0][B2L_node[jnode_box]][0])-1;
             fac_react=exp(charge_i+x[iunk][inode_box]);
             fac_react2=exp(alpha_2-x[iunk][inode_box]);
             resid = -(wt_surfsrc_1el[isten][inormal])*(4.0*PI/Temp_elec)*((Rho_b[0]*pconst/(Rho_b[0]+fac_react))
                                         -(Rho_b[0]*(1.0-pconst)/(Rho_b[0]+fac_react2)));
             if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             resid_sum +=resid;

          if (resid_only_flag==FALSE){
               mat_val=-(wt_surfsrc_1el[isten][inormal])*(4.0*PI/Temp_elec)*((-Rho_b[0]*pconst*fac_react/POW_DOUBLE_INT((Rho_b[0]+fac_react),2))
                                            -(Rho_b[0]*(1.-pconst)*fac_react2/POW_DOUBLE_INT((Rho_b[0]+fac_react2),2)));
             if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk,inode_box,mat_val);
          }
******END REACT SURF TYPE 3********/
                 }
             }
/* End of new code for charge regulating surfaces*/


           }     /* end of loop over all local nodes in this fluid element */
         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */

       /* finally address nodes where there is a constant surface charge */
       if (Nodes_2_boundary_wall[Nlists_HW-1][inode_box] != -1){

           iwall     = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
           if (Type_bc_elec[WallType[iwall]] == CONST_CHARGE ){
/* CORRECT CODE FOR UNIFORM CONSTANT CHARGE BC*/
             charge_i = 0.0;
             for (idim=0; idim<Ndim; idim++) charge_i += Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];
             resid = -4.0*PI*charge_i/Temp_elec;
             if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             resid_sum +=resid;
/**use to comment out... */

/* CASE 3 FOR CHARGE REGULATING SURFACE.  Same as CASE 3, but including all params (not alpha_2=-infty)
 * ASSUME thate Charge_w_sum_els contains instead the parameter alpha. 
 * NOTE THAT TO DO THIS PROPERLY, alpha (or K_d or pK) will be read into a unique storate location.*
             if (WallType[iwall]==0){ 
                  pconst=0.7;
                  alpha_2=-24.0;
             }
             else {
                  pconst=0.7;
                  alpha_2=-24.0;
             }
             charge_i = 0.0;
             for (idim=0; idim<Ndim; idim++) charge_i += Charge_w_sum_els[loc_inode][idim];
             fac_react=exp(charge_i+x[iunk][inode_box]);
             fac_react2=exp(alpha_2-x[iunk][inode_box]);
             resid = -(4.0*PI/Temp_elec)*((Rho_b[0]*pconst/(Rho_b[0]+fac_react))
                                         -(Rho_b[0]*(1.0-pconst)/(Rho_b[0]+fac_react2)));
             if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             resid_sum +=resid;

          if (resid_only_flag==FALSE){
               mat_val=-(4.0*PI/Temp_elec)*((-Rho_b[0]*pconst*fac_react/POW_DOUBLE_INT((Rho_b[0]+fac_react),2))
                                            -(Rho_b[0]*(1.-pconst)*fac_react2/POW_DOUBLE_INT((Rho_b[0]+fac_react2),2)));
             if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[inode_box]+Solver_Unk[iunk]*Nnodes]+=mat_val;
             dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, iunk,inode_box,mat_val);
          }
* END OF FAKE CODE FOR CR SURFACES*/
/* CASE 4 FOR CHARGE REGULATING SURFACE.  Same as CASE 3, but using density rather than electrostatic potential
 * in the boundary conditions - first step to generalize beyond PB electrolytes...*
             if (WallType[iwall]==0){ 
                  pconst=0.7;
                  alpha_2=-24.;
             }
             else {
                  pconst=0.5;
                  alpha_2=-9.0;
             }
             charge_i = 0.0;
             for (idim=0; idim<Ndim; idim++) charge_i += Charge_w_sum_els[loc_inode][idim];
             junk=Phys2Unk_first[DENSITY];
             junk2react=Phys2Unk_first[DENSITY]+1;
             nshift=(int)((0.5+1.e-10)/Esize_x[0]);
       nshift=0; *for PB model*
             if (Surf_normal[Nlists_HW-1][B2L_node[inode_box]][0]>0) {jnode_box=inode_box+nshift;}
             else if (Surf_normal[Nlists_HW-1][B2L_node[inode_box]][0]<0) {jnode_box=inode_box-nshift;}
             else{
                 printf("Problem with Surf_normal in applying reacting boundary conditions\n");
                 printf("inode_box=%d iwall=%d Nlists_HW=%d surf_norm[Nlits_HW-1]=%d surf_normal[0]=%d\n",
                         inode_box,iwall,Nlists_HW,Surf_normal[Nlists_HW-1][B2L_node[inode_box]][0],Surf_normal[0][B2L_node[inode_box]][0]  );
                 exit(-1);
             }
             fac_react=1./(exp(charge_i)+x[junk][jnode_box]);
             fac_react2=1./(exp(alpha_2)+x[junk2react][jnode_box]);
             resid = -(4.0*PI/Temp_elec)*((pconst*x[junk][jnode_box]*fac_react) - ((1.0-pconst)*x[junk2react][jnode_box]*fac_react2));

             if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
             resid_sum +=resid;

          if (resid_only_flag==FALSE){
               mat_val=-(4.0*PI/Temp_elec)*((   pconst*(fac_react)*(1.0-x[junk][jnode_box]*fac_react)));
               if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[junk]*Nnodes]+=mat_val;
               dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, junk,jnode_box,mat_val);
               mat_val2=-(4.0*PI/Temp_elec)*( -((1.0-pconst)*(fac_react2)*(1.0-x[junk2react][jnode_box]*fac_react2)));
             if (Iwrite_files==FILES_DEBUG_MATRIX) Array_test[L2G_node[loc_inode]+Solver_Unk[iunk]*Nnodes][B2G_node[jnode_box]+Solver_Unk[junk2react]*Nnodes]+=mat_val2;
              dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode, junk2react,jnode_box,mat_val2);
          }

* END OF FAKE CODE FOR CR SURFACES*/

        }   /* check for charge b.c. on this wall */
     }      /* end of check for if this is a boundary node */

     if (resid_only_flag==INIT_GUESS_FLAG) resid_sum /= prefactor_sum;
   }
   return(resid_sum);
}
/************************************************************************/
/* load_poisson_bc: Load the Boundary condiditions associated with surfaces
                    of known charge.                                     */
double load_poisson_bc(int iunk,int loc_inode,int inode_box, int resid_only_flag)
{
  int iwall,idim;
  double charge_i,resid=0.0;

     if (Nodes_2_boundary_wall[Nlists_HW-1][inode_box] != -1){

         iwall     = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
         if (Type_bc_elec[WallType[iwall]] == CONST_CHARGE ){

          charge_i = 0.0;
          for (idim=0; idim<Ndim; idim++)
             charge_i += Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];

          resid = -4.0*PI*charge_i/Temp_elec;

          if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

        }   /* check for charge b.c. on this wall */
     }      /* end of check for if this is a boundary node */
  return(resid);
}
/************************************************************************/
