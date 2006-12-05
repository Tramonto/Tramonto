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
 *  FILE: dft_fill_pde.c
 *
 *  This file contains the fills of partial differential equations
 *  for electrostatics and transport.
 */

#include "dft_globals_const.h"
#include "rf_allo.h"

void set_fem_1el_weights(double **, double **, int ***);
void set_fem_weights(double **, double **);

/****************************************************************************/
void set_fem_weights(double **wt_laplace_ptr, double **wt_source_ptr)
/* This routine sets the fem weights for the laplace and source terms */

{
  int num_nodes;
  double dy_x, dx_y, evol, dxy_z, dxz_y, dyz_x;

  num_nodes = POW_INT(3, Ndim);

  if (*wt_laplace_ptr == NULL) {
    *wt_laplace_ptr = (double *) array_alloc(1, num_nodes, sizeof(double));
    *wt_source_ptr  = (double *) array_alloc(1, num_nodes, sizeof(double)); 
  }

  if (Ndim == 1) {
    /*Finite element weights for Laplace terms on rectangular grid */
    (*wt_laplace_ptr)[0] = -1.0/Esize_x[0];
    (*wt_laplace_ptr)[1] =  2.0/Esize_x[0];
    (*wt_laplace_ptr)[2] = -1.0/Esize_x[0];

    /*Finite element weights for source terms on rectangular grid */
    (*wt_source_ptr)[0]  = Esize_x[0]/6.0;
    (*wt_source_ptr)[1]  = Esize_x[0]*2.0/3.0;
    (*wt_source_ptr)[2]  = Esize_x[0]/6.0;
  }

  else  if (Ndim == 2) {
    dx_y = Esize_x[0]/Esize_x[1];
    dy_x = Esize_x[1]/Esize_x[0];

    /*Finite element weights for Laplace terms on rectangular grid */
    /* corners */
    (*wt_laplace_ptr)[0] =  (*wt_laplace_ptr)[2] = (*wt_laplace_ptr)[6] =
      (*wt_laplace_ptr)[8] = - (dy_x + dx_y)/6.0;
    /* top and bottom */
    (*wt_laplace_ptr)[1] = (*wt_laplace_ptr)[7] = (dy_x - 2.0*dx_y)/3.0;
    /* left and right */
    (*wt_laplace_ptr)[3] = (*wt_laplace_ptr)[5] = (-2.0*dy_x + dx_y)/3.0;
    /* center */
    (*wt_laplace_ptr)[4] = 4.0*(dy_x + dx_y)/3.0;;

    /*Finite element weights for source terms on rectangular grid */
    evol = Esize_x[0]*Esize_x[1];
    /* corners */
    (*wt_source_ptr)[0] =  (*wt_source_ptr)[2] = (*wt_source_ptr)[6] =
      (*wt_source_ptr)[8] = evol / 36.0;
    /* top and bottom, left and right */
    (*wt_source_ptr)[1] = (*wt_source_ptr)[7] = (*wt_source_ptr)[3] = 
      (*wt_source_ptr)[5] = evol / 9.0;
    /* center */
    (*wt_source_ptr)[4] = 4.0 * evol / 9.0;

  }

  else {

    /*Finite element weights for Laplace terms on rectangular grid */
    dxy_z = Esize_x[0]*Esize_x[1]/Esize_x[2];
    dyz_x = Esize_x[1]*Esize_x[2]/Esize_x[0];
    dxz_y = Esize_x[2]*Esize_x[0]/Esize_x[1];
    /* center */
    (*wt_laplace_ptr)[13] = 8.0*(dyz_x + dxz_y + dxy_z)/9.0;
    /* +- ix */
    (*wt_laplace_ptr)[12] = (*wt_laplace_ptr)[14] = 2.0*(-2*dyz_x + dxz_y + dxy_z)/9.0;
    /* +- iy */
    (*wt_laplace_ptr)[10] = (*wt_laplace_ptr)[16] = 2.0*(dyz_x - 2.0*dxz_y +dxy_z)/9.0;
    /* +- iz */
    (*wt_laplace_ptr)[4] = (*wt_laplace_ptr)[22] = 2.0*(dyz_x + dxz_y -2.0*dxy_z)/9.0;
    /* +- ix and +- iy */
    (*wt_laplace_ptr)[9] = (*wt_laplace_ptr)[11] = (*wt_laplace_ptr)[15]
        =  (*wt_laplace_ptr)[17] = (-dyz_x - dxz_y + 0.5*dxy_z)/9.0;
    /* +- ix and +- iz */
    (*wt_laplace_ptr)[3] = (*wt_laplace_ptr)[5] = (*wt_laplace_ptr)[21]
        =  (*wt_laplace_ptr)[23] = (-dyz_x + 0.5*dxz_y - dxy_z)/9.0;
    /* +- iy and +- iz */
    (*wt_laplace_ptr)[1] = (*wt_laplace_ptr)[7] = (*wt_laplace_ptr)[19]
        =  (*wt_laplace_ptr)[25] = (0.5*dyz_x - dxz_y - dxy_z)/9.0;
    /* (corners) +- ix and +- iy and +- iz */
    (*wt_laplace_ptr)[ 0] = (*wt_laplace_ptr)[ 2] = (*wt_laplace_ptr)[6]
        =  (*wt_laplace_ptr)[ 8] = (*wt_laplace_ptr)[18] = (*wt_laplace_ptr)[20]
        =  (*wt_laplace_ptr)[24] = (*wt_laplace_ptr)[26] = -(dyz_x + dxz_y + dxy_z)/36.0;

    /*Finite element weights for source terms on rectangular grid */
    evol = Esize_x[0]*Esize_x[1]*Esize_x[2];

    /* center */
    (*wt_source_ptr)[13] = 64.0 * evol/216.0;

    /* One away from center */
    (*wt_laplace_ptr)[12] = (*wt_laplace_ptr)[14] = 
      (*wt_laplace_ptr)[10] = (*wt_laplace_ptr)[16] = 
      (*wt_laplace_ptr)[ 4] = (*wt_laplace_ptr)[22] = 16.0 * evol/216.0;

    /* two away from center */
    (*wt_source_ptr)[ 9] = (*wt_source_ptr)[11] = (*wt_source_ptr)[15]
        =  (*wt_source_ptr)[17] 
        =  (*wt_source_ptr)[ 3] = (*wt_source_ptr)[ 5] = (*wt_source_ptr)[21]
        =  (*wt_source_ptr)[23]
        =  (*wt_source_ptr)[ 1] = (*wt_source_ptr)[ 7] = (*wt_source_ptr)[19]
        =  (*wt_source_ptr)[25] = 4.0 * evol/216.0;

    /* (corners) +- ix and +- iy and +- iz */
    (*wt_source_ptr)[ 0] = (*wt_source_ptr)[ 2] = (*wt_source_ptr)[ 6]
        =  (*wt_source_ptr)[ 8] = (*wt_source_ptr)[18] = (*wt_source_ptr)[20]
        =  (*wt_source_ptr)[24] =  (*wt_source_ptr)[26] = evol/216.0;
  }

}
/****************************************************************************/
void set_fem_1el_weights(double **wt_lp_1el_ptr, double **wt_s_1el_ptr,
                         int ***elem_permute)
/* This routine sets the fem weights for the laplace and source terms */

{
  int iln;
  double dy_x, dx_y, evol, dxy_z, dxz_y, dyz_x;

  if (*wt_lp_1el_ptr == NULL) {
    *wt_lp_1el_ptr = (double *) array_alloc(1, Nnodes_per_el_V,sizeof(double));
    *wt_s_1el_ptr  = (double *) array_alloc(1, Nnodes_per_el_V,sizeof(double));
    *elem_permute  = (int **) array_alloc(2, Nnodes_per_el_V,Nnodes_per_el_V, sizeof(int));
  }

  if (Ndim == 1) {
    /*Finite element weights for Laplace terms on rectangular grid */
    (*wt_lp_1el_ptr)[0] =  1.0/Esize_x[0];
    (*wt_lp_1el_ptr)[1] = -1.0/Esize_x[0];

    /*Finite element weights for source terms on rectangular grid */
    (*wt_s_1el_ptr)[0]  = Esize_x[0]/3.0;
    (*wt_s_1el_ptr)[1]  = Esize_x[0]/6.0;

    (*elem_permute)[0][0]=0;
    (*elem_permute)[1][0]=1;
    (*elem_permute)[0][1]=1;
    (*elem_permute)[1][1]=0;

  }

  else  if (Ndim == 2) {
    dx_y = Esize_x[0]/Esize_x[1];
    dy_x = Esize_x[1]/Esize_x[0];

    /*Finite element weights for Laplace terms on rectangular grid */
    (*wt_lp_1el_ptr)[0] = (dy_x + dx_y) /3.0;
    (*wt_lp_1el_ptr)[1] = (-dy_x + dx_y/2.0) /3.0;
    (*wt_lp_1el_ptr)[2] = ( dy_x/2.0 - dx_y) /3.0;
    (*wt_lp_1el_ptr)[3] = -(dy_x + dx_y) /6.0;

    /*Finite element weights for source terms on rectangular grid */
    evol = Esize_x[0]*Esize_x[1];
    (*wt_s_1el_ptr)[0] = evol/9.0;
    (*wt_s_1el_ptr)[1] = evol/18.0;
    (*wt_s_1el_ptr)[2] = evol/18.0;
    (*wt_s_1el_ptr)[3] = evol/36.0;

    for (iln=0; iln<Nnodes_per_el_V; iln++) {
       (*elem_permute)[iln][0] = iln;
       (*elem_permute)[iln][1] = (13-iln)%4;
       (*elem_permute)[iln][2] = (iln+6)%4;
       (*elem_permute)[iln][3] = 3 - iln;
    }
  }

  else {

    /*Finite element weights for Laplace terms on rectangular grid */
    dxy_z = Esize_x[0]*Esize_x[1]/Esize_x[2];
    dyz_x = Esize_x[1]*Esize_x[2]/Esize_x[0];
    dxz_y = Esize_x[2]*Esize_x[0]/Esize_x[1];

    (*wt_lp_1el_ptr)[0] = (dyz_x + dxz_y + dxy_z) /9.0;
    (*wt_lp_1el_ptr)[1] = (-dyz_x + dxz_y/2.0 + dxy_z/2.0) /9.0;
    (*wt_lp_1el_ptr)[2] = (dyz_x/2.0 - dxz_y + dxy_z/2.0) /9.0;
    (*wt_lp_1el_ptr)[3] = (-dyz_x/2.0 - dxz_y/2.0 + dxy_z/4.0) /9.0;
    (*wt_lp_1el_ptr)[4] = (dyz_x/2.0 + dxz_y/2.0 - dxy_z) /9.0;
    (*wt_lp_1el_ptr)[5] = (-dyz_x/2.0 + dxz_y/4.0 - dxy_z/2.0) /9.0;
    (*wt_lp_1el_ptr)[6] = (dyz_x/4.0 - dxz_y/2.0 - dxy_z/2.0) /9.0;
    (*wt_lp_1el_ptr)[7] = -(dyz_x + dxz_y + dxy_z) /36.0;
    
    /*Finite element weights for source terms on rectangular grid */
    evol = Esize_x[0]*Esize_x[1]*Esize_x[2];
    (*wt_s_1el_ptr)[0] = evol/27.0;
    (*wt_s_1el_ptr)[1] = evol/54.0;
    (*wt_s_1el_ptr)[2] = evol/54.0;
    (*wt_s_1el_ptr)[3] = evol/108.0;
    (*wt_s_1el_ptr)[4] = evol/54.0;
    (*wt_s_1el_ptr)[5] = evol/108.0;
    (*wt_s_1el_ptr)[6] = evol/108.0;
    (*wt_s_1el_ptr)[7] = evol/216.0;

    for (iln=0; iln<4; iln++) {
       (*elem_permute)[iln][0] = (*elem_permute)[iln+4][0+4] = iln;
       (*elem_permute)[iln][1] = (*elem_permute)[iln+4][1+4] = (13-iln)%4;
       (*elem_permute)[iln][2] = (*elem_permute)[iln+4][2+4] = (iln+6)%4;
       (*elem_permute)[iln][3] = (*elem_permute)[iln+4][3+4] = 3-iln; 
    }
    for (iln=4; iln<8; iln++) {
       (*elem_permute)[iln][0] = (*elem_permute)[iln-4][0+4] = iln;
       (*elem_permute)[iln][1] = (*elem_permute)[iln-4][1+4] = (13-iln)%4 + 4;
       (*elem_permute)[iln][2] = (*elem_permute)[iln-4][2+4] = (iln+6)%4 + 4;
       (*elem_permute)[iln][3] = (*elem_permute)[iln-4][3+4] = 7- (iln-4); 
    }
  }
}
/****************************************************************************/
double load_polarize_poissons_eqn(int iunk, int loc_inode, int inode_box, int *ijk_box, double **x)
{

  int iwall,  isten, icomp,ilist,idim,junk,junkP;
  int iln, jln, elem, offset[3], el_box;
  int nodes_volm_el, nodes_surf_el, junk2[3];
  int in_wall,numEntries, nodeIndices[2];
  double values[2];
  double pol_wall,wt;
  double resid,resid_sum=0.0,mat_val;

  /* static variables keep their value for every time the function is called*/
  static double *wt_lp_1el, *wt_s_1el;
  static int   **elem_permute, off_ref[2][2] = { {0,1}, {-1,0}};
  double rho_0, rho_1, psi_0, psi_1, tmp;
  int jnode_box,junk_psi,junk_rho;
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
          resid = x[iunk][inode_box] - Elec_param_w[iwall];
          resid_sum+=resid;
          mat_val = 1.0;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
   }
   else {
       for (iln=0; iln< Nnodes_per_el_V; iln++) {

         /*elem = node_to_elem_v2(node_box_to_node(inode_box), iln);*/
         elem = node_to_elem(node_box_to_node(inode_box), iln,reflect_flag);
         el_box = el_to_el_box(elem);

         if (elem >= 0 ) {

           if (Vol_charge_flag) {
               resid = -Charge_vol_els[el_box]*4.0*PI/(Temp_elec * Nnodes_per_el_V);
               resid_sum+=resid;
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
           }

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
             junkP = Phys2Unk_first[POISSON];
             if (jnode_box >= 0){  /* new flag for boundaries */
                   in_wall=TRUE;
                   for (icomp=0; icomp<Ncomp; icomp++){
                      if (!Zero_density_TF[jnode_box][icomp])in_wall=FALSE;
                   }
                   if (in_wall){
                      /* pol_wall is a constant wall polarization.  set equal
                         to zero, but leave code as place holder. */
                      pol_wall=(KAPPA_H2O-1.0);
                      pol_wall=0.0;
                      resid = (1.0+pol_wall)*wt_lp_1el[isten]*x[junkP][jnode_box];
                      mat_val = (1.0+pol_wall)*wt_lp_1el[isten];
                   }
                   else{
                      resid = wt_lp_1el[isten]*x[junkP][jnode_box];
                      mat_val = wt_lp_1el[isten];
                   }
                   resid_sum+=resid;
                   dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                               junkP,jnode_box,mat_val);
             }

             /* 
              * add in source term (electroneutrality sum) 
              */
             for (junk=Phys2Unk_first[DENSITY]; junk<Phys2Unk_last[DENSITY]; junk++){
               if (Type_poly==WTC) icomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
               else                icomp=junk-Phys2Unk_first[DENSITY];
 
               if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
               else ilist = icomp;
               if (elem !=-2 && (Wall_elems[ilist][el_box] == -1 ||
                   Lsemiperm[WallType[Wall_elems[ilist][el_box]]][icomp]) ){

                  if (Charge_f[icomp] != 0.0) {
                    resid = -wt_s_1el[isten]*KAPPA_H2O*Charge_f[icomp]*x[junk][jnode_box]*
                                                                         4.0*PI/Temp_elec;
                    mat_val = -wt_s_1el[isten]*KAPPA_H2O*Charge_f[icomp]*4.0*PI/Temp_elec;
                    dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                junk,jnode_box,mat_val);
                    
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

           junk_psi = Phys2Unk_first[POISSON];
           for (junk_rho=Phys2Unk_first[DENSITY]; junk_rho<Phys2Unk_last[DENSITY]; junk_rho++){
               icomp=junk-Phys2Unk_first[DENSITY];
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
                    psi_0  = x[junk_psi][nodeIndices[0]];
                    psi_1  = x[junk_psi][nodeIndices[1]];

                    tmp = (rho_0 + rho_1);

                    resid = wt*Pol[icomp]*(psi_0 - psi_1)*tmp;
                    resid_sum+=resid;
                    dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

                    numEntries=2;
                    values[0]=values[1]=wt*Pol[icomp]*(psi_0 - psi_1);
                    dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                      junk_rho, nodeIndices,values,numEntries);
                    values[0]=wt*Pol[icomp]*tmp; values[1]=-values[0];
                    dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                      junk_psi, nodeIndices,values,numEntries);
                  }
               }
             } /*end of test for polarizeable fluid species */
           }   /*end of loop over components for adding in polarization term */


         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */
   }

}
/****************************************************************************/
double load_poissons_eqn(int iunk, int loc_inode, int inode_box, int *ijk_box, double **x)
{

  int iwall,  isten, icomp,ilist,idim;
  int iln, jln, elem, offset[3], el_box;
  int nodes_volm_el, nodes_surf_el, junk2[3];
  /* static variables keep their value for every time the function is called*/
  static double *wt_lp_1el, *wt_s_1el;
  static int   **elem_permute;
  int off_ref[2][2]; /*= { {0,1}, {-1,0}};*/
  int jnode_box,junk,junkP;
  int reflect_flag[3];
  double resid,mat_val,resid_sum=0.0,mat_val_sum;
  off_ref[0][0]=0;
  off_ref[0][1]=1;
  off_ref[1][0]=-1;
  off_ref[1][1]=0;

   mat_val_sum=0.0;
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
          resid = x[iunk][inode_box] - Elec_param_w[iwall];
          mat_val=1.0;
          resid_sum+=resid;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
/*printf("1 load matrix value - diagaonal: iunk=%d loc_inode=%d iunk=%d inode_box=%d mat_val=%11.8f\n",
        iunk,loc_inode,iunk,inode_box,mat_val);*/
   }
   else {
       for (iln=0; iln< Nnodes_per_el_V; iln++) {

         /*elem = node_to_elem_v2(node_box_to_node(inode_box), iln);*/
         elem = node_to_elem(node_box_to_node(inode_box), iln,reflect_flag);
         el_box = el_to_el_box(elem);

         if (elem >= 0 ) {
             if (Vol_charge_flag){ 
               resid = -Charge_vol_els[el_box]*4.0*PI
                              / (Temp_elec * Nnodes_per_el_V);
               resid_sum += resid;
               dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
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
                   resid = Dielec[el_box]*wt_lp_1el[isten]*x[junkP][jnode_box];
                   mat_val = Dielec[el_box]*wt_lp_1el[isten];
                   dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                               junkP,jnode_box,mat_val);
                   mat_val_sum +=mat_val;
/*printf("2 load matrix value : iunk=%d loc_inode=%d junkP=%d jnode_box=%d mat_val=%11.8f mat_val_sum=%11.8f\n",
        iunk,loc_inode,junkP,jnode_box,mat_val,mat_val_sum);*/
                   resid_sum+=resid;
             }

             /* 
              * add in source term (electroneutrality sum) 
              */
             for (junk=Phys2Unk_first[DENSITY]; junk<Phys2Unk_last[DENSITY]; junk++){
               if (Type_poly==WTC) icomp=Unk2Comp[junk-Phys2Unk_first[DENSITY]];
               else                icomp=junk-Phys2Unk_first[DENSITY];

               if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
               else ilist = icomp;

               if (elem !=-2 && (Wall_elems[ilist][el_box] == -1 ||
                   Lsemiperm[WallType[Wall_elems[ilist][el_box]]][icomp]) ){

                  if (Charge_f[icomp] != 0.0) {
                       resid = -wt_s_1el[isten]*Charge_f[icomp]*x[junk][jnode_box]*4.0*PI/Temp_elec;
                       mat_val = -wt_s_1el[isten]*Charge_f[icomp]* 4.0*PI/Temp_elec;
                       dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                       dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                  junk,jnode_box,mat_val);
                   mat_val_sum +=mat_val;
/*printf("3 load matrix value : iunk=%d loc_inode=%d junk=%d jnode_box=%d mat_val=%11.8f  mat_val_sum=%11.8f\n",
        iunk,loc_inode,junk,jnode_box,mat_val,mat_val_sum);*/
                       resid_sum+=resid;
                    
                  } /* End if (Charge_f[icomp] != 0.0) */

               }
             }   /* end of icomp (source term) loop */


           }     /* end of loop over all local nodes in this fluid element */
         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */
   }
   return(resid_sum);
}
/************************************************************************/
/* load_poisson_bc: Load the Boundary condiditions associated with surfaces
                    of known charge.                                     */
double load_poisson_bc(int iunk,int loc_inode,int inode_box )
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
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

        }   /* check for charge b.c. on this wall */
     }      /* end of check for if this is a boundary node */
  return(resid);
}
/************************************************************************/
/* basis_fn_calc: calcs phi and grad_phi for regular hex element */
void basis_fn_calc(double **phi, double ***grad_phi, double *evol)
{
   int iln, igp;
   double gp[2];

   gp[0] = (1.0 - 5.773502691896258e-01) / 2.0;
   gp[1] = (1.0 + 5.773502691896258e-01) / 2.0;

   if (Ndim == 3) {

      *evol = Esize_x[0] * Esize_x[1] * Esize_x[2];
      
      for (iln=0; iln<8; iln++) {
        for (igp=0; igp<8; igp++) {
          phi[iln][igp] = gp[(iln+igp)%2] * gp[((iln/2) + (igp/2))%2]
                                          * gp[((iln/4) + (igp/4))%2];
          grad_phi[iln][igp][0] = gp[((iln/2) + (igp/2))%2]
                                * gp[((iln/4) + (igp/4))%2] / Esize_x[0];
          if (((iln+igp)%2)==0) grad_phi[iln][igp][0] *= -1.0;
          grad_phi[iln][igp][1] = gp[(iln+igp)%2]  
                                * gp[((iln/4) + (igp/4))%2] / Esize_x[1];
          if (((iln/2) + (igp/2))%2==0) grad_phi[iln][igp][1] *= -1.0;
          grad_phi[iln][igp][2] = gp[(iln+igp)%2]  
                                * gp[((iln/2) + (igp/2))%2] / Esize_x[2];
          if (((iln/4) + (igp/4))%2==0) grad_phi[iln][igp][2] *= -1.0;
        }
      }
   }
   else {
     if (Proc==0) printf("ERROR basis_fn_calc: Only for 3D\n");
     exit(-1);
   }
}
/************************************************************************/
double load_nonlinear_transport_eqn(int iunk, int loc_inode, int inode_box,
                                    int *ijk_box, double **x)
{

  int ilist=0, icomp, idim, ijk[3],inode;
  int iln, jln, elem, offset[3];
  int nodes_volm_el, nodes_surf_el, junk2[3];
  int off_ref[2][2] = { {0,1}, {-1,0}};
  int j_box_mu[8], j_box_rho[8], loc_j_mu[8], loc_j_rho[8], jnode_box,flag,count_flag;
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
          mat_val = 1.0;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
   }
   /* check if you are in well-mixed bulk region */
   else if (ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
          resid = x[iunk][inode_box] - Betamu_LBB[icomp];
          resid_sum+=resid;
          mat_val = 1.0;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
   }
   else if (Size_x[Grad_dim]-ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
          resid = x[iunk][inode_box] - Betamu_RTF[icomp];
          resid_sum+=resid;
          mat_val = 1.0;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
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
           junk_rho = Phys2Unk_first[DENSITY];
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
 
                numEntries=2;
                values[0]=values[1]=wt*(mu_0 - mu_1)*(2.0*area_0 + area_1);
                dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                  junk_rho, nodeIndices,values,numEntries);
                values[0] = wt*tmp;  values[1]=-values[0]; 
                dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                  junk_mu, nodeIndices,values,numEntries);

                /* add in a bulk convection term */
                if (iln == 0) tmp = (2.0*area_0 + area_1)/6.0;
                else          tmp = (area_0 + 2.0*area_1)/6.0;
                resid = Velocity * (rho_1 - rho_0) * tmp;
                resid_sum+=resid;
                dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                values[0] = -Velocity*tmp;  values[1]=-values[0];
                dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                                  junk_rho, nodeIndices,values,numEntries);
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
                                 int *ijk_box,double **x)
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
          mat_val = 1.0;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
   }
    if (ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
          resid = x[iunk][inode_box] - Betamu_LBB[icomp];
          resid_sum+=resid;
          mat_val = 1.0;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
   }
   else if (Size_x[Grad_dim]-ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
          resid = x[iunk][inode_box] - Betamu_RTF[icomp];
          resid_sum+=resid;
          mat_val = 1.0;
          dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
          dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                      iunk,inode_box,mat_val);
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
                   mat_val = wt_lp_1el[isten];
                   dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
                   dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,
                                                                iunk,jnode_box,mat_val);
             }

           }     /* end of loop over all local nodes in this fluid element */
         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */
   }
   return(resid_sum);
}
/****************************************************************************/
