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
void load_polarize_poissons_eqn(int i_box, int inode_box, int loc_i, int *ijk_box,
                       double *mat_row,
                       double *resid, double *x, int *bindx_tmp, int fill_flag)
{

  int iwall,  isten, icomp,ilist,idim;
  int iln, jln, elem, offset[3], loc_j,el_box;
  int nodes_volm_el, nodes_surf_el, junk2[3];
  int in_wall;
  double pol_wall,wt;

  /* static variables keep their value for every time the function is called*/
  static double *wt_lp_1el, *wt_s_1el;
  static int   **elem_permute, off_ref[2][2] = { {0,1}, {-1,0}};
  int j_box_psi[8], j_box_rho[8], loc_j_psi[8], loc_j_rho[8];
  double rho_0, rho_1, psi_0, psi_1, tmp;
  int j_box, jnode_box;
  int reflect_flag[3];

  /* First time through, load weights appropriate for this Ndim */

  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;
  
  if (MSR_PREPROCESS){
    set_fem_1el_weights(&wt_lp_1el, &wt_s_1el, &elem_permute);
  }

/* if we are at a boundary node and the boundary condition is constant
 * potential, then solve the equation psi = constant ..... in all other
 * cases loop over the elements near this node and fill the laplace and
 * source terms.  Note that the surface charge boundary condition is
 * handled in the main fill routine. */

   iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
   if (iwall != -1 && Type_bc_elec[WallType[iwall]] == CONST_POTENTIAL) {
       if (fill_flag != MSR_PREPROCESS){
          resid[loc_i] = x[loc_i] - Elec_param_w[iwall];
          mat_row[loc_i] = 1.0;
       }
       else bindx_tmp[i_box]=TRUE;
   }
   else {
       for (iln=0; iln< Nnodes_per_el_V; iln++) {

         /*elem = node_to_elem_v2(node_box_to_node(inode_box), iln);*/
         elem = node_to_elem(node_box_to_node(inode_box), iln,reflect_flag);
         el_box = el_to_el_box(elem);

         if (elem >= 0 ) {

           if (fill_flag != MSR_PREPROCESS){
             if (Vol_charge_flag) 
               resid[loc_i] -= Charge_vol_els[el_box]*4.0*PI
                              / (Temp_elec * Nnodes_per_el_V);
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
             j_box = jnode_box*Nunk_per_node + Ncomp+Nrho_bar;

             /*
              * add in Laplace term
              */
             if (jnode_box >= 0){  /* new flag for boundaries */
                if (fill_flag != MSR_PREPROCESS){
                   loc_j = B2L_unknowns[j_box];
                   in_wall=TRUE;
                   for (icomp=0; icomp<Ncomp; icomp++){
                      if (!Zero_density_TF[jnode_box][icomp])in_wall=FALSE;
                   }
                   if (in_wall){
                      /* pol_wall is a constant wall polarization.  set equal
                         to zero, but leave code as place holder. */
                      pol_wall=(KAPPA_H2O-1.0);
                      pol_wall=0.0;
                      resid[loc_i] += (1.0+pol_wall)*wt_lp_1el[isten]*x[loc_j];
                      mat_row[loc_j] += (1.0+pol_wall)*wt_lp_1el[isten];
                   }
                   else{
                      resid[loc_i] += wt_lp_1el[isten]*x[loc_j];
                      mat_row[loc_j] += wt_lp_1el[isten];
                   }
                }
                else bindx_tmp[j_box]=TRUE;
             }

             /* 
              * add in source term (electroneutrality sum) 
              */
             for (icomp=0; icomp<Ncomp; icomp++) {

               if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
               else ilist = icomp;
               if (elem !=-2 && (Wall_elems[ilist][el_box] == -1 ||
                   Lsemiperm[WallType[Wall_elems[ilist][el_box]]][icomp]) ){

                  if (Charge_f[icomp] != 0.0) {
                    j_box = jnode_box*Nunk_per_node + icomp;
                    if (fill_flag != MSR_PREPROCESS){
                       loc_j = B2L_unknowns[j_box];
                       resid[loc_i] -= wt_s_1el[isten]*KAPPA_H2O*Charge_f[icomp]*x[loc_j]
                                                         *4.0*PI/Temp_elec;
                       mat_row[loc_j] -= wt_s_1el[isten]*KAPPA_H2O*Charge_f[icomp]*
                                                          4.0*PI/Temp_elec;
                    }
                    else bindx_tmp[j_box] = TRUE;
                    
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

           for (icomp=0; icomp<Ncomp; icomp ++){
               if (!Zero_density_TF[inode_box][icomp] && Lpolarize[icomp]){

               for (jln=0; jln< Nnodes_per_el_V; jln++) {
                  for (idim=0; idim<Ndim; idim++){
                      nodes_surf_el = POW_INT(2,idim);
                      nodes_volm_el = POW_INT(2,idim+1);
                      offset[idim] =
                           off_ref[((nodes_volm_el+iln)%nodes_volm_el)/nodes_surf_el]
                                  [((nodes_volm_el+jln)%nodes_volm_el)/nodes_surf_el];
                  }

                  jnode_box = offset_to_node_box(ijk_box, offset, junk2);

                  j_box_psi[jln]  = jnode_box*Nunk_per_node + Ncomp;
                  j_box_rho[jln] = jnode_box*Nunk_per_node + icomp;
                  if (fill_flag != MSR_PREPROCESS){
                    loc_j_psi[jln]  = B2L_unknowns[j_box_psi[jln]];
                    loc_j_rho[jln] = B2L_unknowns[j_box_rho[jln]];
                  }
               }

               in_wall=TRUE;
               for (icomp=0; icomp<Ncomp; icomp++){
                   if (!Zero_density_TF[jnode_box][icomp])in_wall=FALSE;
               }

               if (!in_wall){
               if (fill_flag != MSR_PREPROCESS){

                  if (Ndim==1) {

                    rho_0 = x[loc_j_rho[0]];
                    rho_1 = x[loc_j_rho[1]];
                    psi_0  = x[loc_j_psi[0]];
                    psi_1  = x[loc_j_psi[1]];
                    tmp = (rho_0 + rho_1);

                    resid[loc_i] += wt*Pol[icomp]*(psi_0 - psi_1)*tmp;

                    mat_row[loc_j_rho[0]] += wt*Pol[icomp]*(psi_0 - psi_1);
                    mat_row[loc_j_rho[1]] += wt*Pol[icomp]*(psi_0 - psi_1);
                    mat_row[loc_j_psi[0]] += wt*Pol[icomp]*tmp;
                    mat_row[loc_j_psi[1]]  -= wt*Pol[icomp]*tmp;

                  }
               }
               else {
                 for (jln=0; jln< Nnodes_per_el_V; jln++) {
                   bindx_tmp[j_box_psi[jln]] =TRUE;
                   bindx_tmp[j_box_rho[jln]]=TRUE;
                 }
               }
               }
             } /*end of test for polarizeable fluid species */
           }   /*end of loop over components for adding in polarization term */


         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */
   }

}
/****************************************************************************/
void load_poissons_eqn(int i_box, int inode_box, int loc_i, int *ijk_box,
                       double *mat_row,
                       double *resid, double *x, int *bindx_tmp, int fill_flag)
{

  int iwall,  isten, icomp,ilist,idim;
  int iln, jln, elem, offset[3], loc_j,el_box;
  int nodes_volm_el, nodes_surf_el, junk2[3];
  /* static variables keep their value for every time the function is called*/
  static double *wt_lp_1el, *wt_s_1el;
  static int   **elem_permute, off_ref[2][2] = { {0,1}, {-1,0}};
  int j_box, jnode_box;
  int reflect_flag[3];

  /* First time through, load weights appropriate for this Ndim */

  for (idim=0; idim<Ndim; idim++) reflect_flag[idim]=FALSE;
  
  if (MSR_PREPROCESS){
    set_fem_1el_weights(&wt_lp_1el, &wt_s_1el, &elem_permute);
  }

/* if we are at a boundary node and the boundary condition is constant
 * potential, then solve the equation psi = constant ..... in all other
 * cases loop over the elements near this node and fill the laplace and
 * source terms.  Note that the surface charge boundary condition is
 * handled in the main fill routine. */

   iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
   if (iwall != -1 && Type_bc_elec[WallType[iwall]] == CONST_POTENTIAL) {
       if (fill_flag != MSR_PREPROCESS){
          resid[loc_i] = x[loc_i] - Elec_param_w[iwall];
          mat_row[loc_i] = 1.0;
       }
       else{ 
          bindx_tmp[i_box]=TRUE;
       }
   }
   else {
       for (iln=0; iln< Nnodes_per_el_V; iln++) {

         /*elem = node_to_elem_v2(node_box_to_node(inode_box), iln);*/
         elem = node_to_elem(node_box_to_node(inode_box), iln,reflect_flag);
         el_box = el_to_el_box(elem);

         if (elem >= 0 ) {
           if (fill_flag != MSR_PREPROCESS){
             if (Vol_charge_flag) 
               resid[loc_i] -= Charge_vol_els[el_box]*4.0*PI
                              / (Temp_elec * Nnodes_per_el_V);
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
             if (Type_poly==-1) j_box = jnode_box*Nunk_per_node + Ncomp+Nrho_bar;
             else               j_box = jnode_box*Nunk_per_node + 2*Ncomp+Ngeqn_tot;
             /*
              * add in Laplace term
              */
             if (jnode_box >= 0){  /* new flag for boundaries */
                if (fill_flag != MSR_PREPROCESS){
                   loc_j = B2L_unknowns[j_box];
                   resid[loc_i] += Dielec[el_box]*wt_lp_1el[isten]*x[loc_j];
                   mat_row[loc_j] += Dielec[el_box]*wt_lp_1el[isten];
                }
                else bindx_tmp[j_box]=TRUE;
             }

             /* 
              * add in source term (electroneutrality sum) 
              */
             for (icomp=0; icomp<Ncomp; icomp++) {

               if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
               else ilist = icomp;
               if (elem !=-2 && (Wall_elems[ilist][el_box] == -1 ||
                   Lsemiperm[WallType[Wall_elems[ilist][el_box]]][icomp]) ){

                  if (Charge_f[icomp] != 0.0) {
                    if (Type_poly==-1) j_box = jnode_box*Nunk_per_node + icomp;
                    else               j_box = jnode_box*Nunk_per_node + Ncomp+icomp;
                    if (fill_flag != MSR_PREPROCESS){
                       loc_j = B2L_unknowns[j_box];
                       resid[loc_i] -= wt_s_1el[isten]*Charge_f[icomp]*x[loc_j]
                                                         *4.0*PI/Temp_elec;
                       mat_row[loc_j] -= wt_s_1el[isten]*Charge_f[icomp]*
                                                          4.0*PI/Temp_elec;
                    }
                    else bindx_tmp[j_box] = TRUE;
                    
                  } /* End if (Charge_f[icomp] != 0.0) */

               }
             }   /* end of icomp (source term) loop */


           }     /* end of loop over all local nodes in this fluid element */
         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */
   }

}
/************************************************************************/
/* load_poisson_bc: Load the Boundary condiditions associated with surfaces
                    of known charge.                                     
void load_poisson_bc_old(double *resid)
{
  int loc_inode, inode_box, loc_i, iwall,idim;
  double charge_i,fac;

  if (Type_coul==POLARIZE) fac=KAPPA_H2O;
  else              fac=1.0;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     
     if (Nodes_2_boundary_wall[Nlists_HW-1][inode_box] != -1){

         iwall     = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
         if (Type_bc_elec[WallType[iwall]] == CONST_CHARGE ){

          loc_i = Aztec.update_index[Ncomp+ Nrho_bar+ Nunk_per_node * loc_inode];

          charge_i = 0.0;
          for (idim=0; idim<Ndim; idim++)
             charge_i += Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];
          resid[loc_i] -= 4.0*PI*charge_i*fac/Temp_elec;

        }   
     }
  }
  return;
}*/
/************************************************************************/
/************************************************************************/
/* load_poisson_bc: Load the Boundary condiditions associated with surfaces
                    of known charge.                                     */
void load_poisson_bc(double *resid,int inode_box,int loc_inode,int loc_i)
{
  int iwall,idim;
  double charge_i;

     if (Nodes_2_boundary_wall[Nlists_HW-1][inode_box] != -1){

         iwall     = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
         if (Type_bc_elec[WallType[iwall]] == CONST_CHARGE ){

          charge_i = 0.0;
          for (idim=0; idim<Ndim; idim++)
             charge_i += Charge_w_sum_els[loc_inode][idim]*Area_surf_el[idim];
          resid[loc_i] -= 4.0*PI*charge_i/Temp_elec;

        }   /* check for charge b.c. on this wall */
     }      /* end of check for if this is a boundary node */
  return;
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
void load_nonlinear_transport_eqn(int i_box, int inode_box, int loc_i, int *ijk_box,
                        double *mat_row, double *resid, double *x,
                        int *bindx_tmp, int fill_flag, int iunk)
{

  int ilist=0, icomp, idim, ijk[3],inode;
  int iln, jln, elem, offset[3];
  int nodes_volm_el, nodes_surf_el, junk2[3];
  int off_ref[2][2] = { {0,1}, {-1,0}};
  int j_box_mu[8], j_box_rho[8], loc_j_mu[8], loc_j_rho[8], jnode_box,flag,count_flag;
  double wt=1.0, area_0=0.0, area_1=0.0, rho_0, rho_1, mu_0, mu_1, tmp;

   int dummy_int;

  /* pre calc basis function stuff for 3D */

  int igp;
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

   if (Ipot_ff_c == COULOMB) icomp = iunk - Ncomp - 1;
   else                      icomp = iunk - Ncomp;
   if (Matrix_fill_flag >= 3) icomp -= Nrho_bar;

   if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
   else if (Nlists_HW > 2)  ilist = icomp;

   node_to_ijk(node_box_to_node(inode_box),ijk);

   if (Zero_density_TF[inode_box][icomp]) {
       /* set mu to -VEXT_MAX in walls, and skip transport eqn fill */
       if (fill_flag != MSR_PREPROCESS){
          resid[loc_i] = x[loc_i] + VEXT_MAX;
          mat_row[loc_i] = 1.0;
       }
       else bindx_tmp[i_box]=TRUE;
   }
   /* check if you are in well-mixed bulk region */
   else if (ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
       if (fill_flag != MSR_PREPROCESS){
          resid[loc_i] = x[loc_i] - Betamu_LBB[icomp];
          mat_row[loc_i] = 1.0;
       }
       else bindx_tmp[i_box]=TRUE;
   }
   else if (Size_x[Grad_dim]-ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
       if (fill_flag != MSR_PREPROCESS){
          resid[loc_i] = x[loc_i] - Betamu_RTF[icomp];
          mat_row[loc_i] = 1.0;
       }
       else bindx_tmp[i_box]=TRUE;
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

             j_box_mu[jln]  = jnode_box*Nunk_per_node + iunk;   
             j_box_rho[jln] = jnode_box*Nunk_per_node + icomp;   
             if (fill_flag != MSR_PREPROCESS){
               loc_j_mu[jln]  = B2L_unknowns[j_box_mu[jln]];
               loc_j_rho[jln] = B2L_unknowns[j_box_rho[jln]];
             }
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
           if (fill_flag != MSR_PREPROCESS){

              if (Ndim==1) {   

                rho_0 = x[loc_j_rho[0]];
                rho_1 = x[loc_j_rho[1]];
                mu_0  = x[loc_j_mu[0]];
                mu_1  = x[loc_j_mu[1]];
                tmp = (2.0*rho_0*area_0 + rho_1*area_0 
                       + rho_0*area_1 + 2.0*rho_1*area_1);

                resid[loc_i] += wt*(mu_0 - mu_1)*tmp;

                mat_row[loc_j_rho[0]] += wt*(mu_0 - mu_1)*(2.0*area_0 + area_1);
                mat_row[loc_j_rho[1]] += wt*(mu_0 - mu_1)*(area_0 + 2.0*area_1);
                mat_row[loc_j_mu[0]]  += wt*tmp;
                mat_row[loc_j_mu[1]]  -= wt*tmp;
         
                /* add in a bulk convection term */
                if (iln == 0) tmp = (2.0*area_0 + area_1)/6.0;
                else          tmp = (area_0 + 2.0*area_1)/6.0;
                resid[loc_i] += Velocity * (rho_1 - rho_0) * tmp;
                mat_row[loc_j_rho[0]] -= Velocity * tmp;
                mat_row[loc_j_rho[1]] += Velocity * tmp;

              }
              else if (Ndim==3){     /* Ndim==3 case */

                /* loop over all quad points */

                for (igp=0; igp<8; igp++) { 
  
                  /* calculate quantities at quad point */
                  rho = 0.0; 
                  grad_mu[0] = grad_mu[1] = grad_mu[2] = 0.0;
                  for (jln=0; jln<8; jln++) {
/*if (Proc==4 && inode_box==L2B_node[26]) printf("jln=%d\n",jln);*/
dummy_int=jln;
                     rho += phi[jln][igp] * x[loc_j_rho[jln]];
                     grad_mu[0] += grad_phi[jln][igp][0] * x[loc_j_mu[jln]];
                     grad_mu[1] += grad_phi[jln][igp][1] * x[loc_j_mu[jln]];
                     grad_mu[2] += grad_phi[jln][igp][2] * x[loc_j_mu[jln]];
                  }
                  grad_mu_dot_grad_phi = grad_mu[0]*grad_phi[iln][igp][0]
                                       + grad_mu[1]*grad_phi[iln][igp][1]
                                       + grad_mu[2]*grad_phi[iln][igp][2];
                                      
                  resid[loc_i] += evol * rho * grad_mu_dot_grad_phi;
                  for (jln=0; jln<8; jln++) {
                    mat_row[loc_j_rho[jln]] += 
                                 evol * phi[jln][igp] * grad_mu_dot_grad_phi;
                    mat_row[loc_j_mu[jln]] += evol * rho * (
                                 grad_phi[jln][igp][0] * grad_phi[iln][igp][0]
                               + grad_phi[jln][igp][1] * grad_phi[iln][igp][1]
                               + grad_phi[jln][igp][2] * grad_phi[iln][igp][2]);
                  }
                }
              }
              else {
                printf("load_transport not written for 2D yet\n");
                exit(-1);
              }  
           }  
           else {
             for (jln=0; jln< Nnodes_per_el_V; jln++) {
               bindx_tmp[j_box_mu[jln]] =TRUE;
               bindx_tmp[j_box_rho[jln]]=TRUE;
             }
           }
           }
         }
       }
   }
}
/****************************************************************************/
void load_linear_transport_eqn(int i_box, int inode_box, int loc_i, int *ijk_box,
                        double *mat_row, double *resid, double *x,
                        int *bindx_tmp, int fill_flag, int iunk)
{

  int iwall, isten, idim,ijk[3],icomp;
  int iln, jln, elem, offset[3], loc_j,el_box;
  int nodes_volm_el, nodes_surf_el, junk2[3];
  /* static variables keep their value for every time the function is called*/
  static double *wt_lp_1el, *wt_s_1el;
  static int   **elem_permute, off_ref[2][2] = { {0,1}, {-1,0}};
  int j_box, jnode_box;

  /* First time through, load weights appropriate for this Ndim */
  
  if (MSR_PREPROCESS){
    set_fem_1el_weights(&wt_lp_1el, &wt_s_1el, &elem_permute);
  }

/* if we are at a boundary node and the boundary condition is constant
 * potential, then solve the equation psi = constant ..... in all other
 * cases loop over the elements near this node and fill the laplace and
 * source terms.  Note that the surface charge boundary condition is
 * handled in the main fill routine. */

   /* iunk is nodal unknown number, icomp is component number */

   if (Ipot_ff_c == COULOMB) icomp = iunk - Ncomp - 1;
   else                      icomp = iunk - Ncomp;
   if (Matrix_fill_flag >= 3) icomp -= Nrho_bar;

   node_to_ijk(node_box_to_node(inode_box),ijk);

   iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box];
   if (iwall != -1) {
       if (fill_flag != MSR_PREPROCESS){
          resid[loc_i] = x[loc_i] + 25.0;
          mat_row[loc_i] = 1.0;
       }
       else bindx_tmp[i_box]=TRUE;
   }
   /*else if (ijk[Grad_dim] == 0) {*/
    if (ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
       if (fill_flag != MSR_PREPROCESS){
          resid[loc_i] = x[loc_i] - Betamu_LBB[icomp];
          mat_row[loc_i] = 1.0;
       }
       else bindx_tmp[i_box]=TRUE;
   }
   /*else if (ijk[Grad_dim] == Nodes_x[Grad_dim]-1) {*/
   else if (Size_x[Grad_dim]-ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu+0.0000001) {
       if (fill_flag != MSR_PREPROCESS){
          resid[loc_i] = x[loc_i] - Betamu_RTF[icomp];
          mat_row[loc_i] = 1.0;
       }
       else bindx_tmp[i_box]=TRUE;
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
             j_box = jnode_box*Nunk_per_node + iunk;

             /*
              * add in Laplace term
              */
             if (jnode_box >= 0){  /* new flag for boundaries */
                if (fill_flag != MSR_PREPROCESS){
                   loc_j = B2L_unknowns[j_box];
                   resid[loc_i] += wt_lp_1el[isten]*x[loc_j];
                   mat_row[loc_j] += wt_lp_1el[isten];
                }
                else bindx_tmp[j_box]=TRUE;
             }


           }     /* end of loop over all local nodes in this fluid element */
         }       /* end of test to be sure element is in fluid */
       }         /* end of possible local node positions */
   }

}
/****************************************************************************/
/*void load_stoichiometric_constraint(int i_box, int inode_box, 
                        int loc_i, int *ijk_box, double *mat_row, 
                        double *resid, double *x, int *bindx_tmp, 
                                         int fill_flag, int iunk)
{

    

}*/
/****************************************************************************/
