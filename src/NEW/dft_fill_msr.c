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
 *  FILE: dft_fill_msr.c
 *
 *  This file contains the routines for putting the rows into
 *  the MSR format and determining nonzeros during preprocessing
 */

#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
void put_row_in_msr(int i_box, int loc_i,int *bindx, 
                    int *bindx_tmp, int **bindx_2d,
                    double *val, double *mat_row, 
                    int fill_flag,int *first_time)
{
int i,j,k, nonzeros_in_row, *b2dloci;
/* following needed for ideal gas */
int uniq, new_nonzero;

 if (fill_flag == MSR_PREPROCESS) {

    /* Make last row dense for parallel */
    /* TEMPORARY FIX FOR B2L_unknowns not containing correct unknowns ....
       will probably remove when the preprocessing is clean !! */

     if (*first_time && Num_Proc>1)  {
       for (j=0; j<Nunknowns_box; j++) bindx_tmp[j] = TRUE;
       *first_time = FALSE;
     }
    
     /* Count the number of nonzero off-diagonals */

     bindx_tmp[i_box] = FALSE;  /* don't count diagonal */
     nonzeros_in_row = 0;

     for (j=0; j<Nunknowns_box; j++) 
       if (bindx_tmp[j]) {
         bindx_tmp[nonzeros_in_row] = j;
         if (j>nonzeros_in_row) bindx_tmp[j] = FALSE;
         nonzeros_in_row++;
       }

     /* Allocate this row of nonzeros */
     bindx_2d[loc_i] = (int *) array_alloc(1, nonzeros_in_row, sizeof(int));

     /* Set up nonzero column indicies in ascending order */
     /* bindx_tmp values are box unknowns */

     b2dloci = bindx_2d[loc_i];

     if (!Non_unique_G2B[3]) {
       for (j=0; j<nonzeros_in_row; j++) {
         b2dloci[j] = B2G_unk[bindx_tmp[j]];
         bindx_tmp[j] = FALSE;
       }
     }

     /* if periodic BC cause nonunique G2B map -- use a search */

     else {
       new_nonzero = 0;
       for (j=0; j<nonzeros_in_row; j++) {
         i = B2G_unk[bindx_tmp[j]];
         bindx_tmp[j] = FALSE;

         /* if i is new, put it in b2dloci */

         uniq = TRUE;
         for (k=0; k<new_nonzero; k++) {
           if (b2dloci[k] == i) {uniq=FALSE; continue;}
         }
         if (uniq) {
           b2dloci[new_nonzero++] = i;
         }
       }
       nonzeros_in_row = new_nonzero;
     }

     /* Set nonzero counter for row and running total */

     bindx_2d[Aztec.N_update][loc_i] = nonzeros_in_row;
     Aztec.nonzeros += nonzeros_in_row;

  }
  else {
     val[loc_i] = mat_row[loc_i];
     mat_row[loc_i] = 0.0;
     for (j=bindx[loc_i]; j<bindx[loc_i+1]; j++) {
       val[j] = mat_row[bindx[j]];
       mat_row[bindx[j]] = 0.0;
     }
  }
  return;
}
/****************************************************************************/
void put_zero_in_msr(int loc_i,int **bindx_2d)
{
int nonzeros_in_row;

  /* Count the number of nonzero off-diagonals */
  nonzeros_in_row = 0;

  /* Allocate this row of nonzeros */
  bindx_2d[loc_i] = (int *) array_alloc(1, nonzeros_in_row, sizeof(int));

  /* Set nonzero counter for row and running total */

  bindx_2d[Aztec.N_update][loc_i] = nonzeros_in_row;
  Aztec.nonzeros += nonzeros_in_row;

  return;
}
/****************************************************************************/
/*put_euler_lag_in_msr: put row in msr for euler-lagrange term for which we always 
                    know where the nonzeros are found !!*/
void put_euler_lag_in_msr(int i_box, int loc_i, int **bindx_2d)
{
 int nonzeros_in_row;
 /* following needed for ideal gas */
 int ijk_box[3], icomp;
 int inode_box;
 int imu;

     
   inode_box = i_box/Nunk_per_node;
   icomp = i_box - inode_box*Nunk_per_node;
   node_box_to_ijk_box(inode_box,ijk_box);

                                               nonzeros_in_row = 0;
   if (Lsteady_state)                          nonzeros_in_row++;
   if (Ipot_ff_c == COULOMB)                   nonzeros_in_row++;

   bindx_2d[loc_i] = (int *) array_alloc(1, nonzeros_in_row, sizeof(int));

   /* put in electric potential entry */
   if (Ipot_ff_c == COULOMB) 
       bindx_2d[loc_i][0] = B2G_unk[i_box - icomp + Ncomp];

   /* put in chemical potential entry */
   if (Lsteady_state){
      if (Ipot_ff_c == COULOMB) imu = Ncomp +1;
      else                      imu = Ncomp;
      if (Ipot_ff_c == COULOMB) bindx_2d[loc_i][1] = B2G_unk[i_box + imu];
      else                      bindx_2d[loc_i][0] = B2G_unk[i_box + imu];
   }

   bindx_2d[Aztec.N_update][loc_i] = nonzeros_in_row;
   Aztec.nonzeros += nonzeros_in_row;

 return;

}
/****************************************************************************/
/*put_coarse_in_msr: put coarsened rows in msr ... we always 
                    know where the nonzeros are found !!*/
void put_coarse_in_msr(int mesh_coasen_flag, int i_box, 
                                    int loc_i, int **bindx_2d)
{

  int k, nonzeros_in_row;
  int idim, ijk_box[3], reflect_flag[3]={FALSE,FALSE,FALSE}, 
      icomp, offset[3]={0,0,0}, isten;
  int inode_box, jnode_box;

 inode_box = i_box/Nunk_per_node;
 icomp = i_box - inode_box*Nunk_per_node;
 node_box_to_ijk_box(inode_box,ijk_box);

 idim = -mesh_coasen_flag - 1;

 nonzeros_in_row = 2;
 bindx_2d[loc_i] = (int *) array_alloc(1, nonzeros_in_row, sizeof(int));

 k = 0;

 for (isten = -1; isten <= 1; isten += 2){
      offset[idim] = isten;
      jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
      if (jnode_box >=0 && 
          reflect_flag[0]+reflect_flag[1]+reflect_flag[2]==FALSE)
          bindx_2d[loc_i][k++] = B2G_unk[jnode_box*Nunk_per_node+icomp];
 }

 nonzeros_in_row = k;
 bindx_2d[Aztec.N_update][loc_i] = nonzeros_in_row;
 Aztec.nonzeros += nonzeros_in_row;

 return;
}
/****************************************************************************/
/*put_1Dsolution_in_msr: put 1D solution rows in msr ... we always 
                    know where the nonzeros are found !!*/
void put_1Dsolution_in_msr(int mesh_coasen_flag, int i_box, 
                         int loc_i, int jnode_box,int **bindx_2d)
{

  int nonzeros_in_row;
  int ijk_box[3], icomp;
  int inode_box;

 inode_box = i_box/Nunk_per_node;
 icomp = i_box - inode_box*Nunk_per_node;
 node_box_to_ijk_box(inode_box,ijk_box);


 nonzeros_in_row = 1;
 bindx_2d[loc_i] = (int *) array_alloc(1, nonzeros_in_row, sizeof(int));

 bindx_2d[loc_i][0] = B2G_unk[jnode_box*Nunk_per_node+icomp];
 bindx_2d[Aztec.N_update][loc_i] = nonzeros_in_row;
 Aztec.nonzeros += nonzeros_in_row;

 return;
}
/****************************************************************************/
/*put_transport_in_msr: put row in msr for poisson term for which we always 
                    know where the nonzeros are found !!*/
void put_transport_in_msr(int i_box, int loc_i, int **bindx_2d)
{

 int k, nonzeros_in_row;
 /* following needed for ideal gas */
 int ijk_box[3], reflect_flag[3]={FALSE,FALSE,FALSE}, 
     icomp, offset[3], isten, jsten, ksten;
 int inode_box, jnode_box;
 int minimum[3]={0,0,0}, maximum[3]={0,0,0};
 int irho,ijk[3];

   inode_box = i_box/Nunk_per_node;
   icomp = i_box - inode_box*Nunk_per_node;
   node_box_to_ijk_box(inode_box,ijk_box);

                  minimum[0] = -1; maximum[0] = 1;
   if (Ndim >1)  {minimum[1] = -1; maximum[1] = 1;}
   if (Ndim == 3){minimum[2] = -1; maximum[2] = 1;}

   k = 0;

   if (Linear_transport){
     nonzeros_in_row = POW_INT(3, Ndim) - 1;
     bindx_2d[loc_i] = (int *) array_alloc(1, nonzeros_in_row, sizeof(int));

     for (isten = minimum[0]; isten <= maximum[0]; isten++)
     for (jsten = minimum[1]; jsten <= maximum[1]; jsten++)
     for (ksten = minimum[2]; ksten <= maximum[2]; ksten++) {
          offset[0] = isten; offset[1] = jsten; offset[2] = ksten;
          jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
          if (jnode_box >=0 && 
                 reflect_flag[0]+reflect_flag[1]+reflect_flag[2]==FALSE){
          node_to_ijk(node_box_to_node(inode_box),ijk);
          if (!(ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu) &&
              !(Size_x[Grad_dim]-ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu)) {
              if (!(isten == 0 && jsten==0 && ksten==0))
                 bindx_2d[loc_i][k++] = B2G_unk[jnode_box*Nunk_per_node+icomp];
          }
          }
     }
   }
   else{
    
     nonzeros_in_row = 2*POW_INT(3, Ndim);
     bindx_2d[loc_i] = (int *) array_alloc(1, nonzeros_in_row, sizeof(int));

     for (isten = minimum[0]; isten <= maximum[0]; isten++)
     for (jsten = minimum[1]; jsten <= maximum[1]; jsten++)
     for (ksten = minimum[2]; ksten <= maximum[2]; ksten++) {
          offset[0] = isten; offset[1] = jsten; offset[2] = ksten;
          jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
          if (jnode_box >=0 && 
                 reflect_flag[0]+reflect_flag[1]+reflect_flag[2]==FALSE){
          node_to_ijk(node_box_to_node(inode_box),ijk);
          if (!(ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu) && 
              !(Size_x[Grad_dim]-ijk[Grad_dim]*Esize_x[Grad_dim] <= X_const_mu)) {
              if (Ipot_ff_c == COULOMB) irho = icomp-Ncomp-1;
              else                      irho = icomp-Ncomp;
              bindx_2d[loc_i][k++] = B2G_unk[jnode_box*Nunk_per_node+irho];
              if (!(isten == 0 && jsten==0 && ksten==0))
                 bindx_2d[loc_i][k++] = B2G_unk[jnode_box*Nunk_per_node+icomp];
          }
        }
     }
   }

   nonzeros_in_row = k;
   bindx_2d[Aztec.N_update][loc_i] = nonzeros_in_row;
   Aztec.nonzeros += nonzeros_in_row;

   return;
}
/****************************************************************************/
/*put_poisson_in_msr: put row in msr for poisson term for which we always 
                    know where the nonzeros are found !!*/
void put_poisson_in_msr(int i_box, int loc_i, int **bindx_2d)
{

 int j,k, nonzeros_in_row;
 /* following needed for ideal gas */
 int idim, jdim, ijk_box[3], reflect_flag[3]={FALSE,FALSE,FALSE}, 
     icomp, offset[3], isten, jsten, ksten;
 int inode_box, jnode_box,jtmp;
 int minimum[3]={0,0,0}, maximum[3]={0,0,0};

   inode_box = i_box/Nunk_per_node;
   icomp = i_box - inode_box*Nunk_per_node;
   node_box_to_ijk_box(inode_box,ijk_box);

                  minimum[0] = -1; maximum[0] = 1;
   if (Ndim >1)  {minimum[1] = -1; maximum[1] = 1;}
   if (Ndim == 3){minimum[2] = -1; maximum[2] = 1;}

  nonzeros_in_row = 2*Ndim + (Ncomp+1)*POW_INT(3, Ndim) - 1;
  k = 0;
  bindx_2d[loc_i] = (int *) array_alloc(1, nonzeros_in_row, sizeof(int));

  for (isten = minimum[0]; isten <= maximum[0]; isten++)
  for (jsten = minimum[1]; jsten <= maximum[1]; jsten++)
  for (ksten = minimum[2]; ksten <= maximum[2]; ksten++) {
      offset[0] = isten; offset[1] = jsten; offset[2] = ksten;
      jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
      if (jnode_box >=0 &&
          reflect_flag[0]+reflect_flag[1]+reflect_flag[2]==FALSE)
            for (j=0; j <= Ncomp; j++){
              if (!(isten==0 && jsten==0 && ksten==0 && j==Ncomp )){
                if (Type_poly==-1){
                   if (j==Ncomp) jtmp=Ncomp+Nrho_bar;
                   else jtmp=j;
                }
                else{
                   if (j==Ncomp) jtmp=2*Ncomp+Ngeqn_tot;
                   else jtmp=Ncomp+j;
                }
                bindx_2d[loc_i][k++] = B2G_unk[jnode_box*Nunk_per_node+jtmp];
              }
            }
  }


  /* Jacobian entries for post-processing of derivatives ... these
    are zeros that will be carried along in the solve !! */
  for (idim=0; idim<Ndim; idim++){
     for (jdim=0; jdim<Ndim; jdim++) offset[jdim] = 0;
  
     for (isten=-2; isten<3; isten+=4){
         if (isten != 0) {
            offset[idim] = isten;
            jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
            if (jnode_box >=0 && 
                reflect_flag[0]+reflect_flag[1]+reflect_flag[2]==FALSE)
               if (Type_poly==-1){
                   bindx_2d[loc_i][k++] = B2G_unk[jnode_box*Nunk_per_node+Ncomp+Nrho_bar];
               }
               else{ bindx_2d[loc_i][k++] = B2G_unk[jnode_box*Nunk_per_node+2*Ncomp+Ngeqn_tot];
               }
         }
    }
  }

  nonzeros_in_row = k;
  bindx_2d[Aztec.N_update][loc_i] = nonzeros_in_row;
  Aztec.nonzeros += nonzeros_in_row;
  
}
/****************************************************************************/
