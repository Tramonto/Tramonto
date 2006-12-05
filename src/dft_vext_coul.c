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
 *  FILE: dft_vext.c
 *
 *  This file contains routines to set up an external field for a 
 *  variety of cases.
 */

#include <stdio.h>
#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"



/* Prototypes for functions found in this file */
void comm_loc_to_glob_vec(int *,int *,int *);

#define FLAG_COULOMB 1

void setup_vext_coulomb_vol()
{
   double *charge_el_vol_global, *charge_el_vol_loc, *unk_global,cut_wf_max;
  int iel_box,inode_box,n_loc,reflect_flag[3], *index, inode,i,iel;
  int iwall_type,icomp, idim, loc_inode;
  double elem_pos[3],elem_pos2[3],x,r,rsq;
  double node_pos[3],node_pos_f[3],pos[3];

  double **image_pos,r_center_sq;
  int max_images; 

  int image,image_x,image_xy,iel_y,iel_z;
  int ngp,ngp1,ngp2;
  double *gp,*gw;
  double gp1[12],gw1[12];
  double gp2[12],gw2[12];

   reflect_flag[0]=FALSE; reflect_flag[1]=FALSE; reflect_flag[2]=FALSE;

   /* allocate coulomb external field and external field derivative arrays. */
   Vext_coul = (double **) array_alloc (2,Nnodes_per_proc,Ncomp, sizeof(double));
   /*Vext_dash =  (double ***) array_alloc (3, Nnodes_per_proc,Ncomp*Nlink, Ndim, sizeof(double));*/

   /* allocate temporary local and global vectors to store the charge_els array */
   charge_el_vol_global = (double *) array_alloc (1,Nelements, sizeof(double));
   charge_el_vol_loc = (double *) array_alloc (1,Nnodes_per_proc, sizeof(double));

   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
         charge_el_vol_loc[loc_inode]=0.0;
         for (icomp=0; icomp<Ncomp; icomp++) Vext_coul[loc_inode][icomp]=0.0;
   }
   for (iel=0; iel<Nelements; iel++){
         charge_el_vol_global[iel]=0.0;
   }


   if (Proc==0){
      index = (int *) array_alloc (1, Nnodes, sizeof(int));
      unk_global = (double *) array_alloc (1,Nnodes, sizeof(double));
      for (i=0; i<Nnodes; i++){
         index[i]=0; unk_global[i]=0.0;
      }
   }
   else{ 
     index=NULL;
     unk_global=NULL;
   }

   /* set up the local charge_el array*/

   for (iel_box=0; iel_box<Nelements_box; iel_box++){
       inode_box=element_box_to_node_box(iel_box);
       if (B2L_node[inode_box] >=0) {
            charge_el_vol_loc[B2L_node[inode_box]] = Charge_vol_els[iel_box];
       }
   }

   /* gather the L2G_node array and charge_el_vol_loc array on to processor 0 */

  /* collect the order of unknowns from all the processors */

   MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
              index,Comm_node_proc,Comm_offset_node,
              MPI_INT,0,MPI_COMM_WORLD);

  /* collect the charge_el_array from all the processors */
   MPI_Gatherv(charge_el_vol_loc,Nnodes_per_proc,MPI_DOUBLE,
              unk_global,Comm_node_proc,Comm_offset_node,
              MPI_DOUBLE,0,MPI_COMM_WORLD);

   safe_free((void *) &charge_el_vol_loc);

  /* now processor 0 must sort the charge_el_vol_loc input from other processors
     into a global array that is organized from 1 to Nelements */

  if (Proc == 0){
     for (i=0; i<Nnodes; i++){
         inode = index[i];
         iel = node_to_elem(inode,0,reflect_flag);
         if (iel >= 0){
           charge_el_vol_global[iel]=unk_global[i];
         }
     }
     safe_free((void *) &unk_global);
     safe_free((void *) &index);
  }

  /* now, broadcast the global charge els array so all the processors have the array */
    MPI_Bcast(charge_el_vol_global,Nelements,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* now that each processor has the entire fixed charge per element array, compute
     the external field on each node due to those fixed charges.  Take the center
     of the element as a basis for the 1/r calculation */

    /* need to find the wall images for reflective boundary conditions. 
       Note that while the same code used to find images generally is used
       here, this routine will not be called for periodic systems or multiply
       reflective systems ..... we are not set up for the infinite images
       needed to do periodic Coulomb systems properly */


   max_images=8;  /* element at corner of 3 reflective boundaries */
   image_pos = (double **) array_alloc (2, max_images, Ndim, sizeof(double));

   gp = (double *) array_alloc (1, 12, sizeof(double));
   gw = (double *) array_alloc (1, 12, sizeof(double));
   ngp1  = 6;  ngp2 = 3; 
   ngp1  = 1;  ngp2 = 1; 
   set_gauss_quad(ngp1, gp1, gw1);
   set_gauss_quad(ngp2, gp2, gw2);

   for (iel=0; iel<Nelements; iel++){
        if (fabs(charge_el_vol_global[iel]) > 1.e-8) {

            /* find center of charged element */
            inode = element_to_node(iel);
            node_to_position(inode,node_pos);
            for (idim=0; idim<Ndim; idim++){
                elem_pos[idim]=node_pos[idim]+0.5*Esize_x[idim];
                image_pos[0][idim] = elem_pos[idim];
            }
            image=1;
            /* find all possible images of this charged element */ 
            find_images_coulomb(0,&image,image_pos,elem_pos);
            image_x=image;
            for (iel_y=0; iel_y<image_x; iel_y++){
                for (idim=0; idim<Ndim; idim++)
                    elem_pos2[idim] = image_pos[iel_y][idim];
                find_images_coulomb(1,&image,image_pos,elem_pos2);
            }
            image_xy = image;
             for (iel_z=0; iel_z<image_xy; iel_z++){
                for (idim=0; idim<Ndim; idim++)
                    elem_pos2[idim] = image_pos[iel_z][idim];
                find_images_coulomb(2,&image,image_pos, elem_pos2);
            }

            /*loop over fluid nodes to set up external field for each node */
            for (loc_inode=0; loc_inode< Nnodes_per_proc; loc_inode++){
               inode_box = L2B_node[loc_inode];

               for (icomp=0; icomp<Ncomp; icomp++) {
               if (!Zero_density_TF[inode_box][icomp]) {
                  node_to_position(L2G_node[loc_inode],node_pos_f);

               for (i=0; i<image; i++){
                  r_center_sq = 0.0;
                  for (idim=0; idim<Ndim; idim++) {
                       pos[idim] = image_pos[i][idim]-0.5*Esize_x[idim];
                       r_center_sq += (pos[idim]-node_pos_f[idim])*
                                      (pos[idim]-node_pos_f[idim]);
                  }
                  if (r_center_sq < 4.0) { ngp = ngp1; gp = &gp1[0]; gw = &gw1[0]; }
                  else{                    ngp = ngp2; gp = &gp2[0]; gw = &gw2[0]; }

                  Vext_coul[loc_inode][icomp] += integrate_potential(FLAG_COULOMB,
                          Charge_f[icomp],charge_el_vol_global[iel],
                          1., ngp, 1, gp, NULL, gw, NULL, pos, node_pos_f);
               }
               }
                  if(Type_coul == LIKE_LJ) Vext[loc_inode][icomp] += Vext_coul[loc_inode][icomp];
            }

        }
    }
  }
}
