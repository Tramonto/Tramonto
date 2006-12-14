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

/*#include <stdio.h>
#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"*/

#include "dft_uww.h"

#define FLAG_LJ 0
#define FLAG_COULOMB 1


/* Prototypes for functions found in this file */

/*void setup_lj_atomic(int,int);
void setup_coulomb_atomic(int,int);
void setup_ww_integrated(int,int *,int **,int,int);*/


/******************************************************************************/

void setup_wall_wall_potentials( int **nelems_w_per_w,
                            int ***elems_w_per_w) 

/* Set up any wall-wall potentials we may need for PMF calculations.
 *
 *    Authors:  Laura Frink
 */

{
   /* Local variable declarations */
   
   char *yo = "wall_wall_potentials";
   int ilist,iwall,jwall,itype_w,jtype_w;
   int **elems_w_per_w_global, *nelems_w_per_w_global;
   int L_LJ,L_COULOMB;
   FILE *fp10, *fp11;
  
  /********************** BEGIN EXECUTION ************************************/

  if (Proc==0) printf("\n %s: Setting up Wall-Wall Potentials ... \n",yo);
  /*
   * Allocate and zero the arrays we will calculate here
   */

  if( (fp10  = fopen("dft_uww.dat","w")) == NULL) {
    printf("Can't open file dft_uww.dat\n");
    exit(1);
  }
  if( (fp11  = fopen("dft_uww_link.dat","w")) == NULL) {
    printf("Can't open file dft_uww.dat\n");
    exit(1);
  }

   Uww = (double **) array_alloc (2,Nwall,Nwall, sizeof(double));
   Uww_link = (double **) array_alloc (2,Nlink,Nlink, sizeof(double));

   for (iwall=0; iwall<Nwall; iwall++){
      for (jwall=0; jwall<Nwall; jwall++) {
          Uww[iwall][jwall]=0.0;
          Uww_link[Link[iwall]][Link[jwall]]=0.0;
       }
   }

  /*
   * For each desired case, write a new subroutine,
   * add the case as a choice to the input file,
   * and add a case statement here.    Note that computing wall-wall potentials
   * is only done for the cases where the walls are represented by atomic potentials.
   * So, for atomic 12-6 LJ walls or smeared colloidal LJ walls, or Coulomb walls
   * we can compute these direct interactions.  Note that if a charge per unit
   * area will be applied to a surface, the Uww arrays will be modified after the
   * boundary conditions are set up (see boundary_setup routines).  We can compute
   * volume charge contributions here because we have access to the elems_w_per_w
   * and nelems_per_w arrays.  Some of geometry code is repeated  when the charged
   * boundaries are set up in the case of volume charges.  Clearly both sections
   * must find the same wall elements for each wall or else there will be a
   * discrepancy in Uww computed here and the Uww that would be based on the
   * surface charge distribution computed later.
   */

  for (iwall=0; iwall<Nwall-1; iwall++){
     itype_w=WallType[iwall];
     for (jwall=iwall; jwall<Nwall; jwall++){
        jtype_w=WallType[jwall];
        if (Ipot_ww_n[itype_w][jtype_w]!=NO_WW){

  /* if necessary define and communicate elems_w_per_w_global */
/*  if (Ipot_ww_n[itype_w][jtype_w]==SMEAR_POINT_WW){

       if (Ipot_wf_n[itype_w]==VEXT_3D_INTEGRATED) nwall=iwall; 
       else if (Ipot_wf_n[jtype_w]==VEXT_3D_INTEGRATED) nwall=jwall; 
       else {
          printf("Error setting up wall-wall potentials: can't continue ...\n");
          printf("need to know which surface has smeared out LJ12-6 potentials\n");
          printf("Ipot_ww_n[%d][%d]=%d:  Ipot_wf_n[%d]=%d  and  Ipot_wf_n[%d]=%d\n",
                   itype_w,jtype_w,Ipot_ww_n[itype_w][jtype_w],itype_w,Ipot_wf_n[itype_w],
                   jtype_w,Ipot_wf_n[jtype_w]);
          exit(-1);
       }

       elems_w_per_w_global = (int **) array_alloc(1,Nlists_HW, sizeof(int *));
       nelems_w_per_w_global = (int *) array_alloc(1,Nlists_HW, sizeof(int));

       comm_wall_els(nwall,nelems_w_per_w,elems_w_per_w, nelems_w_per_w_global, elems_w_per_w_global);
  }
*/

  /* compute Lennard-Jones part of wall-wall interaction potentials */
  if (Ipot_ww_n[itype_w][jtype_w]==ATOM_CENTERS_WW){
     setup_lj_atomic(iwall,jwall);
     if (Type_coul >= 0){ 
          setup_coulomb_atomic(iwall,jwall);
     }
  }
/*  else if (Ipot_ww_n[itype_w][jtype_w]==SMEAR_POINT_WW){
       L_LJ=TRUE; L_COULOMB=FALSE;
       setup_ww_integrated(nwall,nelems_w_per_w_global,elems_w_per_w_global,L_LJ,L_COULOMB);
  }
*/

  /* compute Coulomb part of wall-wall interaction potentials. */
/*  if (Type_coul >= 0){ 
     if (Charge_type_atoms == SMEAR_CHARGE){
            L_LJ = FALSE;  L_COULOMB = TRUE;
            setup_ww_integrated(nwall,nelems_w_per_w_global,
                                elems_w_per_w_global,L_LJ,L_COULOMB);
     }
     if (Nlocal_charge >0){
        printf ("need to do some work to compute the energies of local charges\n");
        printf ("that are not associated with any particular surface !!\n");
        exit(-1);
     }
  }
*/

/*  if (Ipot_ww_n[itype_w][jtype_w]==SMEAR_POINT_WW){
       for (ilist=0; ilist<Nlists_HW; ilist++){
             safe_free((void *) &elems_w_per_w_global[ilist]);
       }
       safe_free((void *) &elems_w_per_w_global);
       safe_free((void *) &nelems_w_per_w_global);
  }
*/
  } /* end if not hard or other noninteracting walls */
  } /* end loop over jwall */
  } /* end loop over iwall */

  for (iwall=0; iwall<Nwall-1; iwall++){
     for (jwall=iwall; jwall<Nwall; jwall++){
         if (Link[iwall] != Link[jwall])
         Uww_link[Link[iwall]][Link[jwall]] += Uww[iwall][jwall];
         Uww_link[Link[jwall]][Link[iwall]] += Uww[iwall][jwall];
         fprintf(fp10," %d  %d  %9.6f\n",iwall,jwall,Uww[iwall][jwall]);
     }
  }

  if (Nwall != Nlink){
     for (iwall=0; iwall<Nlink-1; iwall){
         for (jwall=iwall; jwall<Nlink; jwall){
            fprintf(fp11," %d  %d  %9.6f\n",iwall,jwall,Uww_link[iwall][jwall]);
         }
      }
  }

  fclose (fp10);
  fclose (fp11);
  return;
}
/******************************************************************************/
/* setup_lj_atomic: In this routine we base wall-wall potentials on the
                            centers of the atomic particles that make up the
                            surfaces.  */
void setup_lj_atomic(int iwall, int jwall){


  int idim;
  double xi,xj,rsq,r;

  rsq = 0.0;
  for (idim=0; idim<Ndim; idim++){
     xi = WallPos[idim][iwall];
     xj = WallPos[idim][jwall];
     rsq += (xi-xj)*(xi-xj);
  }
  r = sqrt (rsq);

  /* add in Lennard-Jones Contributions */
  Uww[iwall][jwall] += uLJ12_6_cut(r,Sigma_ww[WallType[iwall]][WallType[jwall]],
                                   Eps_ww[WallType[iwall]][WallType[jwall]],
                                   Cut_ww[WallType[iwall]][WallType[jwall]]);
  return;
}
/******************************************************************************/
/* setup_coulomb_atomic: In this routine we base wall-wall potentials on the
                            centers of the atomic particles that make up the
                            surfaces.  */
void setup_coulomb_atomic(int iwall,int jwall){


  int idim,index1,index2;
  double xi,xj,rsq,r,faci,facj;

  faci = 1.0; facj=1.0;
  for (idim = 0; idim<Ndim; idim++) {
     if ((WallPos[idim][iwall] == -0.5*Size_x[idim] && Type_bc[idim][0] == REFLECT) ||
         (WallPos[idim][iwall] == 0.5*Size_x[idim] && Type_bc[idim][1] == REFLECT)) faci *= 2.0;
     if ((WallPos[idim][jwall] == -0.5*Size_x[idim] && Type_bc[idim][0] == REFLECT) ||
         (WallPos[idim][jwall] == 0.5*Size_x[idim] && Type_bc[idim][1] == REFLECT)) facj *= 2.0;
  }

   rsq = 0.0;
   for (idim=0; idim<Ndim; idim++){
       xi = WallPos[idim][iwall];
       xj = WallPos[idim][jwall];
       rsq += (xi-xj)*(xi-xj);
   }
   r = sqrt (rsq);

         /* add in Coulomb Contributions */

   Uww[iwall][jwall] += uCOULOMB(r,faci*Elec_param_w[iwall],facj*Elec_param_w[jwall]);

  return;
}
/******************************************************************************/
/* setup_ww_integrated: In this routine we assume that materials properties
                            are constant in the surfaces of interest, and that
                            surface-surface interactions are described by either 12-6
                            Lennard-Jones or Coulomb potentials.  The total wall-wall potentials
                            are calculated by integrating on an element by element
                            basis through pairs of walls */
void setup_ww_integrated(int iwall,int *nelems_w_per_w,int **elems_w_per_w,int L_LJ, int L_COULOMB)
{
  int iwall_type,ilist,idim, i;
   int jwall,jwall_type=0,jnode_box,jnode_llb,ig,jg,kg,jwall_el;
   int ngpi=0,ngpj=0,ngpk=0;
   int inode_llb,iel,wall_el;
   int max_els,image,image_x,image_xy,iel_y,iel_z;
   int ngp,ngpu,ngp1,ngpu1,ngp2,ngpu2,ngp3,ngpu3;
   int imj,jimage,iside,ijk_j[3];
   double pos_jw[3],**jimage_pos;
   double factor;
   double *gp,*gw,*gpu,*gwu;
   double gp1[12],gw1[12],gpu1[40],gwu1[40];
   double gp2[12],gw2[12],gpu2[12],gwu2[12];
   double gp3[12],gw3[12],gpu3[12],gwu3[12];
   double max_cut,**image_pos,node_pos[3],node_pos_w[3],
          point[3],node_pos_w2[3],u_ww,r_center_sq,
          node_pos_jw[3];


   max_cut = 0.0;
   for (iwall_type=0; iwall_type<Nwall_type; iwall_type++)
     for (jwall_type=0; jwall_type<Nwall_type; jwall_type++)
         if (max_cut < Cut_ww[iwall_type][jwall_type]) 
             max_cut = Cut_ww[iwall_type][jwall_type];
   max_els = POW_DOUBLE_INT(max_cut,Ndim)/Vol_el + 1;


   image_pos = (double **) array_alloc (2, max_els, Ndim, sizeof(double));
   jimage_pos = (double **) array_alloc (2, max_els, Ndim, sizeof(double));

   ilist = Nlists_HW-1; /* only add the contributions of the solid*/

  /* put 12 gauss quadrature points in each element */
   gp = (double *) array_alloc (1, 12, sizeof(double));
   gw = (double *) array_alloc (1, 12, sizeof(double));
   gpu = (double *) array_alloc (1, 12, sizeof(double));
   gwu = (double *) array_alloc (1, 12, sizeof(double));

  ngp1  = 6;  ngp2 = 3; ngp3  = 3;
  ngpu1 = 20; ngpu2 = 12; ngpu3 = 6;

  set_gauss_quad(ngp1, gp1, gw1); 
  set_gauss_quad(ngp2, gp2, gw2); 
  set_gauss_quad(ngp3, gp3, gw3); 

  set_gauss_quad(ngpu1, gpu1, gwu1);
  set_gauss_quad(ngpu2, gpu2, gwu2);
  set_gauss_quad(ngpu3, gpu3, gwu3);

  /* consider all iwall elements and images */
  for (iwall=0; iwall<Nwall-1; iwall++){
    for (wall_el=0; wall_el< nelems_w_per_w[ilist]; wall_el++){
      if (Proc==0 && Iwrite==VERBOSE) 
         printf("iwall: %d  of Nwall: %d wall_el: %d  of nelems_w_per_w: %d\n",
                                               iwall,Nwall,wall_el,nelems_w_per_w[0]);
      iel = elems_w_per_w[ilist][wall_el];
      inode_llb = element_to_node(iel);

      iwall_type = WallType[iwall];

      node_to_position (inode_llb,node_pos);
      for (idim=0; idim<Ndim; idim++) {
        node_pos_w[idim] = node_pos[idim] + 0.5*Esize_x[idim];
        image_pos[0][idim] = node_pos_w[idim];
      }

      /* consider only jwall elements in the solution volumes */
      for (jwall=iwall+1; jwall<Nwall; jwall++){
        if (Link[iwall] != Link[jwall]){
        for (jwall_el=0; jwall_el< nelems_w_per_w[ilist]; jwall_el++){
        jnode_llb = element_to_node(elems_w_per_w[ilist][jwall_el]);
        jnode_box = node_to_node_box(jnode_llb);
        node_to_ijk(jnode_llb,ijk_j);

        /* find images of the jnode if in wall boundaries apply */
        node_to_position(jnode_llb,node_pos_jw);
        for (idim=0; idim<Ndim; idim++) {
          pos_jw[idim] = node_pos_jw[idim] + 0.5*Esize_x[idim];
          jimage_pos[0][idim] = pos_jw[idim];
        }
        jimage = 1;

        iside=-1;
        if (ijk_j[0] == 0 && Type_bc[0][0]==IN_WALL) iside=0;
        else if (ijk_j[0] == Nodes_x[0]-1 && Type_bc[0][1]==IN_WALL) iside=1;
        
        if (iside !=-1) find_images2(0,max_cut,&jimage,jimage_pos,pos_jw,iwall,iside);

     for (imj=0; imj<jimage; imj++){
       
        r_center_sq = 0.0;
        for (idim=0; idim<Ndim; idim++) {
            r_center_sq += (node_pos_w[idim]-jimage_pos[imj][idim])*
                           (node_pos_w[idim]-jimage_pos[imj][idim]);
        }
        if (r_center_sq < 4.0){
            ngp = ngp1; ngpu=ngpu1; 
            gp = &gp1[0]; gpu = &gpu1[0]; 
            gw = &gw1[0]; gwu = &gwu1[0];
         }
         else if (r_center_sq < 16.0){
            ngp = ngp2; ngpu=ngpu2; 
            gp = &gp2[0]; gpu = &gpu2[0]; 
            gw = &gw2[0]; gwu = &gwu2[0];
         }
         else{
            ngp = ngp3; ngpu=ngpu3; 
            gp = &gp3[0]; gpu = &gpu3[0]; 
            gw = &gw3[0]; gwu = &gwu3[0];
         }

        switch(Ndim){
           case 3:  ngpi=ngp; ngpj = ngp; ngpk = ngp; break;
           case 2:  ngpi=ngp; ngpj = ngp; ngpk = 1; break;
           case 1:  ngpi=ngp; ngpj = 1; ngpk = 1; break;
        }

        for (kg=0; kg < ngpk; kg++) {
        for (jg=0; jg < ngpj; jg++) {
        for (ig=0; ig < ngpi; ig++) {
         
           factor = Vol_el; 
           switch(Ndim){
              case 3: point[2] = jimage_pos[imj][2] + gp[kg] * Esize_x[2];
                      factor *= gw[kg];
              case 2: point[1] = jimage_pos[imj][1] + gp[jg] * Esize_x[1];
                      factor *= gw[jg];
              case 1: point[0] = jimage_pos[imj][0] + gp[ig] * Esize_x[0]; 
                      factor *= gw[ig];
           }

           image = 1;
           if (Uww[iwall][jwall] < VEXT_MAX) {

           find_images(0,max_cut,&image,image_pos,node_pos_w,point);

           if (Ndim > 1) {
             image_x = image;
             for (iel_y=0; iel_y<image_x; iel_y++){
                for (idim=0; idim<Ndim; idim++) 
                      node_pos_w2[idim] = image_pos[iel_y][idim];
                find_images(1,max_cut,&image,image_pos,node_pos_w2,point);
             }
          }

          if (Ndim == 3) {
             image_xy = image;
             for (iel_z=0; iel_z<image_xy; iel_z++){
                for (idim=0; idim<Ndim; idim++) 
                    node_pos_w2[idim] = image_pos[iel_z][idim];
                find_images(2,max_cut,&image,image_pos,node_pos_w2,point);
             }
          }

          /* 
           *  now that we have all the images associated with
           *  this particular surface element for icomp...calculate
           *     the wall-wall interactions !!
           */
          
          for (i=0; i<image; i++){
             for (idim=0; idim<Ndim; idim++)
                   node_pos_w2[idim] = image_pos[i][idim];

             if (L_LJ){
                u_ww = integrate_potential(FLAG_LJ,Sigma_ww[iwall_type][jwall_type],
                       Eps_ww[iwall_type][jwall_type], Cut_ww[iwall_type][jwall_type],
                       ngp, ngpu, gp, gpu, gw, gwu, node_pos_w2, point);

                Uww[iwall][jwall] += u_ww*factor;
             }

             if (L_COULOMB && Type_bc_elec[WallType[iwall]]==ATOMIC_CHARGE &&
                              Type_bc_elec[WallType[jwall]]==ATOMIC_CHARGE){
                u_ww = integrate_potential(FLAG_COULOMB, Elec_param_w[iwall],
                                          Elec_param_w[jwall], Cut_ww[iwall_type][jwall_type],
                                          ngp, ngpu, gp, gpu, gw, gwu, node_pos_w2, point);

                Uww[iwall][jwall] += u_ww*factor;
             }

          }  /* end of loop over the surface elements and its images */
          if (Uww[iwall][jwall] >= VEXT_MAX) Uww[iwall][jwall] = VEXT_MAX;

         }   /* if test on Uww*/

       }   /* i - gauss points loop */
       }   /* j - gauss points loop */
       }   /* k - gauss points loop */

    }  /* jwall_image loop */
      
      }    /* jwall_el loop */
      Uww[jwall][iwall] = Uww[iwall][jwall];
      }   /* if iwall and jwall are not on same macrosurface */
      }    /* jwall loop */
    }      /* iwall_el loop */
  }        /* iwall loop */
  safe_free((void *) &image_pos);
  safe_free((void *) &jimage_pos);
  return;
}
/*****************************************************************************/

