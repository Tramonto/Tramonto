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
 *  FILE: dft_mesh_surfaces.c
 *
 *  This file contains the routines used to set up the surface of interest.  All
 *  of the geometry dependent information is found here.  All other routines
 *  are oblivious to surface geometry.
 *  
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

/* Prototypes for functions found in this file */

void els_planar(int, int, int, int **, double *, int **, int ***,int *,int ***,double **);
void els_finite_planar(int,int, int,int **,double *, int **, int ***,int*,int ***,double **);
void els_spheres(int, int, int, int **, double *, int **, int ***, int *,int ***,double **);
void els_cyls_3D(int, int, int, int **, double *, int **, int ***, int *,int ***,double **);
void els_cyls_cos_3D(int, int, int, int **, double *, int **, int ***, int *,int ***,double **);
void els_atomic_centers(int, int, int, int **, double *, int **, int ***, int *,int ***,double **);
void els_cyl_pores(int, int, int, int **, double *, int **, int ***, int *,int ***,double **);
void els_slit_pore_2D(int, int, int, int **,double *,int **,int ***,int *,int ***,double **);
void els_cyl_pore_3D(int, int, int, int **,double *, int **, int ***, int *,int ***,double **);
void els_cone_pore_2D(int, int,int,  int **,double *,int **,int ***,int *,int ***,double **);
void els_cone_pore_3D(int, int,int,  int **,double *,int **,int ***,int *,int ***,double **);
void flag_wall_el(int,int,int,int, int **,int **,int ***,int ***);
void find_wall_images(int, int *, double **, double *);

void setup_surface (FILE *fp2, int *nelems_f,  
                    int **nelems_w_per_w, 
                    int **elems_f, 
                    int ***elems_w_per_w, int *elem_zones,
                    int *fast_fill_elem_TF,int ***el_type)
{
 /* Local variable declarations */
  int itype,iwall,ilist,iel_box,izone, **L_wall,idim,loc_inode,real_wall,*imagetype,flag;
  double *x_min;

  int image,image_x,image_xy,i,nwall_max,image_old,icount;
  double pos[3],**image_pos;
  
  /*  Calculate the number of wall and fluid elements and nodes, and 
      fill arrays that index the wall and fluid elements and nodes.  This
      is where the wall geometry comes in. */ 

  L_wall = (int **) array_alloc (2, Nlists_HW,Nelements_box, sizeof(int));
  x_min  = (double *) array_alloc (1, Nelements_box, sizeof(double));
  Touch_domain_boundary = (int ****) array_alloc (4,Nwall,Nlists_HW,Ndim,2,sizeof(int));
  flag=FALSE;
  for (i=0;i<Nwall_type;i++) if (Ipot_wf_n[i]==VEXT_1D_XMIN) flag=TRUE;
  if (flag){
         X_wall = (double **) array_alloc (2, Nnodes_per_proc, Nwall,sizeof(double));
         X_wall2 = (double **) array_alloc (2, Nnodes_per_proc, Nwall,sizeof(double));
     for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
      for (iwall=0; iwall<Nwall; iwall++)  X_wall[loc_inode][iwall]=-999.0;
      for (iwall=0; iwall<Nwall; iwall++)  X_wall2[loc_inode][iwall]=-999.0;
     }
  }

  for (iel_box=0; iel_box<Nelements_box; iel_box++){
     fast_fill_elem_TF[iel_box] = TRUE;
     if (Coarser_jac != 5) elem_zones[iel_box] = Nzone - 1;
     else                  elem_zones[iel_box] = Nzone - 2;
     x_min[iel_box] = 1000.;
     for (ilist=0; ilist<Nlists_HW; ilist++){
        L_wall[ilist][iel_box] = FALSE;
        Wall_elems[ilist][iel_box] = -1;
        for (iwall=0; iwall<Nwall; iwall++) el_type[iwall][ilist][iel_box] = FLUID_EL;
     }
  }

  for (iwall=0; iwall<Nwall; iwall++)
    for (ilist=0; ilist<Nlists_HW; ilist++){
       nelems_w_per_w[ilist][iwall] = 0;
       for (idim=0; idim<Ndim; idim++) {
         Touch_domain_boundary[iwall][ilist][idim][0]=FALSE;
         Touch_domain_boundary[iwall][ilist][idim][1]=FALSE;
    }
  }

 /* To treat wall-wall boundaries effectively and clean
    up post-processing, we need an array that is false
    if we have a surface that is infinite extent in a
    finite dimension ... e.g. cylinders in z .... */

  for (iwall=0; iwall<Nwall; iwall++){
     itype = WallType[iwall];
     if(Ndim > 1){
        if(Surface_type[itype] == smooth_planar_wall) {
            for (idim=0; idim<Ndim; idim++) {
               if (idim != Orientation[itype]) Xtest_reflect_TF[Link[iwall]][idim] = FALSE;
            } 
        }
        if(Surface_type[itype] == finite_cyl_3D || Surface_type[itype]==cyl_periodic_3D) {
           Xtest_reflect_TF[Link[iwall]][Orientation[itype]] = FALSE;
        }
     }
  }

  /* if we have reflective or periodic boundaries, we need
     to assess contributions of images.
     So first locate all of the images
     of the surfaces in the system, and expand the loops
     over the walls to include these images. */

    
  image_pos = (double **) array_alloc (2, Nwall*POW_INT(3,Ndim), 
                                          Ndim, sizeof(double));
  imagetype = (int *) array_alloc (1, Nwall*POW_INT(3,Ndim), sizeof(double));

  image_old = 0;
  for (iwall=0; iwall<Nwall; iwall++){
     for (idim=0; idim<Ndim; idim++) {
        pos[idim] = WallPos[idim][iwall];
        image_pos[image_old][idim] = pos[idim];
     }
     image = image_old + 1;

     if (Xtest_reflect_TF[Link[iwall]][0])
     find_wall_images(0,&image,image_pos,pos);

     if (Ndim > 1){
        image_x = image-image_old;
        for (i=image_old; i<image_old+image_x; i++){
           for (idim=0; idim<Ndim; idim++)
              pos[idim] = image_pos[i][idim];
              if (Xtest_reflect_TF[Link[iwall]][1])
              find_wall_images(1, &image,image_pos,pos);
        }
     }
     if (Ndim == 3){
        image_xy = image-image_old;
        for (i=image_old; i<image_old+image_xy; i++){
           for (idim=0; idim<Ndim; idim++)
              pos[idim] = image_pos[i][idim];

              if (Xtest_reflect_TF[Link[iwall]][2])
              find_wall_images(2,&image,image_pos,pos);
        }

     }
     for (i=image_old; i<image; i++) imagetype[i]=WallType[iwall];
     image_old=image;
  }
  nwall_max = image;

  icount=0;
  for (iwall=0; iwall<nwall_max; iwall++){
  itype = imagetype[iwall];
  real_wall = iwall/(nwall_max/Nwall);
  switch(Surface_type[itype])
  {
     case smooth_planar_wall:
       els_planar(iwall,real_wall,itype, L_wall,x_min,
                  nelems_w_per_w, elems_w_per_w,
                  fast_fill_elem_TF,el_type,
                  image_pos);
       break;

     case finite_planar_wall:
       els_finite_planar(iwall,real_wall,itype, L_wall,x_min,
                         nelems_w_per_w, elems_w_per_w,
                         fast_fill_elem_TF,el_type,
                         image_pos);
       break; 

     case colloids_cyl_sphere:
     case point_surface:
       els_spheres(iwall,real_wall,itype,L_wall,x_min,
                   nelems_w_per_w, elems_w_per_w,
                   fast_fill_elem_TF,el_type,
                   image_pos);
       break;

     case finite_cyl_3D:
       els_cyls_3D(iwall,real_wall,itype,L_wall,x_min,
                   nelems_w_per_w, elems_w_per_w,
                   fast_fill_elem_TF,el_type,
                   image_pos);
       break;

     case cyl_periodic_3D:
       els_cyls_cos_3D(iwall,real_wall,itype,L_wall,x_min,
                     nelems_w_per_w, elems_w_per_w,
                     fast_fill_elem_TF,el_type,
                     image_pos);
       break;


     case atomic_centers:
       els_atomic_centers(iwall,real_wall,itype,L_wall,x_min,
                   nelems_w_per_w, elems_w_per_w,
                   fast_fill_elem_TF,el_type,
                   image_pos);
       break;
     
     case cyl2D_sphere3D_pore: 
       els_cyl_pores(iwall,real_wall,itype,L_wall,x_min,
                     nelems_w_per_w, elems_w_per_w,
                     fast_fill_elem_TF,el_type,
                     image_pos);
       break;

     case cyl3D_slit2D_pore:            
       if (Ndim == 2)
          els_slit_pore_2D(iwall,real_wall,itype,L_wall,x_min,
                          nelems_w_per_w,elems_w_per_w,
                          fast_fill_elem_TF,el_type,
                          image_pos);

       else if (Ndim == 3)
          els_cyl_pore_3D(iwall,real_wall,itype,L_wall,x_min,
                          nelems_w_per_w,elems_w_per_w,
                          fast_fill_elem_TF,el_type,
                          image_pos);
       break;

     case tapered_pore:            
       if (Ndim == 2)
          els_cone_pore_2D(iwall,real_wall,itype,L_wall,x_min,
                          nelems_w_per_w,elems_w_per_w,
                          fast_fill_elem_TF,el_type,
                          image_pos);
       else if (Ndim == 3)
          els_cone_pore_3D(iwall,real_wall,itype,L_wall,x_min,
                          nelems_w_per_w,elems_w_per_w,
                          fast_fill_elem_TF,el_type,
                          image_pos);
       break;

     default:		/* No surfaces */
       fprintf (fp2,"ERROR:the surface type chosen is not available\n");
       fprintf (fp2,"Surface_type:%d\n",Surface_type[itype]);
       fclose(fp2);
       break;
  }
     if (icount==100 && Iwrite==VERBOSE) {
        printf("Proc: %d in setup surfaces:: iwall=%d of nwall_max=%d Nnodes_wall_box: %d \n",
                                             Proc,iwall,nwall_max,Nnodes_wall_box);
                icount=0;}
     else if (iwall==nwall_max-1 && Iwrite==VERBOSE && Proc==0){
        printf("Proc: %d is done with setup surfaces:: iwall=%d of nwall_max=%d Nnodes_wall_box=%d\n",
                                              Proc,iwall,nwall_max,Nnodes_wall_box);
     }  
     else icount++;
  } /* end of loop over walls */

  /* now finish setting up fluid element list and zones */
  for (ilist=0; ilist<Nlists_HW; ilist++){
     nelems_f[ilist] = 0;
     for (iel_box=0; iel_box<Nelements_box; iel_box++){
        if (L_wall[ilist][iel_box] == FALSE) {
            elems_f[ilist][nelems_f[ilist]] = iel_box;
            nelems_f[ilist]++;
        }
        if (ilist == Nlists_HW-1){
           for (izone = Nzone-2; izone>=0; izone--){
             if (x_min[iel_box] <= Rmax_zone[izone]) 
                               elem_zones[iel_box] = izone;
           }
        }
      }
  }
  safe_free((void *) &imagetype);
  safe_free((void *) &image_pos);
  safe_free((void *) &L_wall);
  safe_free((void *) &x_min);
  
}  /* end of surface element setup*/
/****************************************************************************/
/*find_wall_images: here we take the surface positions from the
 *    central cell and locate all the corresponding images due
 *    to reflective and periodic boundaries.
 */
void find_wall_images(int idim, 
		      int *image, double **image_pos,
		      double *pos)
{
  int iside,jdim;
  double sign,shift=0,node_image[3];

  for (iside = 0; iside < 2; iside++){
     if (iside == 0) sign=-1.0; else sign=1.0;
  
     if (Type_bc[idim][iside] == PERIODIC ||
         Type_bc[idim][iside] == REFLECT) {

        if (Type_bc[idim][iside] == PERIODIC)
           shift = Size_x[idim];
        else if (Type_bc[idim][iside] == REFLECT)
           shift = 2.0*fabs(sign*0.5*Size_x[idim] - pos[idim]);

        if (fabs(shift) > 0.000000001) {
           for (jdim=0; jdim<Ndim; jdim++){
              if (jdim==idim) node_image[jdim] = pos[jdim] + sign*shift;
              else            node_image[jdim] = pos[jdim];
           }
;
           for (jdim=0; jdim<Ndim; jdim++)
              image_pos[*image][jdim] = node_image[jdim];
           (*image)++;
        }

     }  /* end of check for b.c. type */
   }  /* end of left side / right side test */

   return;
}
/****************************************************************************/
/* flag_wall_el: set all appropriate flags if a given element is in a wall */
void flag_wall_el(int inode,int ilist,int iwall,int iel_box, int **L_wall,
             int **nelems_w_per_w, int ***elems_w_per_w, int ***el_type)
{
   int ijk[3],idim,nadd,i,new_wall,new_entry,diff_wall;
   int imax[3],ix[3],iw,ijk_box[3],inode_box,inode_tmp_box,ijk_tmp_box[3],index;
   imax[0]=imax[1]=imax[2]=1;

   new_wall=TRUE;
   for (i=0; i<Nwall_owners[ilist][iel_box]; i++)
      if (Link[iwall]==Link[i]) new_wall=FALSE;
   if (new_wall) Wall_owners[ilist][iel_box][Nwall_owners[ilist][iel_box]++] = Link[iwall];

   if (iwall < Nwall) {
      elems_w_per_w[ilist][nelems_w_per_w[ilist][iwall]][iwall] = iel_box;
      nelems_w_per_w[ilist][iwall]++;
      Wall_elems[ilist][iel_box] = iwall;
      L_wall[ilist][iel_box] = TRUE;
      el_type[iwall][ilist][iel_box] = WALL_EL;


      /* setup wall dielectric constants */
      if (ilist == Nlists_HW-1 && Ipot_ff_c == COULOMB){
          if (Type_dielec != DIELEC_CONST)
                  Dielec[iel_box] = Dielec_wall[WallType[iwall]];
      }

      /* flag cases where a given wall touches the domain boundary - need
         to do some special things if we touch a reflective boundary */
      node_to_ijk(inode,ijk);

      for (idim=0; idim<Ndim; idim++){
      if (ijk[idim] == 0 && !Touch_domain_boundary[iwall][ilist][idim][0]) 
                       Touch_domain_boundary[iwall][ilist][idim][0] = TRUE;
          if (Type_bc[idim][1]==PERIODIC) nadd=0;
          else                            nadd=1;
          if (ijk[idim]+nadd == Nodes_x[idim]-1 && 
                       !Touch_domain_boundary[iwall][ilist][idim][1]) 
                       Touch_domain_boundary[iwall][ilist][idim][1] = TRUE;
      }

      /* make a list of nodes that are hit by wall elements */

      inode_box = element_box_to_node_box(iel_box);
      node_box_to_ijk_box(inode_box,ijk_box);
      for (idim=0; idim<Ndim; idim++) imax[idim]=2;

      for (ix[2]=0; ix[2]<imax[2]; ix[2]++){
        for (ix[1]=0; ix[1]<imax[1]; ix[1]++){
          for (ix[0]=0; ix[0]<imax[0]; ix[0]++){

             for (idim=0; idim<Ndim; idim++){
                 ijk_tmp_box[idim]=ijk_box[idim]+ix[idim]; 
                 if (ijk_tmp_box[idim]==Nodes_x[idim] && Type_bc[idim][1]==PERIODIC && Pflag[idim]){
                      ijk_tmp_box[idim]=0; 
                 }
             }
             inode_tmp_box=ijk_box_to_node_box(ijk_tmp_box);

             if (Index_wall_nodes[inode_tmp_box]==-1) {  /*node not yet in list */
                  Nodes_wall_box[Nnodes_wall_box]=inode_tmp_box;
                  Index_wall_nodes[inode_tmp_box] = Nnodes_wall_box;
                  Nnodes_wall_box++;
             }
             index = Index_wall_nodes[inode_tmp_box];

          /*    determine if this wall is already in the list for this node - if not add it
                also index the list number for use later. */

             new_entry=TRUE; 
             if (Nnodes_wall_box > 0){
                for (iw=0; iw< Nwall_touch_node[index]; iw++){
                   if (Wall_touch_node[index][iw]==iwall && List_wall_node[index][iw]==ilist) new_entry=FALSE;
                }
             }
             if (new_entry) {
                  Wall_touch_node[index][Nwall_touch_node[index]]=iwall;
                  List_wall_node[index][Nwall_touch_node[index]++]=ilist;
             }
          }
        }
      }
   }
   return;
}
/****************************************************************************/
/* els_spheres: Assign each element to either the fluid 
                                      or the fixed surface spheres */

void els_spheres(int iwall, int real_wall, int itype, int **L_wall, double *x_min,
                            int **nelems_w_per_w, int ***elems_w_per_w,
                            int *fast_fill_elem_TF,int ***el_type,
                            double **image_pos)
{
   int iel, iel_box,inode, idim, ilist,L1el_charge;
   double xtest[3], node_pos[3];
   double r12_sq_sum, r12, dx, radius,sten,delr,roff=0.00000000001;
   double vecB[2],angle_block,angle,roughness,shift,cos_theta,angle_deg;
   int iblock,rough_block[2],idim_permute,angle_test;

   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];


   for (ilist=0; ilist<Nlists_HW; ilist++){

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){

         iel=el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);

         switch(Ndim)   /* xtest = position at center of element */
         {
            case 3:  xtest[2] = node_pos[2] + 0.5*Esize_x[2];
            case 2:  xtest[1] = node_pos[1] + 0.5*Esize_x[1];
            default: xtest[0] = node_pos[0] + 0.5*Esize_x[0];   
         }

                                                                /* compute surface for a cylindrical wedge */
         angle_test=TRUE;
         if (Ndim==2 && (fabs(WallParam_2[itype]-0.0)>=1.e-12 || fabs(WallParam_3[itype]-360.0)>1.e-12)){
               vecB[0]=xtest[0]-WallPos[0][iwall]; vecB[1]=xtest[1]-WallPos[1][iwall]; idim_permute=1;
               cos_theta=vecB[0]/(sqrt(vecB[0]*vecB[0]+vecB[1]*vecB[1]));
               angle = acos(cos_theta);
               if (xtest[idim_permute] <WallPos[idim_permute][iwall]) shift=2.0*(PI-angle);
               else shift=0.0;
               angle+=shift;
               angle_deg=180.*angle/PI;
               if(WallParam_3[itype]>WallParam_2[itype]){
                      if(angle_deg <WallParam_2[itype] || angle_deg > WallParam_3[itype]) angle_test=FALSE;
               }
               else if(WallParam_3[itype]<WallParam_2[itype]){
                      if(angle_deg <WallParam_2[itype] && angle_deg > WallParam_3[itype]) angle_test=FALSE;
               }
         }

         if (Lrough_surf[itype] && Ndim==2){
               for (iblock=0;iblock<2;iblock++) rough_block[iblock]=0;

               rough_block[1] = 0;

               if (Rough_length[itype] > WallParam[itype]){
                  printf("ERROR IN LENGTH SCALE FOR ROUGHNESS...Rough length larger than 1/4 of cylinder\n");
                  exit(-1);
               }
               angle_block = asin(Rough_length[itype]/WallParam[itype]);
               vecB[0]=xtest[0];
               vecB[1]=xtest[1];
               idim_permute=1;
               cos_theta=vecB[0]/(sqrt(vecB[0]*vecB[0]+vecB[1]*vecB[1]));
               angle = acos(cos_theta);
               if (vecB[1] <WallPos[idim_permute][iwall]) shift=2.0*(PI-angle);
               else shift=0.0;
               angle+=shift;
               rough_block[0]= (int)(angle/angle_block);

               if (rough_block[0] >= MAX_ROUGH_BLOCK || rough_block[1]>=MAX_ROUGH_BLOCK) {
                    printf("ERROR with rough cylinder - number of rough patches exceeds maximum of MAX_ROUGH_BLOCK=%d\n",
                            MAX_ROUGH_BLOCK);
                    exit(-1);
               }
               roughness = Rough_precalc[itype][rough_block[0]][rough_block[1]];
         }
         else roughness=0.0;

   
         radius = WallParam[itype]+roughness;
         if(Lhard_surf && ilist != Nlists_HW-1 && (Type_func != NONE || Type_poly !=NONE) )  
                     radius += 0.5*Sigma_ff[ilist][ilist];

         radius -= roff; /* adjust for round-off errors */

         r12_sq_sum = 0.0;
         for (idim = 0; idim < Ndim; idim++) {
            dx =  xtest[idim] -  image_pos[iwall][idim];
            r12_sq_sum = r12_sq_sum + dx*dx;
         }
         r12 = sqrt(r12_sq_sum);

         /* for point charge surfaces only allow one element per surface*/
         L1el_charge = FALSE;
         if (Ipot_ff_c==COULOMB && Surface_type[WallType[iwall]]==point_surface &&
                nelems_w_per_w[ilist][iwall]==1) L1el_charge=TRUE;

         if (r12 <= radius && angle_test && !L1el_charge) 
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         else{
            /* in fluid ... flag slow fill regions */
              if (Lhard_surf && r12 <= radius + sten)
                              fast_fill_elem_TF[iel_box] = FALSE; 

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE &&
                        r12 <= radius + Dielec_X) 
                              Dielec[iel_box] = Dielec_pore;
              }
            
         }


         if (ilist == Nlists_HW-1){
             delr = r12-radius;
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }
      }     /* end of loop over elements       */
   }        /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_cyl_3D: Assign each element to either the fluid 
                 or the fixed surface cylinders in 3D ...
                 this routine is set up mainly to be able
                 to test the 3D code against 2D cylinders */

void els_cyls_3D(int iwall, int real_wall, int itype, int **L_wall, double *x_min,
                            int **nelems_w_per_w, int ***elems_w_per_w,
                            int *fast_fill_elem_TF,int ***el_type,
                            double **image_pos)
{
   int iel, iel_box,inode, idim, ilist;
   double xtest[3], node_pos[3],length;
   double r12_sq_sum, r12, dx, radius,sten,delr,roff=0.00000000001;
   double vecB[2],angle_block,angle,roughness,shift,cos_theta,angle_deg;
   int iblock,rough_block[2],idim_permute,angle_test;

   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];

   for (ilist=0; ilist<Nlists_HW; ilist++){

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){

         iel=el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);

         switch(Ndim)   /* xtest = position at center of element */
         {
            case 3:  xtest[2] = node_pos[2] + 0.5*Esize_x[2];
            case 2:  xtest[1] = node_pos[1] + 0.5*Esize_x[1];
            default: xtest[0] = node_pos[0] + 0.5*Esize_x[0];   
         }

         angle_test=TRUE;
         if ((fabs(WallParam_3[itype]-0.0)>=1.e-12 || fabs(WallParam_4[itype]-360.0)>1.e-12)){
               if (Orientation[itype]==0){ vecB[0]=xtest[1]-WallPos[1][iwall];
                                           vecB[1]=xtest[2]-WallPos[2][iwall];
                                           idim_permute=2;}
               else if (Orientation[itype]==1){ vecB[0]=xtest[2]-WallPos[2][iwall];
                                                vecB[1]=xtest[0]-WallPos[0][iwall];
                                                idim_permute=0;}
               else if (Orientation[itype]==2){ vecB[0]=xtest[0]-WallPos[0][iwall]; 
                                                vecB[1]=xtest[1]-WallPos[1][iwall]; 
                                                idim_permute=1;}
               cos_theta=vecB[0]/(sqrt(vecB[0]*vecB[0]+vecB[1]*vecB[1]));
               angle = acos(cos_theta);
               if (xtest[idim_permute] <WallPos[idim_permute][iwall]) shift=2.0*(PI-angle);
               else shift=0.0;
               angle+=shift;
               angle_deg=180.*angle/PI;
               if(WallParam_4[itype]>WallParam_3[itype]){
                      if(angle_deg <WallParam_3[itype] || angle_deg > WallParam_4[itype]) angle_test=FALSE;
               }
               else if(WallParam_4[itype]<WallParam_3[itype]){
                      if(angle_deg <WallParam_3[itype] && angle_deg > WallParam_4[itype]) angle_test=FALSE;
               }
         }

         if (Lrough_surf[itype]){
               for (iblock=0;iblock<2;iblock++) rough_block[iblock]=0;

               rough_block[1] = (int) ((xtest[Orientation[itype]]+0.5*Size_x[Orientation[itype]])/Rough_length[itype]);

               if (Rough_length[itype] > WallParam[itype]){
                  printf("ERROR IN LENGTH SCALE FOR ROUGHNESS...Rough length larger than 1/4 of cylinder\n");
                  exit(-1);
               }
               angle_block = asin(Rough_length[itype]/WallParam[itype]);
               if (Orientation[itype]==0){ vecB[0]=xtest[1];vecB[1]=xtest[2];idim_permute=2;}
               else if (Orientation[itype]==1){ vecB[0]=xtest[2];vecB[1]=xtest[0];idim_permute=0;}
               else if (Orientation[itype]==2){ vecB[0]=xtest[0];vecB[1]=xtest[1];idim_permute=1;}
               cos_theta=vecB[0]/(sqrt(vecB[0]*vecB[0]+vecB[1]*vecB[1]));
               angle = acos(cos_theta);
               if (vecB[1]<WallPos[idim_permute][iwall]) shift=2.0*(PI-angle);
               else shift=0.0;
               angle+=shift;
               rough_block[0]= (int)(angle/angle_block);
               
               if (rough_block[0] >= MAX_ROUGH_BLOCK || rough_block[1]>=MAX_ROUGH_BLOCK) {
                    printf("ERROR with rough cylinder - number of rough patches exceeds maximum of MAX_ROUGH_BLOCK=%d\n",
                            MAX_ROUGH_BLOCK);
                    exit(-1);
               }
               roughness = Rough_precalc[itype][rough_block[0]][rough_block[1]];
         }
         else roughness=0.0;
   
         radius = WallParam[itype]+roughness;
         length = WallParam_2[itype];
         if(Lhard_surf && ilist != Nlists_HW-1)  
                     radius += 0.5*Sigma_wf[ilist][itype];

         radius -= roff; /* adjust for round-off errors */

         r12_sq_sum = 0.0;
         for (idim = 0; idim < Ndim; idim++) {
            dx =  xtest[idim] -  image_pos[iwall][idim];
            if( idim != Orientation[itype])
                   r12_sq_sum = r12_sq_sum + dx*dx;
         }
         r12 = sqrt(r12_sq_sum);

         if (r12 <= radius && 
             fabs(xtest[Orientation[itype]]-image_pos[iwall][Orientation[itype]])<=length && angle_test) 
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         else{
            /* in fluid ... flag slow fill regions */
              if (Lhard_surf && r12 <= radius + sten)
                              fast_fill_elem_TF[iel_box] = FALSE; 

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE &&
                        r12 <= radius + Dielec_X) 
                              Dielec[iel_box] = Dielec_pore;
              }
            
         }


         if (ilist == Nlists_HW-1){
             delr = r12-radius;
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }
      }     /* end of loop over elements       */
   }        /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_cyl_periodic_3D: Assign each element to either the fluid 
                 or the fixed wavy surface cylinders in 3D ...*/

void els_cyls_cos_3D(int iwall, int real_wall, int itype, int **L_wall, double *x_min,
                            int **nelems_w_per_w, int ***elems_w_per_w,
                            int *fast_fill_elem_TF,int ***el_type,
                            double **image_pos)
{
   int iel, iel_box,inode, idim, ilist;
   double xtest[3], node_pos[3];
   double r12_sq_sum, r12, dx, radius,sten,delr,roff=0.00000000001;

   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];

   for (ilist=0; ilist<Nlists_HW; ilist++){

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){

         iel=el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);

         switch(Ndim)   /* xtest = position at center of element */
         {
            case 3:  xtest[2] = node_pos[2] + 0.5*Esize_x[2];
            case 2:  xtest[1] = node_pos[1] + 0.5*Esize_x[1];
            default: xtest[0] = node_pos[0] + 0.5*Esize_x[0];   
         }
  
         idim = Orientation[itype]; 
         radius = WallParam[itype] + WallParam_2[itype]*
                   cos(2*PI*(xtest[idim]-WallPos[idim][iwall])/WallParam_3[itype]);

         if(Lhard_surf && ilist != Nlists_HW-1)  
                     radius += 0.5*Sigma_wf[ilist][itype];

         radius -= roff; /* adjust for round-off errors */

         r12_sq_sum = 0.0;
         for (idim = 0; idim < Ndim; idim++) {
            dx =  xtest[idim] -  image_pos[iwall][idim];
            if( idim != Orientation[itype])
                   r12_sq_sum = r12_sq_sum + dx*dx;
         }
         r12 = sqrt(r12_sq_sum);

         if (r12 <= radius ) 
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         else{
            /* in fluid ... flag slow fill regions */
              if (Lhard_surf && r12 <= radius + sten)
                              fast_fill_elem_TF[iel_box] = FALSE; 

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE &&
                        r12 <= radius + Dielec_X) 
                              Dielec[iel_box] = Dielec_pore;
              }
            
         }


         if (ilist == Nlists_HW-1){
             delr = r12-radius;
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }
      }     /* end of loop over elements       */
   }        /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_atomic_centers: Assign each element to either the fluid 
                                or the fixed atoms */

void els_atomic_centers(int iwall, int real_wall,int itype, int **L_wall, double *x_min,
                            int **nelems_w_per_w, int ***elems_w_per_w,
                            int *fast_fill_elem_TF,int ***el_type,
                            double **image_pos)
{
   int iel, iel_box,inode, idim, ilist;
   double xtest[3], node_pos[3];
   double r12_sq_sum, r12, dx, radius,sten,delr,rsq;
   double roff=0.00000000001;


   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];

   for (ilist=0; ilist<Nlists_HW; ilist++){

       if (Lhard_surf){ 
          radius = 0.5*Sigma_ww[itype][itype];
          if (ilist != Nlists_HW-1 && (Type_func !=NONE || Type_poly !=NONE) )  
                   radius += 0.5*Sigma_ff[ilist][ilist];
       }
       else{
         if (ilist==0){
            rsq = 0.0;
            for (idim=0; idim<Ndim; idim++) rsq += Esize_x[idim]*Esize_x[idim];
            radius = sqrt(rsq);
         }
         else{
            /* this will be used for the proper computation of surface area for LJ atomic systems */
            radius = 0.5*Sigma_ww[itype][itype]+0.5*Sigma_ff[ilist][ilist];
         }
       }
       radius -= roff; /* adjust for round-off errors */

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){

         iel=el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);


         switch(Ndim)   /* xtest = position at center of element */
         {
            case 3:  xtest[2] = node_pos[2] + 0.5*Esize_x[2];
            case 2:  xtest[1] = node_pos[1] + 0.5*Esize_x[1];
            default: xtest[0] = node_pos[0] + 0.5*Esize_x[0];   
         }
 
         r12_sq_sum = 0.0;
         for (idim = 0; idim < Ndim; idim++) {
            dx =  xtest[idim] -  image_pos[iwall][idim];
            r12_sq_sum = r12_sq_sum + dx*dx;
         }
         r12 = sqrt(r12_sq_sum);

         if (r12 <= radius ) {
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         }
         else{
             /* in fluid ... flag slow fill regions */
              if (Lhard_surf && r12 <= radius + sten)
                              fast_fill_elem_TF[iel_box] = FALSE; 

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE &&
                        r12 <= radius + Dielec_X) 
                              Dielec[iel_box] = Dielec_pore;
              }
         }


         if (ilist == Nlists_HW-1){
             delr = r12-radius;
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }
      }     /* end of loop over elements       */
   }        /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_slit_pore_2D: Assign each element to either the fluid 
                     or the pore walls.  This routine is for
                     2D calculations where the Orientation specifies the
                     axis down the length of the slit pore.  */

void els_slit_pore_2D(int iwall, int real_wall, int itype, int **L_wall, double *x_min,
		      int **nelems_w_per_w, int ***elems_w_per_w,
		      int *fast_fill_elem_TF,int ***el_type,
		      double **image_pos)
{
   int iel, iel_box, inode, idim, ilist,dim[3],loc_inode;
   double xtest[3];
   double node_pos[3];
   double r12, radius,sten,x12,length,delr,x12new,r12new;
   double roff=0.00000000001;

   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];

   for (ilist = 0; ilist<Nlists_HW; ilist++){

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){
         
         iel = el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);
   
         switch(Orientation[itype])  /* xtest = center of element */
         {
            case 1:
               xtest[0] = node_pos[0] + 0.5*Esize_x[0]; dim[0] = 0;
               xtest[1] = node_pos[1] + 0.5*Esize_x[1]; dim[1] = 1;
               break;
            case 0:
               xtest[0] = node_pos[1] + 0.5*Esize_x[1]; dim[0] = 1;  
               xtest[1] = node_pos[0] + 0.5*Esize_x[0]; dim[1] = 0;
               break;
         }

         radius = WallParam[itype];
         if (Lhard_surf && ilist != Nlists_HW-1)
                       radius -= 0.5*Sigma_wf[ilist][itype];
         radius += roff;
         r12 = fabs(xtest[0] -  image_pos[iwall][dim[0]]);

         length = 2.0*WallParam_2[itype] - 2.0*roff;
         x12 = fabs(image_pos[iwall][dim[1]] - xtest[1]);

         if (r12 >= radius && x12 <= 0.5*length) {
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         }
         else{
             /* in fluid ... flag slow fill regions */
             if (Lhard_surf && r12 >= radius-sten)
                          fast_fill_elem_TF[iel_box] = FALSE; 

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE &&
                        r12 >= radius - Dielec_X) 
                              Dielec[iel_box] = Dielec_pore;
              }
         }

             /* in fluid ... store the value of x12 if we are doing
                a stepped 9-3 surface */
              if (Ipot_wf_n[iwall]==VEXT_1D_XMIN){
                 switch(Orientation[itype])  /* xtest = llb node of element */
                 {
                   case 1:  xtest[0] = node_pos[0]; xtest[1] = node_pos[1]; break;
                   case 0:  xtest[0] = node_pos[1];  xtest[1] = node_pos[0]; break;
                 }
                 r12new = fabs(xtest[0] -  image_pos[iwall][dim[0]]);
                 x12new = fabs(image_pos[iwall][dim[1]] - xtest[1]);

                 if (x12new <=0.5*length +2.0*roff){
                    loc_inode = B2L_node[node_to_node_box(inode)];
                    if (loc_inode>0){
                        X_wall[loc_inode][real_wall]=radius-r12new;
                        X_wall2[loc_inode][real_wall]=radius+r12new;
                    }
                 }
              }    /* END of VEXT_1D_XMIN case */

         if (ilist == Nlists_HW-1) {
             delr = WallParam[itype] - r12;
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }

      }  /* end of loop over elements */
   }     /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_cyl_pore_3D: Assign each element to either the fluid 
                     or the cylindrical pore walls.  This routine is for
                     3D calculations where the Orientation specifies the
                     axis down the length of the cylinder. */ 

void els_cyl_pore_3D(int iwall, int real_wall, int itype, int **L_wall, double *x_min,
                              int **nelems_w_per_w, int ***elems_w_per_w,
                              int *fast_fill_elem_TF,int ***el_type,
                              double **image_pos)
{
   int iel, iel_box, inode, idim, ilist,dim[3];
   double xtest[3];
   double node_pos[3];
   double r12_sq_sum, r12, dx,radius,sten,x12,length,delr;
   double roff=0.00000000001;

   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];

   for (ilist = 0; ilist<Nlists_HW; ilist++){

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){
         
         iel = el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);
   

         switch(Orientation[itype])  /* xtest = center of element */
         {
            case 2:
               xtest[0] = node_pos[0] + 0.5*Esize_x[0]; dim[0] = 0;
               xtest[1] = node_pos[1] + 0.5*Esize_x[1]; dim[1] = 1;
               xtest[2] = node_pos[2] + 0.5*Esize_x[2]; dim[2] = 2;
               break;
            case 1:
               xtest[0] = node_pos[2] + 0.5*Esize_x[2]; dim[0] = 2;
               xtest[1] = node_pos[0] + 0.5*Esize_x[0]; dim[1] = 0;
               xtest[2] = node_pos[1] + 0.5*Esize_x[1]; dim[2] = 1;
               break;
            case 0:
               xtest[0] = node_pos[1] + 0.5*Esize_x[1]; dim[0] = 1;  
               xtest[1] = node_pos[2] + 0.5*Esize_x[2]; dim[1] = 2;
               xtest[2] = node_pos[0] + 0.5*Esize_x[0]; dim[2] = 0;
               break;
         }


         radius = WallParam[itype];
         if (Lhard_surf && ilist != Nlists_HW-1 )  
              radius -= 0.5*Sigma_wf[ilist][itype];
         radius += roff;

         r12_sq_sum = 0.0;
         for (idim = 0; idim < Ndim-1; idim++) {
            dx =  xtest[idim] -  image_pos[iwall][dim[idim]];
            r12_sq_sum = r12_sq_sum + dx*dx;
         }
         r12 = sqrt(r12_sq_sum);

         length = 2.0*WallParam_2[itype] - 2.0*roff;
         x12 = fabs(image_pos[iwall][dim[2]] - xtest[2]);

         if (r12 >= radius && x12 <= 0.5*length) 
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         else{
            /* in fluid ... flag slow fill regions */
             if (Lhard_surf && r12 >= radius-sten)
                           fast_fill_elem_TF[iel_box] = FALSE; 

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE &&
                        r12 >= radius - Dielec_X) 
                              Dielec[iel_box] = Dielec_pore;
              }
         }


         if (ilist == Nlists_HW-1) {
             delr = x12 - WallParam_2[itype];
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }

      }  /* end of loop over elements */
   }     /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_cone_pore_2D: Assign each element to either the fluid 
                     or the cylindrical pore walls.  This routine is for
                     3D calculations where the Orientation specifies the
                     axis down the length of the cylinder.  It is assumed
                     that the ends of the box are periodic.  */

void els_cone_pore_2D(int iwall, int real_wall, int itype, int **L_wall, double *x_min,
		      int **nelems_w_per_w, int ***elems_w_per_w,
		      int *fast_fill_elem_TF,int ***el_type,
		      double **image_pos)
{
   int iel, iel_box, inode, idim, ilist,dim[3],loc_inode,operp;
   double xtest[3];
   double node_pos[3];
   double r12, delr,
          radius,radius1,radius2,sten,
          x12,length,r12new,x12new,radiusnew;
   double roff=0.00000000001;

   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];

   for (ilist = 0; ilist<Nlists_HW; ilist++){

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){
         
         iel = el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);
   
         switch(Orientation[itype])  /* xtest = center of element */
         {
            case 1:
               xtest[0] = node_pos[0] + 0.5*Esize_x[0]; dim[0] = 0;
               xtest[1] = node_pos[1] + 0.5*Esize_x[1]; dim[1] = 1;
               operp=0;
               break;
            case 0:
               xtest[0] = node_pos[1] + 0.5*Esize_x[1]; dim[0] = 1;  
               xtest[1] = node_pos[0] + 0.5*Esize_x[0]; dim[1] = 0;
               operp=1;
               break;
         }

         length = 2.0*WallParam_3[itype] - 2.0*roff;
         x12 = xtest[1] - image_pos[iwall][dim[1]];

         radius1 = WallParam[itype];
         radius2 = WallParam_2[itype];
         if (Lhard_surf && ilist != Nlists_HW-1) {
            radius1 -= 0.5*Sigma_wf[ilist][itype];
            radius2 -= 0.5*Sigma_wf[ilist][itype];
         }
         radius = radius2 + (radius2-radius1)*(x12-0.5*length)/length;
         radius += roff;

         r12 = fabs(xtest[0] - image_pos[iwall][dim[0]]);
         x12 = fabs(x12);

         if (r12 >= radius && x12 <= 0.5*length) 
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         else{
             /* in fluid ... flag slow fill regions */
             if (Lhard_surf && r12 >= radius - sten)
                             fast_fill_elem_TF[iel_box] = FALSE;  

             /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE && 
                      r12 >= radius - Dielec_X) Dielec[iel_box] = Dielec_pore;
              }
         }

             /* in fluid ... store the value of x12 if we are doing
                a stepped 9-3 surface */
              if (Ipot_wf_n[iwall]==VEXT_1D_XMIN){
                 switch(Orientation[itype])  /* xtest = llb node of element */
                 {
                   case 1:  xtest[0] = node_pos[0]; xtest[1] = node_pos[1]; break;
                   case 0:  xtest[0] = node_pos[1];  xtest[1] = node_pos[0]; break;
                 }
                 r12new = fabs(xtest[0] -  image_pos[iwall][dim[0]]);
                 x12new = xtest[1] - image_pos[iwall][dim[1]];
                 radiusnew = radius2 + (radius2-radius1)*(x12new-0.5*length)/length;
                 radiusnew += roff;
                 x12new = fabs(x12new);

                 if (x12new <=0.5*length +2.*roff){
                    loc_inode = B2L_node[node_to_node_box(inode)];
                    if (loc_inode>0){
                        X_wall[loc_inode][real_wall]=radiusnew-r12new;
                        X_wall2[loc_inode][real_wall]=radiusnew+r12new;
                    }
                 }
              }    /* END of VEXT_1D_XMIN case */


         if (ilist == Nlists_HW-1) {
             delr = radius - r12;
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }

      }  /* end of loop over elements */
   }     /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_cone_pore_3D: Assign each element to either the fluid 
                     or the cylindrical pore walls.  This routine is for
                     3D calculations where the Orientation specifies the
                     axis down the length of the cylinder.  It is assumed
                     that the ends of the box are periodic.  */

void els_cone_pore_3D(int iwall, int real_wall, int itype, int **L_wall, double *x_min,
                              int **nelems_w_per_w, int ***elems_w_per_w,
                              int *fast_fill_elem_TF,int ***el_type,
                              double **image_pos)
{
   int iel, iel_box, inode, idim, ilist,dim[3];
   double xtest[3];
   double node_pos[3];
   double r12_sq_sum, r12, delr,dx, 
          radius,radius1,radius2,sten,
          x12,length;
   double roff=0.00000000001;

   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];

   for (ilist = 0; ilist<Nlists_HW; ilist++){

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){
         
         iel = el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);
   
         switch(Orientation[itype])  /* xtest = center of element */
         {
            case 2:
               xtest[0] = node_pos[0] + 0.5*Esize_x[0]; dim[0] = 0;
               xtest[1] = node_pos[1] + 0.5*Esize_x[1]; dim[1] = 1;
               xtest[2] = node_pos[2] + 0.5*Esize_x[2]; dim[2] = 2;
               break;
            case 1:
               xtest[0] = node_pos[2] + 0.5*Esize_x[2]; dim[0] = 2;
               xtest[1] = node_pos[0] + 0.5*Esize_x[0]; dim[1] = 0;
               xtest[2] = node_pos[1] + 0.5*Esize_x[1]; dim[2] = 1;
               break;
            case 0:
               xtest[0] = node_pos[1] + 0.5*Esize_x[1]; dim[0] = 1;  
               xtest[1] = node_pos[2] + 0.5*Esize_x[2]; dim[1] = 2;
               xtest[2] = node_pos[0] + 0.5*Esize_x[0]; dim[2] = 0;
               break;
         }


         length = 2.0*WallParam_3[itype] - 2.0*roff;
         x12 = xtest[2] - image_pos[iwall][dim[2]];

         radius1 = WallParam[itype];
         radius2 = WallParam_2[itype];
         if (Lhard_surf && ilist != Nlists_HW-1){
            radius1 -= 0.5*Sigma_wf[ilist][itype];
            radius2 -= 0.5*Sigma_wf[ilist][itype];
         }
         radius = radius2 + (radius2-radius1)*(x12-0.5*length)/length;
         radius += roff;

         r12_sq_sum = 0.0;
         for (idim = 0; idim < Ndim-1; idim++) {
            dx =  xtest[idim] -  image_pos[iwall][dim[idim]];
            r12_sq_sum = r12_sq_sum + dx*dx;
         }
         r12 = sqrt(r12_sq_sum);
         x12 = fabs(x12);

         if (r12 >= radius && x12 <= 0.5*length) {
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         }
         else{
             /* in fluid ... flag slow fill regions */
             if (Lhard_surf && r12 >= radius - sten)
                               fast_fill_elem_TF[iel_box] = FALSE;   

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE &&
                        r12 >= radius - Dielec_X) Dielec[iel_box] = Dielec_pore;
              }

         }
         if (ilist == Nlists_HW-1) {
             delr = radius - r12;
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }

      }  /* end of loop over elements */
   }     /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_cyl_pores: Assign each element to either the fluid 
                                      or the cylindrical pore walls */

void els_cyl_pores(int iwall,int real_wall, int itype, int **L_wall, double *x_min,
                              int **nelems_w_per_w, int ***elems_w_per_w,
                              int *fast_fill_elem_TF,int ***el_type,
                              double **image_pos)
{
   int iel, iel_box, inode, idim, ilist,count;
   double xtest[3],node_pos[3];
   double r12_sq_sum, r12, delr,dx,radius,sten;
   double roff=0.00000000001;

   sten=0.0;
   for (idim=0; idim<Ndim; idim++)
       if (Esize_x[idim]*Sten_length_hs[idim] > sten)
          sten = Esize_x[idim]*Sten_length_hs[idim];

   count = 0;

   for (ilist = 0; ilist<Nlists_HW; ilist++){

      for (iel_box = 0; iel_box < Nelements_box; iel_box++){
         
         iel = el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);
   
         switch(Ndim)                 /* xtest = position at center of element */
         {
            case 3:  xtest[2] = node_pos[2] + 0.5*Esize_x[2];
            case 2:  xtest[1] = node_pos[1] + 0.5*Esize_x[1];
            default: xtest[0] = node_pos[0] + 0.5*Esize_x[0];   
         }
   
         radius = WallParam[itype];
         if (Lhard_surf && ilist != Nlists_HW-1) 
             radius -= 0.5*Sigma_ff[ilist][ilist];
         radius += roff;

         r12_sq_sum = 0.0;
         for (idim = 0; idim < Ndim; idim++) {
            dx =  xtest[idim] -  image_pos[iwall][idim];
            r12_sq_sum = r12_sq_sum + dx*dx;
         }
         r12 = sqrt(r12_sq_sum);  

         if (r12 >= radius) 
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         else{
            /* in fluid ... flag slow fill regions */
              if (Lhard_surf && r12 >= radius - sten)
                              fast_fill_elem_TF[iel_box] = FALSE;   

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE &&
                        r12 < radius - Dielec_X) 
                              Dielec[iel_box] = Dielec_pore;
              }
         }


         if (ilist == Nlists_HW-1) {
             delr = WallParam[itype] - r12;
             if (delr > 0.0 ) x_min[iel_box] = AZ_MIN(delr,x_min[iel_box]);
             else             x_min[iel_box] = 0.0;
         }

      }  /* end of loop over elements */
   }     /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_planar: Assign each element to either the fluid 
                                            or the planar surface */

void els_planar(int iwall, int real_wall, int itype, 
                int **L_wall, double *x_min,
                int **nelems_w_per_w, int ***elems_w_per_w,
                int *fast_fill_elem_TF,int ***el_type,
                double **image_pos)
{
  int iel, iel_box,inode, idim, ilist;
  double x, x12, wall_thick,sten,diam,delx;
  double node_pos[3];
  double roff=0.00000000001;
  double roughness;
  int iblock,rough_block[2];

   /*
    * Find the number of fluid and wall elements based on the 
    * center of the element 
    */
   for (ilist=0; ilist<Nlists_HW; ilist++){
      diam = Sigma_ff[ilist][ilist];
   
      for (iel_box = 0; iel_box < Nelements_box; iel_box++){

         iel = el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);

         if (Lrough_surf[itype]){
              if (Ndim==1) roughness=0.0; /* sanity check */
              else{
                  for (iblock=0;iblock<2;iblock++) rough_block[iblock]=0;
                  iblock=0;
                  for (idim=0;idim<Ndim;idim++){
                     if (idim != Orientation[itype]) {
                         rough_block[iblock] = (int) ((node_pos[idim]+0.5*Size_x[idim])/Rough_length[itype]);
                         printf("idim=%d  node_pos[idim]=%9.6f  Rough_length=%9.6f  iblock=%d  rough_block=%d  roughness=%9.6f\n",idim,node_pos[idim],
                                 Rough_length[itype],iblock,rough_block[iblock],Rough_precalc[itype][rough_block[0]][0]);
                         if (rough_block[iblock] >= MAX_ROUGH_BLOCK) {
                             printf("ERROR with rough surfacess - number of rough patches exceeds maximum of MAX_ROUGH_BLOCK=%d\n",
                                      MAX_ROUGH_BLOCK);
                             exit(-1);
                         }
                         iblock++;
                     }
                  }
                  roughness = Rough_precalc[itype][rough_block[0]][rough_block[1]];
              }
         }
         else roughness=0.0;

   
         wall_thick = 2.0*(WallParam[itype]+roughness) - 2.0*roff;
         if (Lhard_surf && ilist !=Nlists_HW-1) wall_thick += diam;

         if      (Orientation[itype] == 0) idim = 0;
         else if (Orientation[itype] == 1) idim = 1;
         else                              idim = 2;

         if (Lhard_surf) 
	   {
	     sten = Esize_x[idim]*Sten_length_hs[idim];
	   }
	 else 
	   {
	     sten = 0;
	   }

         /* get center of element in appropriate dimension */
         x = node_pos[idim] + 0.5*Esize_x[idim];
         x12 = fabs(image_pos[iwall][idim] - x);
         if (x12 <= 0.5*wall_thick){ 
             flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                          nelems_w_per_w, elems_w_per_w,el_type);
         }
         else{
            /* in fluid ... flag fast fill regions */
              if (Lhard_surf && x12 <= (0.5*wall_thick + sten))
                                      fast_fill_elem_TF[iel_box] = FALSE; 

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){

                  if (Type_dielec == DIELEC_WF_PORE &&
                        x12 < 0.5*wall_thick + Dielec_X) 
                              Dielec[iel_box] = Dielec_pore;
              }
         }
      

         if (ilist == Nlists_HW-1){
           delx = x12 - WallParam[itype];
           if (delx > 0.0) x_min[iel_box] = AZ_MIN(delx,x_min[iel_box]);
           else            x_min[iel_box] = 0.0;
         }

      }    /* end of loop over elements */
   }       /* end of loop over boundary lists */
}
/****************************************************************************/
/* els_finite_planar: Assign each element to either the fluid 
                                                   or the planar surface */

void els_finite_planar(int iwall, int real_wall, int itype,
                       int **L_wall, double *x_min,
                       int **nelems_w_per_w, int ***elems_w_per_w,
                       int *fast_fill_elem_TF,int ***el_type,
                       double **image_pos)
{
   int iel, iel_box,inode, idim, 
       ilist,in_wall,in_sten,near_wall,npos;
   double xtest[3], xsave[3],x12[3],length[3],diam,sten[3];
   double node_pos[3],x[3],x_this_wall=0.0,roff=0.00000000001;
   

   for (idim=0; idim<Ndim; idim++)
        sten[idim] = Esize_x[idim]* Sten_length_hs[idim];

   /*
    * Find the number of fluid and wall elements based on the 
    * center of the element 
    */
   for (ilist=0; ilist<Nlists_HW; ilist++){
      diam = Sigma_wf[ilist][itype];
   
      for (iel_box = 0; iel_box < Nelements_box; iel_box++){
   
         iel = el_box_to_el(iel_box);
         inode = element_to_node(iel);
         node_to_position(inode,node_pos);

                       length[0] = 2.0*WallParam[itype];
                       length[1] = 2.0*WallParam_2[itype];
         if (Ndim==3)  length[2] = 2.0*WallParam_3[itype];
         if (Lhard_surf && ilist != Nlists_HW-1){
                       length[0] += diam;
                       length[1] += diam;
          if (Ndim==3) length[2] += diam;
         }
         for (idim=0; idim<Ndim; idim++)
            length[idim] -= 2.0*roff;
         

         in_wall = TRUE;
         in_sten = TRUE;
         near_wall = TRUE;

         for (idim=0; idim<Ndim; idim++){
             x[idim] = node_pos[idim] + 0.5*Esize_x[idim];
             x12[idim] = fabs(image_pos[iwall][idim] - x[idim]);
             if (x12[idim] > 0.5*length[idim])    in_wall = FALSE;
             if (x12[idim] > 0.5*length[idim]+sten[idim]) in_sten = FALSE;
             if (Ipot_ff_c==COULOMB && x12[idim] > 
                  0.5*length[idim]+ Dielec_X)  near_wall = FALSE;
         }

         if (in_wall) 
            flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,
                                   nelems_w_per_w, elems_w_per_w,el_type);
         else{
            /* in fluid ... flag fast fill regions */
              if (in_sten && Lhard_surf)
                      fast_fill_elem_TF[iel_box] = FALSE; /*flag slow fills*/

            /* in fluid ... set up dielectric constants */
              if (ilist == Nlists_HW-1 && Ipot_ff_c==COULOMB){
                  if (Type_dielec == DIELEC_WF_PORE && near_wall)
                        Dielec[iel_box] = Dielec_pore;
              }

         }


         if (ilist == Nlists_HW-1){
             npos = 0;
             for (idim=0; idim<Ndim; idim++){
                 if (idim == 0) 
                       xtest[idim] = x12[idim] - WallParam[itype];
                 else if (idim == 1) 
                       xtest[idim] = x12[idim] - WallParam_2[itype];
                 else  xtest[idim] = x12[idim] - WallParam_3[itype];

                 if (xtest[idim] > 0) xsave[npos++] = xtest[idim];
             }
           
             if (npos == 0) 
	       {
		 x_this_wall = 0.0;
	       }
             else if (npos == 1) 
	       {
		 x_this_wall = xsave[0];
	       }
             else if (npos >= 2)
	       {
                 x_this_wall = AZ_MAX(xsave[0],xsave[1]);
                 if (npos == 3) 
		   {
		     x_this_wall = AZ_MAX(x_this_wall,xsave[2]);
		   }
	       }
             x_min[iel_box] = AZ_MIN(x_min[iel_box],x_this_wall); 
         }


      }    /* end of loop over elements */
   }       /* end of loop over boundary lists */
}
/****************************************************************************/
