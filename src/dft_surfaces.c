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
 *  FILE: dft_surfaces.c
 *
 *  This file contains the routines used to set up the surface of interest.  All
 *  of the geometry dependent information is found here.  All other routines
 *  are oblivious to surface geometry.
 *  
 */

#include "dft_surfaces.h"
/**********************************************************************/
void setup_surface (FILE *fp2, int *nelems_f,  
                    int **nelems_w_per_w, int **elems_f, 
                    int ***elems_w_per_w, int *elem_zones, int ***el_type)
{
 /* Local variable declarations */
  int itype,iwall,ilist,iel_box,izone, **L_wall,idim,inode_box;
  int loc_inode,real_wall,*imagetype,flag;
  int iel,inode,angle_test;
  double xtest[3],dist_roughness,dist_periodic, dist_linear,node_pos[3];
  double *x_min;

  int image,image_x,image_xy,i,nwall_max,image_old,icount,orientation,surfaceTypeID;
  double pos[3],**image_pos;
  struct SurfaceGeom_Struct *sgeom_iw;
  double delx,dist_adjustments;
  int logical_inwall=FALSE, logical_nearWallDielec=FALSE;

  double (*fp_roughness)(double *,int,int);
  void (*fp_inSurfaceTest)(int, int,int,int,double *, double **, double, double *, int *, int *);
  double (*fp_periodic)(double *,int,int);
  double (*fp_linear)(double *,int,int);
  int (*fp_angleCutout)(int, int,double *);

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
     for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) 
        for (iwall=0; iwall<Nwall; iwall++)  X_wall[loc_inode][iwall]=-999.0;
  }
  for (idim=0;idim<Ndim;idim++) xtest[idim]=0.0;

  for (iel_box=0; iel_box<Nelements_box; iel_box++){
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
     sgeom_iw = &(SGeom[itype]);
     orientation = sgeom_iw->orientation;
     surfaceTypeID = sgeom_iw->surfaceTypeID;
     if(Ndim > 1){
        if(surfaceTypeID == smooth_planar_wall) {
            for (idim=0; idim<Ndim; idim++) {
               if (idim != orientation) Xtest_reflect_TF[Link[iwall]][idim] = FALSE;
            } 
        }
        if(surfaceTypeID == finite_cyl_3D) {
           Xtest_reflect_TF[Link[iwall]][orientation] = FALSE;
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

     if (Xtest_reflect_TF[Link[iwall]][0]) find_wall_images(0,&image,image_pos,pos);

     if (Ndim > 1){
        image_x = image-image_old;
        for (i=image_old; i<image_old+image_x; i++){
           for (idim=0; idim<Ndim; idim++)
              pos[idim] = image_pos[i][idim];
              if (Xtest_reflect_TF[Link[iwall]][1]) find_wall_images(1, &image,image_pos,pos);
        }
     }
     if (Ndim == 3){
        image_xy = image-image_old;
        for (i=image_old; i<image_old+image_xy; i++){
           for (idim=0; idim<Ndim; idim++)
              pos[idim] = image_pos[i][idim];
              if (Xtest_reflect_TF[Link[iwall]][2]) find_wall_images(2,&image,image_pos,pos);
        }

     }
     for (i=image_old; i<image; i++) imagetype[i]=WallType[iwall];
     image_old=image;
  }
  nwall_max = image;

  icount=0;
  for (iwall=0; iwall<nwall_max; iwall++){
     itype = imagetype[iwall];
     sgeom_iw=&(SGeom[itype]);
     surfaceTypeID = sgeom_iw->surfaceTypeID;
  
     real_wall = iwall/(nwall_max/Nwall);
     MPI_Barrier(MPI_COMM_WORLD);     /* this may be a bug, but without this barrier statement we have issues with the surface areas */

  /* NEW CODE STARTS */

/*  would like to move switch to function, but it doesn't return the correct function pointers for some reason */
/* put switch here for now */
/*     setup_surface_fnc_pointers(sgeom_iw->surfaceTypeID,
                   fp_roughness, fp_periodic, fp_angleCutout, fp_inSurfaceTest);*/
  fp_roughness=NULL;
  fp_angleCutout=NULL;
  fp_inSurfaceTest=NULL;
  fp_periodic==NULL;
  switch(surfaceTypeID)
  {
     case smooth_planar_wall:
       fp_inSurfaceTest=&surface_planar_inSurfaceTest;
       fp_roughness=&surface_planar_roughness;
       if (Ndim==2){
          fp_angleCutout=&surface_angleCutout2D;
          fp_periodic=&surface_periodic_offset;
          fp_linear=&surface_linear_offset;
       }
       if (Ndim==3){
          fp_angleCutout=&surface_angleCutout3D_cyl;
          fp_periodic=&surface_periodic_offset;
          fp_linear=&surface_linear_offset;
       }
       break;

     case finite_planar_wall:
       fp_inSurfaceTest=&surface_block_inSurfaceTest;
       if (Ndim==2){
          fp_angleCutout=&surface_angleCutout2D;
          fp_periodic=&surface_periodic_offset;
          fp_linear=&surface_linear_offset;
       }
       if (Ndim==3){
          fp_angleCutout=&surface_angleCutout3D_cyl;
          fp_periodic=&surface_periodic_offset;
          fp_linear=&surface_linear_offset;
       }
       break;

     case colloids_cyl_sphere:
       if (Ndim==3) fp_inSurfaceTest=&surface_sphere_inSurfaceTest;
       else if (Ndim==2){
          fp_inSurfaceTest=&surface_cylinder2D_inSurfaceTest;
          fp_roughness=&surface_cylinder2D_roughness;
          fp_angleCutout=&surface_angleCutout2D;
          fp_linear=&surface_linear_offset;
       }
       break;

     case point_surface:
       fp_inSurfaceTest=&surface_1elemSurface_inSurfaceTest;
       break; 

     case finite_cyl_3D:
       fp_inSurfaceTest=&surface_cylinder3D_inSurfaceTest;
       fp_roughness=&surface_cylinder3D_roughness;
       fp_angleCutout=&surface_angleCutout3D_cyl;
       fp_periodic=&surface_periodic_offset;
       fp_linear=&surface_linear_offset;
       break;

     case atomic_centers:
       fp_inSurfaceTest=&surface_atoms_inSurfaceTest;
       break;

     case cyl2D_sphere3D_pore:
       if (Ndim==3) fp_inSurfaceTest=&surface_sphericalCavity3D_inSurfaceTest;
       else if (Ndim==2){
         fp_inSurfaceTest=&surface_cylindricalPore2D_inSurfaceTest;
         fp_roughness=&surface_cylinder2D_roughness;
         fp_angleCutout=&surface_angleCutout2D;
       }
       break;

     case cyl3D_slit2D_pore:
       if (Ndim==3){
         fp_inSurfaceTest=&surface_cylindricalPore3D_inSurfaceTest;
         fp_angleCutout=&surface_angleCutout3D_cyl;
         fp_roughness=&surface_cylinder3D_roughness;
         fp_linear=&surface_linear_offset;
       }
       else if (Ndim==2){
         fp_inSurfaceTest=&surface_slitPore2D_inSurfaceTest;
         fp_roughness=&surface_cylinder2D_roughness;
         fp_angleCutout=&surface_angleCutout2D;
         fp_linear=&surface_linear_offset;
       }
       break;

     default:           
       printf ("ERROR:the surface type chosen is not available\n");
       printf ("Surface type:%d\n",surfaceTypeID);
       exit(-1);
       break;
  }

     for (ilist=0; ilist<Nlists_HW; ilist++){
        for (iel_box = 0; iel_box < Nelements_box; iel_box++){

           iel=el_box_to_el(iel_box);
           inode = element_to_node(iel);
           node_to_position(inode,node_pos);
           inode_box = element_box_to_node_box(iel_box);

           loc_inode=B2L_node[inode_box];

           for (idim=0;idim<Ndim;idim++) xtest[idim]=node_pos[idim]+0.5*Esize_x[idim];

           /* compute adjustments to basic geometries*/
           dist_roughness=0.0; dist_periodic=0.0; dist_linear=0.0; angle_test=TRUE;
           if (sgeom_iw->Lrough_surface && fp_roughness!=NULL &&Ndim>1) dist_roughness=(*fp_roughness)(xtest,itype,iwall); 
           if (sgeom_iw->Lperiodic_overlay && fp_periodic!=NULL) dist_periodic=(*fp_periodic)(xtest,itype,iwall); 
           if (sgeom_iw->Llinear_overlay && fp_linear!=NULL) dist_linear=(*fp_linear)(xtest,itype,iwall); 
           if (sgeom_iw->Lwedge_cutout && fp_angleCutout!=NULL)  angle_test=(*fp_angleCutout)(real_wall,itype,xtest); 

           if (Lhard_surf && ilist !=Nlists_HW-1) dist_adjustments=0.5*Sigma_ff[ilist][ilist]; 
           else dist_adjustments=0.0;
           dist_adjustments+=(dist_roughness+dist_periodic+dist_linear);

           (*fp_inSurfaceTest)(iwall,itype,loc_inode,flag,xtest,image_pos,dist_adjustments,
                                                      &delx,&logical_inwall,&logical_nearWallDielec);

           if (logical_inwall==TRUE && angle_test==TRUE )  
                flag_wall_el(inode,ilist,real_wall,iel_box,L_wall,nelems_w_per_w, elems_w_per_w,el_type);

           else if (logical_nearWallDielec) Dielec[iel_box] = Dielec_pore;

           if (ilist==Nlists_HW-1){
              if (delx > 0.0) x_min[iel_box] = AZ_MIN(delx,x_min[iel_box]);
              else            x_min[iel_box] = 0.0;
           }

        }
     }
   
     if (icount==100 && Iwrite==VERBOSE) {
        printf("Proc: %d in setup surfaces:: iwall=%d of nwall_max=%d Nnodes_wall_box: %d \n",
                                             Proc,iwall,nwall_max,Nnodes_wall_box);
        icount=0;
     }
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
   int ijk[3],idim,nadd,i,new_wall,new_entry;
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
