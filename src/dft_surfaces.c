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
  int itype,iwall,ilist,iel_box,izone, **L_wall,idim,inode_box,inode_box_tmp,save_dim,inode_global;
  int loc_inode,real_wall,flag,flag_X_to_center;
  int iel,inode,angle_test,jdim,save_el_wall;
  int edge_element, ijk[3],imax[3],ix[3],Linwall,Lfind_images[3],ijk_global[3],ijk_box_tmp[3];
  double xtest[3],dist_roughness,dist_periodic, dist_linear,node_pos[3],xtest_tmp[3];
  double *x_min;

  int image,image_x,image_xy,i,image_old,icount,orientation,surfaceTypeID;
  double pos[3],**image_pos;
  struct SurfaceGeom_Struct *sgeom_iw;
  double delx_vext,delx_zone,dist_adjustments;
  int logical_inwall=FALSE, logical_nearWallDielec=FALSE;

  double (*fp_roughness)(double *,int,int);
  void (*fp_inSurfaceTest)(int, int,double *, double **, double, int, double *, double *, int *, int *);
  double (*fp_periodic)(double *,int,int);
  double (*fp_linear)(double *,int,int);
  int (*fp_angleCutout)(int, int,double *,double **);

  /*  Calculate the number of wall and fluid elements and nodes, and 
      fill arrays that index the wall and fluid elements and nodes.  This
      is where the wall geometry comes in. */ 

  L_wall = (int **) array_alloc (2, Nlists_HW,Nelements_box, sizeof(int));
  x_min  = (double *) array_alloc (1, Nelements_box, sizeof(double));
  Touch_domain_boundary = (int ****) array_alloc (4,Nwall,Nlists_HW,Ndim,2,sizeof(int));

  flag=FALSE;
  for (i=0;i<Nwall_type;i++) if (Ipot_wf_n[i]==VEXT_DIST_TO_SURF || Ipot_wf_n[i]==VEXT_DIST_TO_CENTER)  flag=TRUE;
  for (idim=0;idim<Ndim;idim++) xtest[idim]=0.0;

  for (iel_box=0; iel_box<Nelements_box; iel_box++){
     if (Coarser_jac != JAC_ZONES_SETFIXED_ESIZE) elem_zones[iel_box] = Nzone - 1;
     else                  elem_zones[iel_box] = Nzone - 2;
     x_min[iel_box] = 1000.;
     for (ilist=0; ilist<Nlists_HW; ilist++){
        L_wall[ilist][iel_box] = FALSE;
        Wall_elems[ilist][iel_box] = -1;
        for (iwall=0; iwall<Nwall; iwall++) el_type[iwall][ilist][iel_box] = FLUID_EL;
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

  WallPos_Images = (double **) array_alloc (2, Nwall*POW_INT(3,Ndim),Ndim, sizeof(double));
  WallType_Images = (int *) array_alloc (1, Nwall*POW_INT(3,Ndim), sizeof(int));
  RealWall_Images = (int *) array_alloc (1, Nwall*POW_INT(3,Ndim), sizeof(int));
  Image_IDCheck = (int *) array_alloc (1, Nwall*POW_INT(3,Ndim), sizeof(int));
  /*Image_Link = (int *) array_alloc (1, Nwall*POW_INT(3,Ndim), sizeof(int));*/

  image=0;
  image_old = 0;
  for (iwall=0; iwall<Nwall; iwall++){  /* put real walls into image list first */
      for (idim=0;idim<Ndim;idim++) WallPos_Images[iwall][idim]=WallPos[idim][iwall];
      WallType_Images[iwall]=WallType[iwall];
      Image_IDCheck[iwall]=iwall;
      RealWall_Images[iwall]=iwall;
  }
  image=Nwall;

  for (iwall=0; iwall<Nwall; iwall++){

     for (idim=0;idim<Ndim;idim++) Lfind_images[idim]=TRUE;
     itype = WallType[iwall];
     sgeom_iw = &(SGeom[itype]);
     orientation = sgeom_iw->orientation;
     surfaceTypeID = sgeom_iw->surfaceTypeID;
     if (surfaceTypeID==smooth_planar_wall) {
         for (idim=0;idim<Ndim;idim++){
              if (idim!=orientation) Lfind_images[idim]=FALSE;
         }
     }
   
     image_old=image;
     for (i=0; i<image_old; i++){
         if (RealWall_Images[i]==iwall){
        for (idim=0; idim<Ndim; idim++)  pos[idim] = WallPos_Images[iwall][idim];
        if (Lfind_images[0]) find_wall_images(0,&image,WallPos_Images,pos);
        }
     }
     for (i=image_old;i<image;i++){
          WallType_Images[i]=WallType[iwall];
          RealWall_Images[i]=iwall;
          Image_IDCheck[i]=-1;
     }

     if (Ndim>1){
     image_old=image;
     for (i=0; i<image_old; i++){
         if (RealWall_Images[i]==iwall){
           for (idim=0; idim<Ndim; idim++) pos[idim] = WallPos_Images[i][idim];
           if (Lfind_images[1]) find_wall_images(1, &image,WallPos_Images,pos);
         }
     }
     for (i=image_old;i<image;i++){
          WallType_Images[i]=WallType[iwall];
          RealWall_Images[i]=iwall;
          Image_IDCheck[i]=-1;
     }
     }

     if (Ndim==3){ 
     image_old=image;
     for (i=0; i<image_old; i++){
         if (RealWall_Images[i]==iwall){
           for (idim=0; idim<Ndim; idim++) pos[idim] = WallPos_Images[i][idim];
           if (Lfind_images[2]) find_wall_images(2, &image,WallPos_Images,pos);
         }
     }
     for (i=image_old;i<image;i++){
          WallType_Images[i]=WallType[iwall];
          RealWall_Images[i]=iwall;
          Image_IDCheck[i]=-1;
     }
     }
  }
  Nwall_Images = image;

  for (iwall=0; iwall<Nwall_Images; iwall++)
    for (ilist=0; ilist<Nlists_HW; ilist++){
       nelems_w_per_w[ilist][iwall] = 0;
       if (Image_IDCheck[iwall]>=0){
         for (idim=0; idim<Ndim; idim++) {
           Touch_domain_boundary[Image_IDCheck[iwall]][ilist][idim][0]=FALSE;
           Touch_domain_boundary[Image_IDCheck[iwall]][ilist][idim][1]=FALSE;
         }
       }
  }

  if (Proc==0 && Iwrite==VERBOSE){
      printf("Number of surfaces requested=%d   Number of additional images=%d\n",Nwall,Nwall_Images-Nwall);
      for (iwall=0;iwall<Nwall_Images;iwall++){
         printf("image=%d  Image_IDCheck=%d  position=%g %g\n",iwall,Image_IDCheck[iwall],WallPos_Images[iwall][0],WallPos_Images[iwall][1]);
      }
  }

  if (flag){  /* if Vext is to be computed based on nearest x or r allocate array for walls and all images of walls */
     X_wall = (double **) array_alloc (2, Nnodes_box, Nwall_Images,sizeof(double));
     Xwall_delUP = (double ***) array_alloc (3, Nnodes_box, Nwall_Images,Ndim,sizeof(double));
     Xwall_delDOWN = (double ***) array_alloc (3, Nnodes_box, Nwall_Images,Ndim,sizeof(double));
     for (inode_box=0; inode_box<Nnodes_box; inode_box++) {
        for (iwall=0; iwall<Nwall_Images; iwall++){
                X_wall[inode_box][iwall]=1000.0;
                for (idim=0;idim<Ndim;idim++) {
                  Xwall_delUP[inode_box][iwall][idim]=1000.0;
                  Xwall_delDOWN[inode_box][iwall][idim]=1000.0;
                }
      }   }
  }

  icount=0;
  for (iwall=0; iwall<Nwall_Images; iwall++){
     itype = WallType_Images[iwall];
     sgeom_iw=&(SGeom[itype]);
     surfaceTypeID = sgeom_iw->surfaceTypeID;
  
     real_wall = RealWall_Images[iwall];
     if(Ipot_wf_n[itype]==VEXT_DIST_TO_CENTER) flag_X_to_center=TRUE;
     else                              flag_X_to_center=FALSE;
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
       if (Ndim>1){
          fp_angleCutout=&surface_angleCutout2D;
          fp_periodic=&surface_periodic_offset;
          fp_linear=&surface_linear_offset;
          fp_roughness=&surface_planar_roughness;
       }
       break;

     case finite_planar_wall:
       fp_inSurfaceTest=&surface_block_inSurfaceTest;
       if (Ndim>=1){
          fp_roughness=&surface_planar_roughness;
          fp_angleCutout=&surface_angleCutout2D;
          fp_periodic=&surface_periodic_offset;
          fp_linear=&surface_linear_offset;
       }
       break;

     case colloids_cyl_sphere:
       if (Ndim==3){
          fp_angleCutout=&surface_angleCutout3D_cyl;
          fp_inSurfaceTest=&surface_sphere_inSurfaceTest;
       }
       else if (Ndim==2){
          fp_inSurfaceTest=&surface_cylinder2D_inSurfaceTest;
          fp_roughness=&surface_cylinder2D_roughness;
          fp_angleCutout=&surface_angleCutout2D;
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
       if (Ndim==3){
         fp_angleCutout=&surface_angleCutout3D_cyl;
         fp_inSurfaceTest=&surface_sphericalCavity3D_inSurfaceTest;
       }
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

    
              /* find center of element - and set up all surface logical, wall elements, and x_min based on element center */
           node_to_position(inode,node_pos);
           inode_box = element_box_to_node_box(iel_box);
           loc_inode=B2L_node[inode_box];
           for (idim=0;idim<Ndim;idim++) xtest[idim]=node_pos[idim]+0.5*Esize_x[idim];

           /* compute adjustments to basic geometries*/
           dist_roughness=0.0; dist_periodic=0.0; dist_linear=0.0; angle_test=TRUE;
           if (sgeom_iw->Lrough_surface && fp_roughness!=NULL &&Ndim>1) dist_roughness=(*fp_roughness)(xtest,itype,iwall); 
           if (sgeom_iw->Lperiodic_overlay && fp_periodic!=NULL) dist_periodic=(*fp_periodic)(xtest,itype,iwall); 
           if (sgeom_iw->Llinear_overlay && fp_linear!=NULL) dist_linear=(*fp_linear)(xtest,itype,iwall); 
           if (sgeom_iw->Lwedge_cutout && fp_angleCutout!=NULL)  angle_test=(*fp_angleCutout)(iwall,itype,xtest,WallPos_Images); 

           /*if (Lhard_surf && ilist !=Nlists_HW-1) dist_adjustments=0.5*Sigma_wf[ilist][itype]; */
           if (Lhard_surf && ilist !=Nlists_HW-1) dist_adjustments=0.5*Sigma_ff[ilist][ilist]; 
           else dist_adjustments=0.0;
           dist_adjustments+=(dist_roughness+dist_periodic+dist_linear);

           (*fp_inSurfaceTest)(iwall,itype,xtest,WallPos_Images,dist_adjustments,flag_X_to_center,
                                                      &delx_vext,&delx_zone,&logical_inwall,&logical_nearWallDielec);

           if (logical_inwall==TRUE && angle_test==TRUE )  {
               if (surfaceTypeID != point_surface || nelems_w_per_w[ilist][iwall]==0){
                flag_wall_el(inode,ilist,iwall,iel_box,L_wall,nelems_w_per_w, elems_w_per_w,el_type);
               }
           }

           else if (logical_nearWallDielec) Dielec[iel_box] = Dielec_pore;

           if (ilist==Nlists_HW-1){
              if (delx_zone > 0.0) x_min[iel_box] = AZ_MIN(delx_zone,x_min[iel_box]);
              else            x_min[iel_box] = 0.0;
           }

               /* now set up X_wall array find minimum surface distances for nodes in the fluid */
               /*  first identify if this is an edge element where all nodes must be considered */
           if (flag==TRUE) {  
              node_to_ijk(inode,ijk);
              edge_element=FALSE;
              imax[0]=imax[1]=imax[2]=1;     
              for (idim=0; idim<Ndim; idim++){
                   if (ijk[idim]==Nodes_x[idim]-2 && Type_bc[idim][1] !=PERIODIC) imax[idim]=2;
              }

              for (ix[2]=0; ix[2]<imax[2]; ix[2]++){
                 for (ix[1]=0; ix[1]<imax[1]; ix[1]++){
                     for (ix[0]=0; ix[0]<imax[0]; ix[0]++){

                        for (idim=0;idim<Ndim;idim++){
                              xtest[idim]=node_pos[idim]+ix[idim]*Esize_x[idim];
                        }
                        (*fp_inSurfaceTest)(iwall,itype,xtest,WallPos_Images,dist_adjustments,flag_X_to_center,
                                                      &delx_vext,&delx_zone,&logical_inwall,&logical_nearWallDielec);
                         if (ilist==Nlists_HW-1){
                           inode_global=position_to_node(xtest);
                           node_to_ijk(inode_global,ijk_global);
                           ijk_to_ijk_box(ijk_global,ijk_box_tmp);
                           inode_box_tmp=ijk_box_to_node_box(ijk_box_tmp);

                           if (delx_vext > 0.0) X_wall[inode_box_tmp][iwall] = AZ_MIN(delx_vext,X_wall[inode_box_tmp][iwall]);
                           else            X_wall[inode_box_tmp][iwall] = 0.0;
                         }

                         /* now set up X_wall_delUP and X_wall_delDOWN to facilitate more accurate estimates of Vdash */
                         for (idim=0;idim<Ndim;idim++){
                            for (jdim=0;jdim<Ndim;jdim++){
                                 if (jdim==idim) xtest_tmp[jdim]=xtest[jdim]+VDASH_DELTA;
                                 else            xtest_tmp[jdim]=xtest[jdim];
                            }
                            (*fp_inSurfaceTest)(iwall,itype,xtest_tmp,WallPos_Images,dist_adjustments,flag_X_to_center,
                                                      &delx_vext,&delx_zone,&logical_inwall,&logical_nearWallDielec);
                             if (ilist==Nlists_HW-1){
                               if (delx_vext > 0.0) Xwall_delUP[inode_box_tmp][iwall][idim] = AZ_MIN(delx_vext,Xwall_delUP[inode_box_tmp][iwall][idim]);
                               else            Xwall_delUP[inode_box_tmp][iwall][idim] = 0.0;
                             }
                         }
                         for (idim=0;idim<Ndim;idim++){
                            for (jdim=0;jdim<Ndim;jdim++){
                                 if (jdim==idim) xtest_tmp[jdim]=xtest[jdim]-VDASH_DELTA; 
                                 else            xtest_tmp[jdim]=xtest[jdim];
                            }
                            (*fp_inSurfaceTest)(iwall,itype,xtest_tmp,WallPos_Images,dist_adjustments,flag_X_to_center,
                                                      &delx_vext,&delx_zone,&logical_inwall,&logical_nearWallDielec);
                             if (ilist==Nlists_HW-1){
                               if (delx_vext > 0.0) Xwall_delDOWN[inode_box_tmp][iwall][idim] = AZ_MIN(delx_vext,Xwall_delDOWN[inode_box_tmp][iwall][idim]);
                               else            Xwall_delDOWN[inode_box_tmp][iwall][idim] = 0.0;
                             }
                         }


                      }
                  }
               }


             } /* logical_local_node */
        }  /* iel_box loop */
     }    /* ilist loop */
   
     if (icount==100 && Iwrite==VERBOSE) {
        printf("Proc: %d in setup surfaces:: iwall=%d of Nwall_Images=%d Nnodes_wall_box: %d \n",
                                             Proc,iwall,Nwall_Images,Nnodes_wall_box);
        icount=0;
     }
     else if (iwall==Nwall_Images-1 && Iwrite==VERBOSE && Proc==0){
        printf("Proc: %d is done with setup surfaces:: iwall=%d of Nwall_Images=%d Nnodes_wall_box=%d\n",
                                              Proc,iwall,Nwall_Images,Nnodes_wall_box);
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
  safe_free((void *) &L_wall);
  safe_free((void *) &x_min);
  
}  /* end of surface element setup*/
/****************************************************************************/
/*find_wall_images: here we take the surface positions from the
 *    central cell and locate all the corresponding images due
 *    to reflective and periodic boundaries.
 */
void find_wall_images(int idim, int *image, double **image_pos, double *pos)
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

/*   if (iwall < Nwall) {*/
      elems_w_per_w[ilist][nelems_w_per_w[ilist][iwall]][iwall] = iel_box;
      nelems_w_per_w[ilist][iwall]++;
      Wall_elems[ilist][iel_box] = iwall;
      L_wall[ilist][iel_box] = TRUE;
      el_type[iwall][ilist][iel_box] = WALL_EL;


      /* setup wall dielectric constants */
      if (ilist == Nlists_HW-1 && Ipot_ff_c == COULOMB){
          if (Type_dielec != DIELEC_CONST)
                  Dielec[iel_box] = Dielec_wall[WallType_Images[iwall]];
      }

      /* flag cases where a given wall touches the domain boundary - need
         to do some special things if we touch a reflective boundary */
      node_to_ijk(inode,ijk);

      if (RealWall_Images[iwall]<Nwall){
         for (idim=0; idim<Ndim; idim++){
         if (ijk[idim] == 0 && !Touch_domain_boundary[RealWall_Images[iwall]][ilist][idim][0]) 
                       Touch_domain_boundary[RealWall_Images[iwall]][ilist][idim][0] = TRUE;
             if (Type_bc[idim][1]==PERIODIC) nadd=0;
             else                            nadd=1;
             if (ijk[idim]+nadd == Nodes_x[idim]-1 && 
                       !Touch_domain_boundary[RealWall_Images[iwall]][ilist][idim][1]) 
                       Touch_domain_boundary[RealWall_Images[iwall]][ilist][idim][1] = TRUE;
         }
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
  /* }*/
   return;
}
/****************************************************************************/
