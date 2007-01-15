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
 *  FILE: dft_images.c
 *
 *  This file contains routines to find images of a given surface node or
 *  so that external fields will correctly include all images.
 */

#include "dft_images.h"

/******************************************************************/
/*find_images: here we take all the element positions in image_pos
 *          and then find the images of these elements in the
 *          direction idim.  The counter image tells us how many
 *          entries there are in image_pos.
 *
 *        Note: the variable node_ref refers to the fixed node
 *              in the domain against which images (found in 
 *              node_image) will be compared to determine if 
 *              the image should be included !!
 */
void find_images(int idim, double cut, 
                   int *image, double **image_pos,
                   double *node_image, double *node_ref) 
{
  int iside,done,i,jdim;
  double sign,shift,node_image2[3],rsq,r,r_wf,r_wf_sq=0.0;

  for (jdim=0; jdim<Ndim; jdim++) 
     r_wf_sq += (node_image[jdim]-node_ref[jdim])
               *(node_image[jdim]-node_ref[jdim]);
  r_wf = sqrt(r_wf_sq);

  for (iside = 0; iside < 2; iside++){
     if (iside == 0) sign=-1.0; else sign=1.0;

     if (Type_bc[idim][iside] == PERIODIC ||
         (Type_bc[idim][iside] == REFLECT) || 
         (Type_bc[idim][iside] == IN_WALL  && 
              (fabs(node_image[idim]-sign*(0.5*Size_x[idim]-0.5*Esize_x[idim]))< 0.0001))){

        done = FALSE;
        i=1;

        if (Type_bc[idim][iside] == PERIODIC)
           shift = Size_x[idim];
        else if (Type_bc[idim][iside] == REFLECT)
           shift = 2.0*fabs(sign*0.5*Size_x[idim] - node_image[idim]);
        else 
           shift = Esize_x[idim];
                     

        while (!done){
           for (jdim=0; jdim<Ndim; jdim++) node_image2[jdim] = node_image[jdim];
           node_image2[idim] = node_image[idim] + i*sign*shift;

           rsq=0.0;
           for (jdim=0; jdim<Ndim; jdim++) 
              rsq += (node_ref[jdim]-node_image2[jdim])*
                     (node_ref[jdim]-node_image2[jdim]);
           r = sqrt(rsq);

           if (r<=cut){
              for (jdim=0; jdim<Ndim; jdim++)
                 image_pos[*image][jdim] = node_image2[jdim];
              (*image)++;
           }
           else done = TRUE;

           if (Type_bc[idim][iside] == REFLECT) done = TRUE;
           i++;    
        }

     }  /* end of check for b.c. type */
   }  /* end of left side / right side test */

   return;
}
/******************************************************************/
/*find_images_1D: a special case of find_images where we have
 *                a 2D or 3D domain, but only a 1D external field.
 */
void find_images_1D(int idim, double cut, 
                   int *image, double **image_pos,
                   double *node_image, double *node_ref) 
{
  int iside,done,i,jdim;
  double sign,shift,node_image2[3],r;

  for (iside = 0; iside < 2; iside++){
     if (iside == 0) sign=-1.0; else sign=1.0;

     if (Type_bc[idim][iside] == PERIODIC || Type_bc[idim][iside] == REFLECT  ){

        done = FALSE;
        i=1;

        if (Type_bc[idim][iside] == PERIODIC)  shift = Size_x[idim];
        else if (Type_bc[idim][iside] == REFLECT)
           shift = 2.0*fabs(sign*0.5*Size_x[idim] - node_image[idim]);

        while (!done){
           for (jdim=0; jdim<Ndim; jdim++) node_image2[jdim] = node_image[jdim];
           node_image2[idim] = node_image[idim] + i*sign*shift;

           r = fabs((node_ref[idim]-node_image2[idim]));

           if (r<=cut){
              for (jdim=0; jdim<Ndim; jdim++)
                 image_pos[*image][jdim] = node_image2[jdim];
              (*image)++;
           }
           else done = TRUE;

           if (Type_bc[idim][iside] == REFLECT) done = TRUE;
           i++;    
        }

     }  /* end of check for b.c. type */
   }  /* end of left side / right side test */

   return;
}
/******************************************************************/
/*find_images_coulomb: 
 *   find images of a given charged element across reflective 
 *   boundaries only.  There are no cutoffs because Coulomb
 *   potentials are of infinite extent.  We only get here with
 *   reflective/bulk or all bulk boundary conditions on the domain.
 */
void find_images_coulomb(int idim, int *image, double **image_pos,double *node_image) 
{
  int iside,done,i,jdim;
  double sign,shift,node_image2[3],rsq,r;

  for (iside = 0; iside < 2; iside++){
     if (iside == 0) sign=-1.0; else sign=1.0;
                
     if (Type_bc[idim][iside] == REFLECT){
        done = FALSE;
        i=1;


        shift = 2.0*fabs(sign*0.5*Size_x[idim] - node_image[idim]);
                     
        while (!done){
           for (jdim=0; jdim<Ndim; jdim++) node_image2[jdim] = node_image[jdim];
           node_image2[idim] = node_image[idim] + i*sign*shift;

           for (jdim=0; jdim<Ndim; jdim++){
               image_pos[*image][jdim] = node_image2[jdim];
           }
           (*image)++;

           if (Type_bc[idim][iside] == REFLECT) done = TRUE;
           i++;    
        }

     }  /* end of check for b.c. type */
   }  /* end of left side / right side test */

   return;
}
/***********************************************************************/
 /*find_images2: here we take all the element positions in image_pos
 *          and then find the images of these elements in the
 *          direction idim.  The counter image tells us how many
 *          entries there are in image_pos.
 *
 *        In this routine we consider only IN_WALL images ... this
 *        is needed for calculations of wall-wall interactions.
 *
 *        In this case, the reference position is set by the position
 *        of iwall.
 */
void find_images2(int idim, double cut, 
                  int *image, double **image_pos,
                  double *node_image, 
                  int iwall,int iside) 
{
  int done,i,jdim;
  double sign,shift,node_image2[3],rsq,r,r_ww,r_ww_sq=0.0;

  r_ww_sq += (node_image[idim]-WallPos[idim][iwall])
               *(node_image[0]-WallPos[idim][iwall]);
  r_ww = sqrt(r_ww_sq);

  if (iside == 0) sign=-1.0; else sign=1.0;
                
  if ( Type_bc[idim][iside] == IN_WALL  && 
       fabs(node_image[idim]+sign*0.5*Esize_x[idim] 
                            -sign*0.5*Size_x[idim]) < 0.0001  ){

    done = FALSE;
    i=1;

    shift = Esize_x[idim];
                     
    while (!done){
        for (jdim=0; jdim<Ndim; jdim++) node_image2[jdim] = node_image[jdim];
        node_image2[idim] = node_image[idim] + i*sign*shift;

        rsq=0.0;
        rsq += (WallPos[idim][iwall]-node_image2[idim])*
               (WallPos[idim][iwall]-node_image2[idim]);
        r = sqrt(rsq);

        if ( r-WallParam[WallType[iwall]] <= cut) {
              for (jdim=0; jdim<Ndim; jdim++)
                 image_pos[*image][jdim] = node_image2[jdim];
              (*image)++;
           }
           else done = TRUE;

           i++;    
        }

     }  /* end of check for b.c. type */

   return;
}
/***********************************************************************/
