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
 *  FILE: dft_surfacesAdjust_roughness.c
 *
 *  This file contains the routines used to introduce roughness for various surface
 *  types.
 *
 */

#include "dft_surfacesAdjust_roughness.h"
/****************************************************************************/
double surface_planar_roughness(double *fluidEl_center,int iwall_type,int iwall)
{
  int iblock,rough_block[2],idim,orientation;
  double roughness,rough_length;
  struct SurfaceGeom_Struct *sgeom_iw;

  sgeom_iw=&(SGeom[iwall_type]);
  rough_length=sgeom_iw->roughness_length;
  orientation=sgeom_iw->orientation;

  for (iblock=0;iblock<2;iblock++) rough_block[iblock]=0;
  iblock=0;
  for (idim=0;idim<Ndim;idim++){
     if (idim != orientation) {
         rough_block[iblock] = (int) ((fluidEl_center[idim]+0.5*Size_x[idim])/rough_length);
         if (rough_block[iblock] >= MAX_ROUGH_BLOCK) {
              printf("ERROR with rough surfacess - number of rough patches exceeds MAX_ROUGH_BLOCK=%d\n",
              MAX_ROUGH_BLOCK);
              exit(-1);
         }
         iblock++;
      }
   }
   roughness = Rough_precalc[iwall_type][rough_block[0]][rough_block[1]];
   return roughness;
}
/****************************************************************************/
double surface_cylinder2D_roughness(double *fluidEl_center,int iwall_type,int iwall)
{
  int iblock,rough_block[2],idim;
  double roughness,rough_length,rsq,radius,shift;
  double vecB[2],r,angle,cos_angle,angle_block,angle_deg,angle_reflect;
  struct SurfaceGeom_Struct *sgeom_iw;

  sgeom_iw=&(SGeom[iwall_type]);
  rough_length=sgeom_iw->roughness_length;
  radius=sgeom_iw->radius;


  for (iblock=0;iblock<2;iblock++) rough_block[iblock]=0;
  rough_block[1] = 0;

  if (rough_length > radius){
    printf("ERROR IN LENGTH SCALE FOR ROUGHNESS...Rough length larger than radius\n");
    exit(-1);
  }
  
  angle_block = asin(rough_length/radius);

  rsq=0.0;
  for (idim=0;idim<Ndim;idim++){
     vecB[idim]=fluidEl_center[idim]-WallPos[idim][iwall];
    rsq+=vecB[idim]*vecB[idim];
  }
  r=sqrt(rsq);

  cos_angle=vecB[0]/r;
  angle = acos(cos_angle);
  angle_deg=180.*angle/PI;
  if (fluidEl_center[1] <WallPos[1][iwall]){
      angle_reflect=360-angle_deg;
      angle_deg=angle_reflect;
  }
  angle=angle_deg*PI/180.;

  rough_block[0]= (int)(angle/angle_block);

  if (rough_block[0] >= MAX_ROUGH_BLOCK || rough_block[1]>=MAX_ROUGH_BLOCK) {
         printf("ERROR with rough cylinder - number of rough patches exceeds MAX_ROUGH_BLOCK %d\n", MAX_ROUGH_BLOCK);
         exit(-1);
  }
  roughness = Rough_precalc[iwall_type][rough_block[0]][rough_block[1]];
  return (roughness);
}
/****************************************************************************/
double surface_cylinder3D_roughness(double *fluidEl_center,int iwall_type,int iwall)
{
  double roff=0.00000000001;
  int iblock,rough_block[2],idim,orientation,dim[3],idim_testCos;
  double roughness,rough_length,rsq,radius,halflength;
  double vecB[3],r,angle,cos_angle,angle_block,xtest[3],angle_deg,angle_reflect;
  struct SurfaceGeom_Struct *sgeom_iw;

  sgeom_iw=&(SGeom[iwall_type]);
  rough_length=sgeom_iw->roughness_length;
  radius=sgeom_iw->radius;
  orientation=sgeom_iw->orientation;
  halflength=sgeom_iw->halflength-roff;

  switch(orientation)  {
       case 2:
          xtest[0] = fluidEl_center[0]; dim[0] = 0;
          xtest[1] = fluidEl_center[1]; dim[1] = 1;
          xtest[2] = fluidEl_center[2]; dim[2] = 2;
          idim_testCos=1;
          break;
       case 1:
          xtest[0] = fluidEl_center[2]; dim[0] = 2;
          xtest[1] = fluidEl_center[0]; dim[1] = 0;
          xtest[2] = fluidEl_center[1]; dim[2] = 1;
          idim_testCos=0;
          break;
       case 0:
          xtest[0] = fluidEl_center[1]; dim[0] = 1;
          xtest[1] = fluidEl_center[2]; dim[1] = 2;
          xtest[2] = fluidEl_center[0]; dim[2] = 0;
          idim_testCos=2;
          break;
  }

  for (iblock=0;iblock<2;iblock++) rough_block[iblock]=0;
  rough_block[1] = 0;

  if (rough_length > halflength/2.){
    printf("ERROR IN LENGTH SCALE FOR ROUGHNESS...Rough length larger than 1/4 of the length of cylinder\n");
    exit(-1);
  }
  if (rough_length > radius){
    printf("ERROR IN LENGTH SCALE FOR ROUGHNESS...Rough length larger than the cylinder radius\n");
    exit(-1);
  }
 
  angle_block = asin(rough_length/radius);

  rsq=0.0;
  for (idim=0;idim<Ndim-1;idim++){
    vecB[idim]=xtest[idim]-WallPos[dim[idim]][iwall];
    rsq+=vecB[idim]*vecB[idim];
  }
  r=sqrt(rsq); 

  cos_angle=vecB[0]/r;
  angle = acos(cos_angle);
  angle_deg=180.*angle/PI;

  if (fluidEl_center[idim_testCos] <WallPos[idim_testCos][iwall]){
      angle_reflect=360-angle_deg;
      angle_deg=angle_reflect;
  }
  angle=PI*angle_deg/180;

  rough_block[0]= (int)(angle/angle_block);
  rough_block[1] = (int) ((fluidEl_center[orientation]+0.5*Size_x[orientation])/rough_length);

  if (rough_block[0] >= MAX_ROUGH_BLOCK || rough_block[1]>=MAX_ROUGH_BLOCK) {
         printf("ERROR with rough cylinder - number of rough patches exceeds MAX_ROUGH_BLOCK %d\n", MAX_ROUGH_BLOCK);
         exit(-1);
  }
  roughness = Rough_precalc[iwall_type][rough_block[0]][rough_block[1]];
  return (roughness);
}
/****************************************************************************/

