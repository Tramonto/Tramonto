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
 *  FILE: dft_surfacesAdjust_angleCutout.c
 *
 *  This file contains the routines used to make adjustments to cylindrical
 *  surfaces and pores.  These routines allow the definition of a partial
 *  cylinder rather than a full cylinder.  Matching two or more of these
 *  partial cylinders together allows for chemical variation around the 
 *  circumference of the cylinder.
 *
 */

#include "dft_surfacesAdjust_angleCutout.h"
/****************************************************************************/
int surface_angleCutout2D(int iwall,int iwall_type,double *fluidEl_center, double **wall_pos)
{
   int angle_test,idim;
   double vecB[2],angle,cos_angle,shift,angle_deg,rsq,r;
   double angle_wedge_start,angle_wedge_end,angle_reflect;
   struct SurfaceGeom_Struct *sgeom_iw;

   sgeom_iw=&(SGeom[iwall_type]);
   angle_wedge_start=sgeom_iw->angle_wedge_start;
   angle_wedge_end=sgeom_iw->angle_wedge_end;

   rsq=0.0;
   for (idim=0;idim<Ndim;idim++){
         vecB[idim]=fluidEl_center[idim]-wall_pos[idim][iwall];
         rsq+=vecB[idim]*vecB[idim];
   }
   r=sqrt(rsq);

   idim=1;
   cos_angle=vecB[0]/r;
   angle = acos(cos_angle);
   angle_deg=180.*angle/PI;
   if (fluidEl_center[1] <wall_pos[1][iwall]){
       angle_reflect=360-angle_deg;
       angle_deg=angle_reflect;
   }

   if(angle_wedge_end>angle_wedge_start){
          angle_test=TRUE;
          if(angle_deg <angle_wedge_start || angle_deg > angle_wedge_end) angle_test=FALSE;
   }
   else if(angle_wedge_end<angle_wedge_start){
          angle_test=FALSE;
          if(angle_deg >angle_wedge_start || angle_deg < angle_wedge_end) angle_test=TRUE;
   }
   return (angle_test);
}
/****************************************************************************/
int surface_angleCutout3D_cyl(int iwall,int iwall_type,double *fluidEl_center,double **wall_pos)
{
  int angle_test,idim,orientation,dim[3],idim_testCos;
  double vecB[3],xtest[3],angle,cos_angle,angle_deg,rsq,r;
  double angle_wedge_start,angle_wedge_end,angle_reflect;
  struct SurfaceGeom_Struct *sgeom_iw;


  sgeom_iw=&(SGeom[iwall_type]);
  angle_wedge_start=sgeom_iw->angle_wedge_start;
  angle_wedge_end=sgeom_iw->angle_wedge_end;
  orientation=sgeom_iw->orientation;

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

   rsq=0.0;
   for (idim=0;idim<Ndim-1;idim++){
         vecB[idim]=xtest[idim]-wall_pos[dim[idim]][iwall];
         rsq+=vecB[idim]*vecB[idim];
   }
   r=sqrt(rsq);

   cos_angle=vecB[0]/r;
   angle = acos(cos_angle);
   angle_deg=180.*angle/PI;

   if (fluidEl_center[idim_testCos] <wall_pos[idim_testCos][iwall]){
       angle_reflect=360-angle_deg;
       angle_deg=angle_reflect;
   }

   if(angle_wedge_end>angle_wedge_start){
          angle_test=TRUE;
          if(angle_deg <angle_wedge_start || angle_deg > angle_wedge_end) angle_test=FALSE;
   }
   else if(angle_wedge_end<angle_wedge_start){
          angle_test=FALSE;
          if(angle_deg >angle_wedge_start || angle_deg < angle_wedge_end) angle_test=TRUE;
   }
   return (angle_test);
}
/****************************************************************************/

