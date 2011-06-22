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
 *  FILE: dft_surfaces_cylinder3D.c
 *
 *  This file contains the routines used to set up cylindrical surfaces in
 *  3D computational problems. The cylinders may be infinite or finite in length.
 *  A sine function may be added to the cylindrical surface for morphology studies.
 *
 */
#include "dft_surfaces_cylinder3D.h"
/****************************************************************************/
void surface_cylinder3D_inSurfaceTest(int iwall,int iwall_type,
                                  double *fluid_testpos, double **wall_pos,
                                  double dist_adjustments, int flag_X_to_center, double *delr_vext, double *delr_zone,
                                  int *logical_inwall, int *logical_nearWallDielec)
{

  double dx,x12,r12_sq_sum,r12,radius;
  double roff=0.00000000001,halflength,delr_sq,delr_sq_center;
  double xtest[3];
  int  dim[3];
  int npos,idim;
  int    orientation;
  struct SurfaceGeom_Struct *sgeom_iw;

  *logical_nearWallDielec=FALSE;
  *logical_inwall=FALSE;

  sgeom_iw=&(SGeom[iwall_type]);
  if (Lapply_offset[0]==TRUE) radius=sgeom_iw->radius+dist_adjustments-roff;   
  else                       radius=sgeom_iw->radius-roff;   
  if (Lapply_offset[1]==TRUE) halflength=sgeom_iw->halflength+dist_adjustments-roff;
  else                        halflength=sgeom_iw->halflength-roff;
  orientation=sgeom_iw->orientation;

  switch(orientation)  {
       case 2:
          xtest[0] = fluid_testpos[0]; dim[0] = 0;
          xtest[1] = fluid_testpos[1]; dim[1] = 1;
          xtest[2] = fluid_testpos[2]; dim[2] = 2;
          break;
       case 1:
          xtest[0] = fluid_testpos[2]; dim[0] = 2;
          xtest[1] = fluid_testpos[0]; dim[1] = 0;
          xtest[2] = fluid_testpos[1]; dim[2] = 1;
          break;
       case 0:
          xtest[0] = fluid_testpos[1]; dim[0] = 1;
          xtest[1] = fluid_testpos[2]; dim[1] = 2;
          xtest[2] = fluid_testpos[0]; dim[2] = 0;
          break;
  }

  r12_sq_sum = 0.0;
  for (idim = 0; idim < Ndim-1; idim++) {
     dx =  xtest[idim] -  wall_pos[iwall][dim[idim]];
     r12_sq_sum = r12_sq_sum + dx*dx;
  }
  r12 = sqrt(r12_sq_sum);
  x12 = fabs(wall_pos[iwall][dim[2]] - xtest[2]);

  if (r12 <= radius && x12 <= halflength) *logical_inwall=TRUE;
  if (Ipot_ff_c==COULOMB && r12 < radius+Dielec_X && x12<halflength+Dielec_X) *logical_nearWallDielec  = TRUE;

  delr_sq=0.0;
  delr_sq_center=0.0;
     /*dist to surface */
  if (r12 > radius) delr_sq+=(r12-radius)*(r12-radius);
  if (x12 > halflength) delr_sq+=(x12-halflength)*(x12-halflength);

     /*dist to center */
  if (r12 > radius) delr_sq_center+=(r12)*(r12);
  if (x12 > halflength) delr_sq_center+=(x12)*(x12);

  if (flag_X_to_center==TRUE) *delr_vext=sqrt(delr_sq_center);
  else                        *delr_vext=sqrt(delr_sq);
  *delr_zone=sqrt(delr_sq);

  return;
}

