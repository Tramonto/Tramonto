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
 *  FILE: dft_surfaces_sphere.c
 *
 *  This file contains the routines used to set up a spherical surface
 *  in a 3-dimensional problem. Use this to study colloidal fluids.
 *
 */

#include "dft_surfaces_sphere.h"
/****************************************************************************/
void surface_sphere_inSurfaceTest(int iwall,int iwall_type,int loc_inode,int flag_setup_Xwall,
                                  double *fluidEl_center, double **image_pos,
                                  double dist_adjustments, double *delr,
                                  int *logical_inwall, int *logical_nearWallDielec)
{
  double x12,radius,r12sq_sum,r12;
  double roff=0.00000000001,halfwidth;
  int idim;
  struct SurfaceGeom_Struct *sgeom_iw;

  sgeom_iw=&(SGeom[iwall_type]);

  if (Ipot_ff_c==COULOMB) *logical_nearWallDielec=FALSE;
  *logical_inwall=FALSE;

  radius=sgeom_iw->radius+dist_adjustments-roff;

  r12sq_sum=0.0;
  for (idim=0; idim<Ndim; idim++){
     x12 = image_pos[iwall][idim] - fluidEl_center[idim];
     r12sq_sum += x12*x12;
  }
  r12=sqrt(r12sq_sum);

  if (r12<=radius) *logical_inwall=TRUE;  
  if (Ipot_ff_c==COULOMB && r12 < radius+Dielec_X) *logical_nearWallDielec  = TRUE;
  *delr=r12-radius;

  if (flag_setup_Xwall){
      printf("error : the Vext_1D_xmin option is not available for spherical surfaces\n");
      printf("try the Vext_integrated surface instead.\n");
      exit(-1);
  }
  return;
}
/****************************************************************************/
