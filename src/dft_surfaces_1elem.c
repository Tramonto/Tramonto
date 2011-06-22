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
 *  FILE: dft_surfaces_1elem.c
 *
 *  This file contains the routines used to set up point surfaces.
 *  These surfaces are approximated as one element surfaces in Ndim=1,2,3.
 *  so they are really approxmates to thin membranes, lines, and points.
 *
 */

#include "dft_surfaces_1elem.h"
/****************************************************************************/
void surface_1elemSurface_inSurfaceTest(int iwall,int iwall_type,
                                  double *fluid_testpos, double **wall_pos,
                                  double dist_adjustments, int flag_X_to_center,double *delr_vext,double *delr_zone,
                                  int *logical_inwall, int *logical_nearWallDielec)
{
  double x12,radius,r12sq_sum,r12;
  double roff=0.00000000001;
  int idim;
  struct SurfaceGeom_Struct *sgeom_iw;
  static int count_wallEl=0;


  sgeom_iw=&(SGeom[iwall_type]);

  if (Ipot_ff_c==COULOMB) *logical_nearWallDielec=FALSE;
  *logical_inwall=FALSE;

  radius=sgeom_iw->radius-roff;

  r12sq_sum=0.0;
  for (idim=0; idim<Ndim; idim++){
     x12 = wall_pos[iwall][idim] - fluid_testpos[idim];
     r12sq_sum += x12*x12;
  }
  r12=sqrt(r12sq_sum);

  if (r12<=radius+roff && count_wallEl==0){ 
     *logical_inwall=TRUE;  
     count_wallEl++;      
  }

  if (Ipot_ff_c==COULOMB && r12 < radius+Dielec_X) *logical_nearWallDielec  = TRUE;
  *delr_zone=r12-radius;
  if (flag_X_to_center==TRUE) *delr_vext=r12;
  else                        *delr_vext=r12-radius;

  return;
}
/****************************************************************************/
