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
 *  FILE: dft_block_planar.c
 *
 *  This file contains the routines used to set up planar surfaces of infinite
 *  extent in 2 dimensions
 *
 */

#include "dft_surfaces_block.h"
/****************************************************************************/
void surface_block_inSurfaceTest(int iwall,int iwall_type,
                                  double *fluid_testpos, double **wall_pos,
                                  double dist_adjustments, int flag_X_to_center, double *delx_vext, double *delx_zone,
                                  int *logical_inwall, int *logical_nearWallDielec)
{
  double x12,rsqsum_to_surface,rsqsum_to_center,xtest;
  double roff=0.00000000001,halfwidth;
  int npos,idim;
  struct SurfaceGeom_Struct *sgeom_iw;

  sgeom_iw=&(SGeom[iwall_type]);

  if (Ipot_ff_c==COULOMB) *logical_nearWallDielec=TRUE;
  *logical_inwall=TRUE;


  npos = 0;
  rsqsum_to_surface=0.0;
  rsqsum_to_center=0.0;
  for (idim=0; idim<Ndim; idim++){
     if (Lapply_offset[idim]==TRUE) halfwidth=sgeom_iw->halfwidth[idim]+dist_adjustments-roff;
     else                           halfwidth=sgeom_iw->halfwidth[idim]-roff;
     x12 = fabs(wall_pos[iwall][idim] - fluid_testpos[idim]);
     if (x12 > halfwidth)                                *logical_inwall = FALSE;
     if (Ipot_ff_c==COULOMB && x12 > halfwidth+Dielec_X) *logical_nearWallDielec  = FALSE;

     xtest=x12-halfwidth;
     if (xtest > 0) rsqsum_to_surface+=xtest*xtest;
     if (xtest > 0) rsqsum_to_center+=x12*x12;
  }
  if (*logical_inwall==TRUE) {
          *delx_vext = 0.0;   /* element in wall */
          *delx_zone = 0.0;   /* element in wall */
  }
  else{      
        if (flag_X_to_center==TRUE) *delx_vext = sqrt(rsqsum_to_center);  /*distance to center*/
        else                        *delx_vext = sqrt(rsqsum_to_surface);  /*distance to surface*/
        *delx_zone = sqrt(rsqsum_to_surface);  /*distance to surface*/
  }

  return;
}
/****************************************************************************/
