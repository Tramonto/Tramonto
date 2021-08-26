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
 *  FILE: dft_surfaces_planar.c
 *
 *  This file contains the routines used to set up planar surfaces of infinite
 *  extent in 2 dimensions
 *
 */

#include "dft_surfaces_planar.h"
/****************************************************************************/
void surface_planar_inSurfaceTest(int iwall,int iwall_type,
                                  double *fluid_testpos, double **wall_pos,
                                  double dist_adjustments, int flag_X_to_center, double *delx_vext, double *delx_zone,
                                  int *logical_inwall, int *logical_nearWallDielec)
{
  int orientation;
  double halfwidth,x12;
  double roff=0.00000000001;
  struct SurfaceGeom_Struct *sgeom_iw;

  *logical_inwall=FALSE;
  *logical_nearWallDielec=FALSE;

  sgeom_iw=&(SGeom[iwall_type]);

  orientation=sgeom_iw->orientation;
  halfwidth=sgeom_iw->halfwidth[orientation];

  halfwidth +=dist_adjustments-roff;


  x12 = fabs(wall_pos[iwall][orientation] - fluid_testpos[orientation]);

  if (x12 <= halfwidth) *logical_inwall=TRUE;
  if (Type_dielec==DIELEC_WF_PORE && x12<halfwidth + Dielec_X)  *logical_nearWallDielec=TRUE;

  if (flag_X_to_center==TRUE) *delx_vext=x12;
  else                        *delx_vext=x12-halfwidth;
  *delx_zone=x12-halfwidth;

  return;
}
/****************************************************************************/
