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
 *  FILE: dft_surfaces_cylindricalPore2D.c
 *
 *  This file contains the routines used to set up cylindrical pores of infinite
 *  extent (2 dimensional computational problem).
 *
 */

#include "dft_surfaces_cylindricalPore2D.h"
/****************************************************************************/
void surface_cylindricalPore2D_inSurfaceTest(int iwall,int iwall_type,
                                  double *fluid_testpos, double **wall_pos,
                                  double dist_adjustments, int flag_X_to_center, double *delr_vext, double *delr_zone,
                                  int *logical_inwall, int *logical_nearWallDielec)
{
  double x12,r12sq,xtest,r12,radius;
  double roff=0.00000000001,halfwidth;
  int npos,idim;
  struct SurfaceGeom_Struct *sgeom_iw;

  sgeom_iw=&(SGeom[iwall_type]);

  *logical_nearWallDielec=FALSE;
  *logical_inwall=FALSE;

   radius=sgeom_iw->radius+dist_adjustments-roff;

  r12sq=0.0;
  for (idim=0; idim<Ndim; idim++){
     x12 = wall_pos[iwall][idim] - fluid_testpos[idim];
     r12sq += x12*x12;
  }
  r12=sqrt(r12sq);

  if (r12>=radius) *logical_inwall=TRUE;  
  if (Ipot_ff_c==COULOMB && r12 > radius-Dielec_X) *logical_nearWallDielec  = TRUE;
  *delr_vext=radius-r12;
  *delr_zone=radius-r12;

  return;
}
/****************************************************************************/
