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
 *  FILE: dft_surfacesAdjust_linear.c
 *
 *  This file contains the routines used to superimpose a linear function
 *  on a basic surface shape.
 *
 */
#include "dft_surfacesAdjust_linear.h"
/****************************************************************************/
double surface_linear_offset(double *fluidEl_center,int iwall_type,int iwall)
{
  struct SurfaceGeom_Struct *sgeom_iw;
  double slope,offset,origin,endpoint;
  int    orientation,nlinear_fnc,i;

  sgeom_iw = &(SGeom[iwall_type]);
  nlinear_fnc=sgeom_iw->Nlinear_overlay;
  offset=0.0;

  for(i=0;i<nlinear_fnc;i++){
     orientation=sgeom_iw->orientation_linear[i];
     slope=sgeom_iw->slope[i];
     origin=sgeom_iw->origin_LinearFunc[i];
     endpoint=sgeom_iw->endpoint_LinearFunc[i];

     if ( (origin<fluidEl_center[orientation] && fluidEl_center[orientation]<=endpoint) ||
          (origin>fluidEl_center[orientation] && fluidEl_center[orientation]>=endpoint)  )
           offset+=slope*(fluidEl_center[orientation]-origin);
  }
  return(offset); 
}
/****************************************************************************/

