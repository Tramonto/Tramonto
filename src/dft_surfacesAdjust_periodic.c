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
 *  FILE: dft_surfacesAdjust_periodic.c
 *
 *  This file contains the routines used to superimpose a periodic function
 *  on a basic surface shape.
 *
 */
#include "dft_surfacesAdjust_periodic.h"
/****************************************************************************/
double surface_periodic_offset(double *fluidEl_center,int iwall_type,int iwall)
{
  struct SurfaceGeom_Struct *sgeom_iw;
  double amplitude,wavelength,offset,origin;
  int    orientation,nperiodic_fnc,i;

  sgeom_iw = &(SGeom[iwall_type]);
  nperiodic_fnc=sgeom_iw->Nperiodic_overlay;
  offset=0.0;

  for(i=0;i<nperiodic_fnc;i++){
     orientation=sgeom_iw->orientation_periodic[i];
     amplitude=sgeom_iw->amplitude[i];
     wavelength=sgeom_iw->wavelength[i];
     origin=sgeom_iw->origin_PeriodicFunc[i];

     offset+=amplitude*cos(2*PI*(fluidEl_center[orientation]-origin)/wavelength);
  }
  return(offset); 
}
/****************************************************************************/

