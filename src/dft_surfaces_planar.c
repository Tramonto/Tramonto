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
void surface_planar_inSurfaceTest(int iwall,int iwall_type,int loc_inode, int flag_setup_Xwall,
                                  double *fluidEl_center, double **image_pos,
                                  double dist_adjustments, double *delx,
                                  int *logical_inwall, int *logical_nearWallDielec)
{
  int orientation;
  double halfwidth,x12;
  double roff=0.00000000001;
  struct SurfaceGeom_Struct *sgeom_iw;
  int imax[3],ix[3],idim,inode_box,ijk[3],ijk_tmp[3],ijk_tmp_box[3];

  sgeom_iw=&(SGeom[iwall_type]);

  orientation=sgeom_iw->orientation;
  halfwidth=sgeom_iw->halfwidth[orientation];

  halfwidth +=dist_adjustments-roff;

  /* get distance from surface to center of element in appropriate dimension */
  x12 = fabs(image_pos[iwall][orientation] - fluidEl_center[orientation]);
  if (x12 <= halfwidth) *logical_inwall=TRUE;
  else *logical_inwall=FALSE;

  if (Type_dielec==DIELEC_WF_PORE && x12<halfwidth + Dielec_X)  *logical_nearWallDielec=TRUE;
  else *logical_nearWallDielec=FALSE;
  *delx=x12-halfwidth;

           /* need to set up X_wall here because it depends on specific surface geometries. */
  if (loc_inode >= 0 && flag_setup_Xwall) {
     imax[0]=imax[1]=imax[2]=1;
     node_to_ijk(L2G_node[loc_inode],ijk);
     for (idim=0; idim<Ndim; idim++) imax[idim]=2;

     for (ix[2]=0; ix[2]<imax[2]; ix[2]++){
      for (ix[1]=0; ix[1]<imax[1]; ix[1]++){
       for (ix[0]=0; ix[0]<imax[0]; ix[0]++){
           for (idim=0;idim<Ndim;idim++) ijk_tmp[idim]=ijk[idim]+ix[idim];                
               ijk_to_ijk_box(ijk_tmp,ijk_tmp_box);
               inode_box=ijk_box_to_node_box(ijk_tmp_box);
               if (B2L_node[inode_box]>=0){
  
               if (*delx>0){
                    if (fluidEl_center[orientation]>image_pos[iwall][orientation]){  
                       if (ix[orientation]==0) X_wall[B2L_node[inode_box]][iwall]=*delx-0.5*Esize_x[orientation];              
                       else{
                           X_wall[B2L_node[inode_box]][iwall]=*delx+0.5*Esize_x[orientation];
                       }
                    }
                    else{
                       if (ix[orientation]==0) X_wall[B2L_node[inode_box]][iwall]=*delx+0.5*Esize_x[orientation];              
                       else X_wall[B2L_node[inode_box]][iwall]=*delx-0.5*Esize_x[orientation];
                    }
               }
               else X_wall[B2L_node[inode_box]][iwall]=0.0;
               }
     }}}
  }

  return;
}
/****************************************************************************/
