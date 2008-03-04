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
 *  FILE: dft_mesh_lib.c
 *
 *  This file contains a collection of functions that may be useful for 
 *  negotiating a regular mesh, converting between nodes and elements, 
 *  positions, and locating node numbers from stencil offsets.  
 *  These routines were originally found in dft_mesh.c. Included are:
 *
 * round_to_int:      round a double precision number the nearest integer. 
 * offset_to_node:    given a stencil offset, find the node number.
 *                    (this calls a routine "boundary_condition" which applies
 *                    the correct -periodic or constant- boundary condition)
 * element_to_node:   given an element, find the lower left back node number
 * node_to_elem:      given a node, and a corner, find the associated element
 *                    corners 0,1 (left right ... 1D)
 *                    corners 0,1,2,3 (lower-left,lower-right,upper-left,upper-right-2D)
 *                    corners 0,1,2,3,4,5,6,7 (llb,lrb,ulb,urb,llf,lrf,ulf,urf-3D)
 * node_to_position:  given a node number calculate its position.
 * position_to_node:  given a position, find the node number.
 * map_0th_plane:     for 3D problems, map the node or element # to the node or
 *                    element # in the 0th x1-x2 plane.
 * node_to_ijk:       given a node, find the i,j, and k indexes in the x1,x2,x3
 *                    directions.
 * ijk_to_node:       given the ijk indexes, find the node number.
 */

#include "dft_mesh_lib.h"

/****************************************************************************/
/* round_to_int:  Take a double precision number and round to the nearest
                  Integer value       */
int round_to_int(double x)
{
double temp;
 
temp = floor(x);
if ((x - temp) < 0.5)
   return( (int) temp );
else
   return( (int) ceil(x) );
}

/****************************************************************************/
/* offset_to_node:   Given an offset from the stencil, and the current node
                     in terms of ijk, find the node indicated by the stencil*/ 

/* NOTE: This same logic is in 2 places in function load_Jacobian. The code*/
/* runs a bunch faster unrolling this routine into that one. So any changes*/
/* made here should be repeated twice in load_Jacobian                     */
/* reflect_flag tags each dimension which reflected, for HW_Boundary       */

int offset_to_node(int *inode_ijk, int *offset_node,
                   int *reflect_flag)
{
   int i, inode_sten[3];

   for (i=0; i<Ndim; i++) {
      inode_sten[i] = inode_ijk[i] + offset_node[i];
      reflect_flag[i] = FALSE;

      /* the following lines used to be function "boundary_condition" */
      /* This logic is repeated twice in load_Jacobian */
      if (inode_sten[i] < 0) {
        switch (Type_bc[i][0]) {
          case IN_BULK:  
               if ( (Nwall == 0 && Iliq_vap == 3) || 
                     Lsteady_state == TRUE) return(-3);
               else                         return(-1);                     
          case IN_WALL:                     return(-2);                       
          case PERIODIC: 
               inode_sten[i] += Nodes_x[i]; break; 
          case REFLECT:   
               reflect_flag[i] = TRUE;
               inode_sten[i] = -inode_sten[i];
               if (inode_sten[i] >= Nodes_x[i]) {
                   if      (Type_bc[i][1] == IN_WALL) return(-2);
                   else if (Type_bc[i][1] == IN_BULK) return(-1);
                   else {
                     printf("ERROR: You are reflecting beyond the domain\n");
                     exit(-1);
                   }
               }
          case LAST_NODE:
          case LAST_NODE_RESTART:
               inode_sten[i]=0;break;
        }
      }
      else if (inode_sten[i] >= Nodes_x[i]) {
        switch (Type_bc[i][1]) {
          case IN_BULK:  
               if ( (Nwall == 0 && Iliq_vap == 3) ||
                     Lsteady_state == TRUE)    return(-4);
               else                            return(-1);                       
          case IN_WALL:                        return(-2);                        
          case PERIODIC: 
               inode_sten[i] -= Nodes_x[i]; break; 
          case REFLECT:   
               reflect_flag[i] = TRUE;
               inode_sten[i] = 2*(Nodes_x[i]-1)-inode_sten[i]; 
               if (inode_sten[i] < 0) {
                  if      (Type_bc[i][0] == IN_WALL) return(-2);
                  else if (Type_bc[i][0] == IN_BULK) return(-1);
                  else {
                     printf("ERROR: You are reflecting beyond the domain\n");
                     exit(-1);
                  }
               }
          case LAST_NODE:
          case LAST_NODE_RESTART:
               inode_sten[i]=Nodes_x[i]-1;break;
        }
      }
   }
   
   /* convert back to a global node number */
   /* It would be slightly faster to inline that routine here */
   return(ijk_to_node(inode_sten));

}
/****************************************************************************/
/* ijk_to_ijk_box: converts ijk to ijk_box, with special attention to
                   periodic boundary conditions */
void ijk_to_ijk_box(int *ijk, int *ijk_box)
{
  int idim;

  for (idim=0; idim<Ndim; idim++) {
    ijk_box[idim] = ijk[idim] - Min_IJK_box[idim];

    if (Type_bc[idim][0] == PERIODIC && ijk_box[idim] < 0)
      ijk_box[idim] += Nodes_x[idim];
    if (Type_bc[idim][1] == PERIODIC && ijk_box[idim] >= Nodes_x_box[idim])
      ijk_box[idim] -= Nodes_x[idim];
  }
}
/***************************************************************************/
/* ijk_to_ijk_box_no_bound: converts ijk to ijk_box, with no special attention to
                   periodic boundary conditions */
void ijk_to_ijk_box_no_bound(int *ijk, int *ijk_box)
{
  int idim;

  for (idim=0; idim<Ndim; idim++) 
    ijk_box[idim] = ijk[idim] - Min_IJK_box[idim];

  return;
}
/***************************************************************************/
/* ijk_box_to_ijk: converts ijk_box to ijk, with special attention to
                   periodic boundary conditions */
void ijk_box_to_ijk(int *ijk_box, int *ijk)
{
  int idim;

  for (idim=0; idim<Ndim; idim++) {
    ijk[idim] = ijk_box[idim] + Min_IJK_box[idim];

    if (Type_bc[idim][0] == PERIODIC && ijk[idim] < 0)
       ijk[idim] += Nodes_x[idim];
    if (Type_bc[idim][1] == PERIODIC && ijk[idim] >= Nodes_x[idim])
       ijk[idim] -= Nodes_x[idim];
  }
}
/***************************************************************************/
/* offset_to_node_box:   Given an offset from the stencil, and the current node
                     in terms of ijk, find the node indicated by the stencil*/ 

/* reflect_flag tags each dimension which reflected, for HW_Boundary       */
/* since periodic BC are included seemlessly in the box units,
   no special attentian needed for periodic BC */
/* Since the boxes are made big enough to cover all stencils except when
   they would pass the edge of the domain, it is good enough to check
   if the stencil point is outside the box to know if it outside the
   global domain */

int offset_to_node_box(int *ijk_box, int *offset,
                        int *reflect_flag)
{
   int i, ijk_sten_box[3];

   for (i=0; i<Ndim; i++) {
      ijk_sten_box[i] = ijk_box[i] + offset[i];
      reflect_flag[i] = FALSE;

      /* the following lines used to be function "boundary_condition" */
      /* This logic is repeated twice in load_Jacobian */
      if (ijk_sten_box[i] < 0) {
        switch (Type_bc[i][0]) {
          case IN_BULK:
               if ( (Nwall == 0 && Iliq_vap == 3) ||
                     Lsteady_state == TRUE)        return(-3);
               else                                return(-1);
          case IN_WALL:                            return(-2);
          case REFLECT:
               reflect_flag[i] = TRUE;
               ijk_sten_box[i] = -ijk_sten_box[i];
               if  (ijk_sten_box[i] >= Nodes_x_box[i] && 
                      Max_IJK_box[i] == Max_IJK[i]) {
                 switch (Type_bc[i][1]) {
                    case IN_BULK:  
                         if ( (Nwall == 0 && Iliq_vap == 3) ||
                               Lsteady_state == TRUE)        return(-4);
                         else                                return(-1);  
                    case IN_WALL:                            return(-2); 
                    case REFLECT:   
                        printf("Problems with multiple reflections\n");
                        printf("(1) Check domain size and maximum stencil size \n");
                        printf("ijk_sten_box[i=%d]=%d ijk_box=%d offset=%d Nodes_x_box=%d Max_IJK_box=%d  Max_IJK=%d\n",i,ijk_sten_box[i],ijk_box[i],offset[i],Nodes_x_box[i],Max_IJK_box[i],Max_IJK[i]);
                        exit(-1);
                 }
               }
               break;
          case PERIODIC:
               if (Pflag[i]) ijk_sten_box[i] += Nodes_x[i]; break; 
          case LAST_NODE:
          case LAST_NODE_RESTART:
               ijk_sten_box[i]=0;break;
        }
      }

      else if (ijk_sten_box[i] >= Nodes_x_box[i]) {
        switch (Type_bc[i][1]) {
          case IN_BULK:  
               if ( (Nwall == 0 && Iliq_vap == 3) ||
                     Lsteady_state == TRUE)        return(-4);
               else                                return(-1);                       
          case IN_WALL:                            return(-2);                        
          case REFLECT:   
               reflect_flag[i] = TRUE;
               ijk_sten_box[i] = 2*(Nodes_x_box[i]-1)-ijk_sten_box[i]; 
               if (ijk_sten_box[i] < 0 && Min_IJK_box[i] == Min_IJK[i]) {
                  switch (Type_bc[i][0]) {
                    case IN_BULK:
                    if ( (Nwall == 0 && Iliq_vap == 3) ||
                         Lsteady_state == TRUE)        return(-3);
                    else                               return(-1);
                    case IN_WALL:                      return(-2);
                    case REFLECT:
                        printf("Problems with multiple reflections\n");
                        printf("(2) Check domain size and maximum stencil size \n");
                        printf("ijk_sten_box[i=%d]=%d (min=0) Min_IJK_box=%d  Min_IJK=%d\n",i,ijk_sten_box[i],Min_IJK_box[i],Min_IJK[i]);
                        exit(-1);
                  }
               }
               break;
          case PERIODIC:
               if (Pflag[i]) ijk_sten_box[i] -= Nodes_x[i]; break; 
          case LAST_NODE:
          case LAST_NODE_RESTART:
               ijk_sten_box[i]=Nodes_x_box[i]-1;break;
        }
      }
   }
   
   /* convert back to node number in extended local domain */
   /* It is slightly faster to inline ijk_box_to_node_box here */
   /* return(ijk_box_to_node_box(ijk_sten_box)); */

  switch (Ndim){
    case 1:   return (*ijk_sten_box);
    case 2:   return (*ijk_sten_box + ijk_sten_box[1] * (*Nodes_x_box));
    case 3:   return (*ijk_sten_box + ijk_sten_box[1] * (*Nodes_x_box)
                              + ijk_sten_box[2] * Nodes_plane_box);
  }
  return 0;
}
/****************************************************************************/
/* element_to_node:   Given to element, return an associated node. It will
                      be the node to the left(1D), the lower left corner (2D), 
                      or the lower left back corner (3D).   */ 
int element_to_node(int ielement)
{
  int inode=0, iplane, Element_0;

  if (Ndim == 1)
     inode = ielement;
  else if (Ndim == 2)
     inode = ielement + (ielement/Elements_x[0])*(Nodes_x[0]-Elements_x[0]);
  else if (Ndim == 3){
     iplane = ielement/Elements_plane;
     Element_0 = map_0th_plane(ielement,Elements_plane);
     inode = iplane*Nodes_plane + Element_0
             + (Element_0/Elements_x[0])*(Nodes_x[0]-Elements_x[0]);;
  }
  return inode;
}
/****************************************************************************/
/* element_to_nodes:  Given to element, return all of the associated nodes. */
void element_to_nodes(int ielement,int *nodes)
{
  int inode, idim,iplane, Element_0,ijk[3];

  if (Ndim == 1) {
     nodes[0] = ielement;
     nodes[1] = nodes[0]+1;
  }
  else if (Ndim == 2) {
     nodes[0] = ielement + (ielement/Elements_x[0])*(Nodes_x[0]-Elements_x[0]);
     nodes[1] = nodes[0]+1;
     nodes[2] = nodes[0]+Nodes_x[0];
     nodes[3] = nodes[0]+Nodes_x[0]+1;
  }
  else if (Ndim == 3){
     iplane = ielement/Elements_plane;
     Element_0 = map_0th_plane(ielement,Elements_plane);
     nodes[0] = iplane*Nodes_plane + Element_0
             + (Element_0/Elements_x[0])*(Nodes_x[0]-Elements_x[0]);;
     nodes[1] = nodes[0]+1;
     nodes[2] = nodes[0]+Nodes_x[0];
     nodes[3] = nodes[0]+Nodes_x[0]+1;
     nodes[4] = nodes[0]+Nodes_plane;
     nodes[5] = nodes[0]+Nodes_plane+1;
     nodes[6] = nodes[0]+Nodes_plane+Nodes_x[0];
     nodes[7] = nodes[0]+Nodes_plane+Nodes_x[0]+1;
  }
  for (inode=0; inode<Nnodes_per_el_V; inode++){
     node_to_ijk(nodes[inode],ijk);
     for (idim=0; idim<Ndim; idim++){
        if (Type_bc[idim][1]==PERIODIC && ijk[idim] >= Nodes_x[idim]){
            ijk[idim] = ijk[idim]-Nodes_x[idim];
            nodes[inode]=ijk_to_node(ijk);
        }
     }
  }
  return;
}
/****************************************************************************/
int node_to_elem(int inode_all, int local_node, int *reflect_flag)
/*
 * node_to_elem gives the element that has node inode_all as
 * its local node local_node
 * local node numbering given by "set_offsets_and_weights"
 */
{
  int idim, ijk[NDIM_MAX], iel;

  node_to_ijk(inode_all, ijk);
  iel = 0;

  /*
   * ijk[Ndim] is the offsets in each dimension to the element
   * that has this node as node 0. Now adjust based on local node #.
   */

  if ((local_node/2 == (local_node+1)/2) == reflect_flag[0])   (ijk[0])--;

  if (Ndim > 1)
    if ((local_node-4*(local_node/4) < 2) == reflect_flag[1])  (ijk[1])--;

  if (Ndim == 3)
    if ((local_node < 4) == reflect_flag[2])                   (ijk[2])--;

  /*
   * Now adjust if nodal values are out of range
   * Error checking should be removed in fast production version.
   */
  for (idim=0; idim<Ndim; idim++) {
    if (ijk[idim] == -1) {
       if (Type_bc[idim][0] == PERIODIC) ijk[idim] = Elements_x[idim] - 1;
       else if (Type_bc[idim][0] == REFLECT) ijk[idim] = 0;
       else if (Type_bc[idim][0] == IN_BULK)  iel = -2;
       else if (Type_bc[idim][0] == IN_WALL)  iel = -1;
       else if (Type_bc[idim][0] == LAST_NODE || Type_bc[idim][0] == LAST_NODE_RESTART) ijk[idim] = 0;
    }
    else if (ijk[idim] == Elements_x[idim]) {
      if (Type_bc[idim][1] == PERIODIC) ijk[idim] = 0; 
      else if (Type_bc[idim][1] == REFLECT) ijk[idim] = Elements_x[idim]-1; 
      else if (Type_bc[idim][1] == IN_BULK)  iel = -2;
      else if (Type_bc[idim][1] == IN_WALL)  iel = -1;
      else if (Type_bc[idim][1] == LAST_NODE || Type_bc[idim][1] == LAST_NODE_RESTART) ijk[idim] = Elements_x[idim]-1;
    }
  }

  /*
   * Now, compute element number.
   */
  if (iel >= 0 ){
     iel = ijk[0];
     if (Ndim > 1)  iel += ijk[1] * Elements_x[0];
     if (Ndim == 3) iel += ijk[2] * Elements_plane;
  }

  return(iel);
}
/****************************************************************************/
int node_to_elem_return_dim(int inode_all, int local_node, int *reflect_flag,
                             int *idim_return, int *iside, int *periodic_flag)
/*
 * node_to_elem gives the element that has node inode_all as
 * its local node local_node
 * local node numbering given by "set_offsets_and_weights"
 */
{
  int idim, ijk[NDIM_MAX], iel;

  node_to_ijk(inode_all, ijk);
  iel = 0;

  /*
   * ijk[Ndim] is the offsets in each dimension to the element
   * that has this node as node 0. Now adjust based on local node #.
   */

  if ( (local_node/2 != (local_node+1)/2) && (reflect_flag[0]==FALSE) )
    {
      (ijk[0])--;
    }
  else if ( (local_node/2 == (local_node+1)/2) && (reflect_flag[0]==TRUE) )
    {
      (ijk[0])--;
    }

  if (Ndim > 1) 
    {
      if ( (local_node-4*(local_node/4) >= 2) && (reflect_flag[1]==FALSE) )
	{
	  (ijk[1])--;
	}
      else if ( (local_node-4*(local_node/4) < 2) && (reflect_flag[1]==TRUE) )
	{
	  (ijk[1])--;
	}
    }

  if (Ndim == 3)
    {
      if ( (local_node >= 4) && (reflect_flag[2]==FALSE) )
	{
	  (ijk[2])--;
	}
      else if ( (local_node < 4) && (reflect_flag[2]==TRUE) )
	{
	  (ijk[2])--;
	}
    }

  /*
   * Now adjust if nodal values are out of range
   * Error checking should be removed in fast production version.
   */

  *periodic_flag=FALSE;
  for (idim=0; idim<Ndim; idim++) {
    if (ijk[idim] == -1) {
       if (Type_bc[idim][0] == PERIODIC) {
           ijk[idim] = Elements_x[idim] - 1;
           *periodic_flag = TRUE;
       }
       else if (Type_bc[idim][0] == REFLECT) ijk[idim] = 0;
       else if (Type_bc[idim][0] == IN_BULK)  iel = -2;
       else if (Type_bc[idim][0] == IN_WALL)  {
           iel = -1; *idim_return = idim; *iside = 0;
       }
       else if (Type_bc[idim][0] == LAST_NODE || Type_bc[idim][0] == LAST_NODE_RESTART) ijk[idim]=0; 
    }
    else if (ijk[idim] == Elements_x[idim]) {
      if (Type_bc[idim][1] == PERIODIC) {
          ijk[idim] = 0;
          *periodic_flag=TRUE;
      }
      else if (Type_bc[idim][1] == REFLECT) ijk[idim] = Elements_x[idim]-1;
      else if (Type_bc[idim][1] == IN_BULK)  iel = -2;
      else if (Type_bc[idim][1] == IN_WALL) {
           iel = -1; *idim_return = idim; *iside = 1;
      }
      else if (Type_bc[idim][1] == LAST_NODE || Type_bc[idim][1] == LAST_NODE_RESTART) ijk[idim]=Elements_x[idim]-1; 
    }
  }

  /*
   * Now, compute element number.
   */
  if (iel >= 0){
     iel = ijk[0];
     if (Ndim > 1)  iel += ijk[1] * Elements_x[0];
     if (Ndim == 3) iel += ijk[2] * Elements_plane;
  }

  return(iel);
}
/****************************************************************************/
int node_to_elem_v2(int inode_all, int local_node)
/*
 * node_to_elem gives the element that has node inode_all as
 * its local node local_node
 * local node numbering given by "set_offsets_and_weights"
 * version 2 returns flag for reflections instead of reflected element.
 */
{
  int idim, ijk[NDIM_MAX], iel;

  node_to_ijk(inode_all, ijk);
  iel = 0;

  /*
   * ijk[Ndim] is the offsets in each dimension to the element
   * that has this node as node 0. Now adjust based on local node #.
   */

  if (local_node/2 != (local_node+1)/2)     (ijk[0])--;

  if (Ndim > 1) 
    if (local_node-4*(local_node/4) >= 2)    (ijk[1])--;

  if (Ndim == 3)
    if (local_node >= 4)                   (ijk[2])--;

  /*
   * Now adjust if nodal values are out of range, which should
   * only happen for periodic BC's. (There shouldn't be a rho_bulk
   * BC in an element touching a hard wall.)
   *
   * Error checking should be removed in fast production version.
   */

  for (idim=0; idim<Ndim; idim++) {
    if (ijk[idim] == -1) {
       if (Type_bc[idim][0] == PERIODIC) ijk[idim] = Elements_x[idim] - 1;
       else if (Type_bc[idim][0] == REFLECT)  iel = -2;
       else if (Type_bc[idim][0] == IN_BULK)  iel = -2;
       else if (Type_bc[idim][0] == IN_WALL)  iel = -1;
       else if (Type_bc[idim][0] == LAST_NODE || Type_bc[idim][0] == LAST_NODE_RESTART)  iel = -2;
    }
    else if (ijk[idim] == Elements_x[idim]) {
      if (Type_bc[idim][1] == PERIODIC) ijk[idim] = 0;
      else if (Type_bc[idim][1] == REFLECT)  iel = -2;
      else if (Type_bc[idim][1] == IN_BULK)  iel = -2;
      else if (Type_bc[idim][1] == IN_WALL)  iel = -1;
      else if (Type_bc[idim][1] == LAST_NODE || Type_bc[idim][1] == LAST_NODE_RESTART)  iel = -2;
    }
  }

  /*
   * Now, compute element number.
   */
  if (iel >= 0){
     iel = ijk[0];
     if (Ndim > 1)  iel += ijk[1] * Elements_x[0];
     if (Ndim == 3) iel += ijk[2] * Elements_plane;
  }

  return(iel);
}
/****************************************************************************/
/* node_to_position:  Calculate the position of a node based on the node
                      number.  Note, this node number must be a reference
                      to a rectangular grid in the entire domain */ 
void node_to_position(int inode, double *NodePos)
{
  int Node_0th_plane;

  /* write column headings to file: for checking routine */

/*
 *fprintf(fp2,"\n   ********* TEST THE MESH GENERATION ********\n");
 *if (Ndim ==3){
 *   fprintf(fp2,"\tinode\tNodePos[0]\tNodePos[1]\tNodePos[2]\tinode_calc\n");
 *}
 *else if (Ndim == 2) {
 *   fprintf(fp2,"\tinode \tNodePos[0] \tNodePos[1] \tinode_calc\n");
 *}
 *else if (Ndim ==1) {
 *   fprintf(fp2,"\tinode \tNodePos[0] \tinode_calc\n"); 
 *} 
 */
  /* fill an array of the positions of the nodes */

  if (Ndim == 3) {

     Node_0th_plane = map_0th_plane(inode, Nodes_plane);

          
     NodePos[0] = -0.5*Size_x[0] + Esize_x[0]*( Node_0th_plane 
                  - Nodes_x[0]*(Node_0th_plane/Nodes_x[0]) ) ;

     NodePos[1] = -0.5*Size_x[1] + Esize_x[1]*(
                                 Node_0th_plane/Nodes_x[0] );

     NodePos[2] = -0.5*Size_x[2] + Esize_x[2]*( inode /
                                               Nodes_plane );
  }
  else if (Ndim == 2){
     NodePos[0] = -0.5*Size_x[0] + Esize_x[0]*( inode -
                        (inode/Nodes_x[0])*Nodes_x[0] );

     NodePos[1] = -0.5*Size_x[1] + Esize_x[1]*( inode /
                                           Nodes_x[0] );
  }
  else if (Ndim == 1){
     NodePos[0] = -0.5*Size_x[0] + Esize_x[0]*( inode );
  }

  /* write positions and calculated node number to the file */

/*
 *if (Ndim ==3){
 *   inode_calc = position_to_node(NodePos);
 *   fprintf(fp2,"%7d \t%8.4f \t%8.4f \t%8.4f \t%7d",inode,NodePos[0],
 *                                  NodePos[1],NodePos[2],inode_calc); 
 *}
 *else if (Ndim == 2) {
     inode_calc = position_to_node(NodePos);
 *   fprintf(fp2,"%7d \t%8.4f \t%8.4f \t%7d",inode,NodePos[0],
 *                                     NodePos[1],inode_calc); 
 *}
 *else if (Ndim ==1) {
     inode_calc = position_to_node(NodePos);
     fprintf(fp2,"%7d \t%8.4f \t%7d\n",inode,NodePos[0],inode_calc); 
 *} 
 */

}
/****************************************************************************/

/* position_to_node:  Calculate a node number based on the position of the
		      node.  Use to check the reverse calculation !!! */
int position_to_node(double *NodePos)
{
   int inode_calc;

   if (Ndim == 3){
      inode_calc = ((NodePos[2]+0.5*Size_x[2])/Esize_x[2]) * Nodes_x[0] * Nodes_x[1]
                          + ((NodePos[1]+0.5*Size_x[1])/Esize_x[1])* Nodes_x[0] 
                             + (NodePos[0] + 0.5*Size_x[0])/Esize_x[0] + 0.001;
      return inode_calc;
   }   
   else if (Ndim == 2){
      inode_calc = ((NodePos[1]+0.5*Size_x[1])/Esize_x[1])*Nodes_x[0] 
                              + (NodePos[0] + 0.5*Size_x[0])/Esize_x[0] + 0.001;
      return inode_calc;
   }   
   else if (Ndim == 1){
      inode_calc =  (NodePos[0] + 0.5*Size_x[0])/Esize_x[0] + 0.001;
      return inode_calc;
   }   
  return 0;
}
/****************************************************************************/
/* map_0th_plane:  Given an element or node 
                   return the corresponding element or node in the 0th plane */
int map_0th_plane(int i, int Nplane)
{
    int mapping = 0;
    if (Ndim ==2) 
         mapping = i;
    else if (Ndim == 3) 
         mapping = i - (i/Nplane)*Nplane;
    return mapping;
}
/****************************************************************************/
/* node_to_ijk:  Given a node number, find the i,j,k node numbers, the offsets
                   into each dimension */
void node_to_ijk(int node, int *ijk)
{

  switch (Ndim) {
     /* If 3D, put the plane number in ijk[2] */
     case 3:     ijk[2] = node/(Nodes_plane);
                 node -= ijk[2] * Nodes_plane;

     /* Now in a plane, put the line number in ijk[1] */
     case 2:     ijk[1] = node/(*Nodes_x);
                 node -= ijk[1] * (*Nodes_x);

     /* Now its reduced to 1D, which is trivial */
     case 1:     ijk[0] = node;
  }
}
/****************************************************************************/
/* node_box_to_ijk_box:  Given a node number, find the i,j,k node numbers, the offsets
                   into each dimension */
void node_box_to_ijk_box(int node_box, int *ijk_box)
{

  switch (Ndim) {
     /* If 3D, put the plane number in ijk[2] */
     case 3:     ijk_box[2] = node_box/(Nodes_plane_box);
                 node_box -= ijk_box[2] * Nodes_plane_box;

     /* Now in a plane, put the line number in ijk[1] */
     case 2:     ijk_box[1] = node_box/(*Nodes_x_box);
                 node_box -= ijk_box[1] * (*Nodes_x_box);

     /* Now its reduced to 1D, which is trivial */
     case 1:     ijk_box[0] = node_box;
  }
}
/****************************************************************************/
/* ijk_to_node:  Given an i,j,k position of a node, return its node number */

int ijk_to_node(int *ijk)
{
  switch (Ndim){
    case 1:   return (*ijk);
    case 2:   return (*ijk + ijk[1] * (*Nodes_x));
    case 3:   return (*ijk + ijk[1] * (*Nodes_x) + ijk[2] * Nodes_plane);
  }
  return 0;
}
/****************************************************************************/
/* ijk_box_to_node_box:  Given an i,j,k position of a node, return its node number */

int ijk_box_to_node_box(int *ijk_box)
{
  switch (Ndim){
    case 1:   return (*ijk_box);
    case 2:   return (*ijk_box + ijk_box[1] * (*Nodes_x_box));
    case 3:   return (*ijk_box + ijk_box[1] * (*Nodes_x_box) 
                              + ijk_box[2] * Nodes_plane_box);
  }
  return 0;
}
/****************************************************************************/
/* node_box_to_node: */

int node_box_to_node(int inode_box)
{
  int ijk_box[3], ijk[3];
  node_box_to_ijk_box(inode_box, ijk_box);
  ijk_box_to_ijk(ijk_box, ijk);
  return(ijk_to_node(ijk));
}
/****************************************************************************/
/* node_to_node_box: */

int node_to_node_box(int inode)
{
  int ijk_box[3], ijk[3],inode_box;
  node_to_ijk(inode,ijk);
  ijk_to_ijk_box(ijk, ijk_box);
  inode_box =ijk_box_to_node_box(ijk_box);
  return(inode_box);
}
/****************************************************************************/
/* node_to_node_box_no_bound: don't do boundary checking*/

int node_to_node_box_no_bound(int inode)
{
  int ijk_box[3], ijk[3],inode_box,idim;

  node_to_ijk(inode,ijk);
  ijk_to_ijk_box_no_bound(ijk,ijk_box);

  for (idim=0; idim<Ndim; idim++)
     if(ijk_box[idim] < 0 || ijk_box[idim] >= Nodes_x_box[idim])
        return(-1);

  inode_box =ijk_box_to_node_box(ijk_box);
  return(inode_box);
}
/****************************************************************************/
/* unk_box_to_unk: */

int unk_box_to_unk(int i_box)
{
   int inode_box, inode, icomp;

   inode_box = i_box/Nunk_per_node;
   icomp = i_box - inode_box*Nunk_per_node;
   inode = node_box_to_node(inode_box);
   return(inode*Nunk_per_node + icomp);
}
/****************************************************************************/
/* unk_to_unk_box: */

int unk_to_unk_box(int i)
{
   int inode_box, inode, icomp;

   inode = i/Nunk_per_node;
   icomp = i - inode*Nunk_per_node;
   inode_box = node_to_node_box(inode);
   return(inode_box*Nunk_per_node + icomp);
}
/****************************************************************************/
/* el_box_to_el: */

int el_box_to_el(int iel_box)
{
  int inode_box,inode,iel;
  int reflect_flag[3];

  reflect_flag[0] = reflect_flag[1] = reflect_flag[2] = FALSE;
  inode_box = element_box_to_node_box (iel_box);
  inode = node_box_to_node (inode_box);

  /*always look for element w/ 0th node at inode */
  iel = node_to_elem(inode,0,reflect_flag);
  return(iel);
}
/****************************************************************************/
/* element_box_to_node_box:   
        Given to element, return an associated node. It will
        be the node to the left(1D), the lower left corner (2D), 
        or the lower left back corner (3D).   */ 
int element_box_to_node_box(int iel_box)
{
  int iplane, Element_0,inode_box=0;

  if (Ndim == 1)
     inode_box = iel_box;
  else if (Ndim == 2)
     inode_box = iel_box + (iel_box/Elements_x_box[0])*
                           (Nodes_x_box[0]-Elements_x_box[0]);
  else if (Ndim == 3){
     iplane = iel_box/Elements_plane_box;
     Element_0 = map_0th_plane(iel_box,Elements_plane_box);
     inode_box = iplane*Nodes_plane_box + Element_0
             + (Element_0/Elements_x_box[0])*
               (Nodes_x_box[0]-Elements_x_box[0]);;
  }
  return inode_box;
}
/****************************************************************************/
/* el_to_el_box: */

int el_to_el_box(int iel)
{
  int inode_box,inode,ijk[3],idim,ijk_box[3],iel_box;

  inode = element_to_node(iel);

  node_to_ijk(inode,ijk);
  ijk_to_ijk_box(ijk, ijk_box);

  for (idim=0; idim<Ndim; idim++){
      if (ijk_box[idim] > Elements_x_box[idim]-1){
         if (Type_bc[idim][1] == PERIODIC) ijk_box[idim] -= Nodes_x[idim];
         else return(-1);
      }
/*      if (ijk_box[idim] < 0) return(-1);*/
      else if (ijk_box[idim] < 0) {
         if (Type_bc[idim][0] == PERIODIC) ijk_box[idim] += Nodes_x[idim];
         else return(-1);
      }
  }

  inode_box = ijk_box_to_node_box(ijk_box);

  /*always look for element w/ 0th node at inode_box */
  iel_box = node_box_to_elem_box(inode_box,0); 

  return(iel_box);
}
/****************************************************************************/
int node_box_to_elem_box(int inode_box, int local_node)
/*
 * node_box_to_elem_box: gives the element that touches inode_box 
 *                       at the local_node corner.  All in box units.
 * we don't check reflections in this routine .... this may be needed.
 */
{
  int ijk_box[NDIM_MAX], iel_box;

  node_box_to_ijk_box(inode_box, ijk_box);

  /*
   * ijk_box[Ndim] is the offsets in each dimension to the element
   * that has this node as node 0. Now adjust based on local node #.
   */

  if (local_node/2 != (local_node+1)/2)       (ijk_box[0])--;
  if (Ndim > 1) 
    if (local_node-4*(local_node/4) >= 2 )    (ijk_box[1])--;
  if (Ndim == 3)
    if (local_node >= 4)                      (ijk_box[2])--;

  /*
   * Now, compute element number.
   */

  iel_box = ijk_box[0];
  if (Ndim > 1)  iel_box += ijk_box[1] * Elements_x_box[0];
  if (Ndim == 3) iel_box += ijk_box[2] * Elements_plane_box;

  return(iel_box);
}
/****************************************************************************/
int node_box_to_elem_box_reflect(int inode_box, int local_node, int *reflect_flag)
/*
 * node_box_to_elem_box: gives the element that touches inode_box
 *                       at the local_node corner.  All in box units.
 */
{
  int ijk_box[3], iel_box=0;


  node_box_to_ijk_box(inode_box, ijk_box);

  /*
   * ijk_box[Ndim] is the offsets in each dimension to the element
   * that has this node as node 0. Now adjust based on local node #.
   */

  switch (Ndim) {

     case 3:
       if ((local_node < 4) == reflect_flag[2])  {
          (ijk_box[2])--;
          if (ijk_box[2] == -1) {
            if (Type_bc[2][0] == REFLECT) ijk_box[2] = 0;
            else if (Type_bc[2][0] == LAST_NODE || Type_bc[2][0] == LAST_NODE_RESTART) ijk_box[2] = 0;
            else if (Type_bc[2][0] == PERIODIC && Pflag[2]) 
                                          ijk_box[2] = Elements_x[2]-1;
            else     return (-1); /* out of domain -- no contribution */
          }
       }
       else if (ijk_box[2] == Elements_x_box[2]) {
            if (Type_bc[2][1] == REFLECT) ijk_box[2] = Elements_x_box[2]-1;
            else if (Type_bc[2][1] == LAST_NODE || Type_bc[2][1] == LAST_NODE_RESTART) ijk_box[2] = Elements_x_box[2]-1;
            else if (Type_bc[2][1] == PERIODIC && Pflag[2]) ijk_box[2] = 0;
            else     return (-1); /* out of domain -- no contribution */
       }
       iel_box += ijk_box[2] * Elements_plane_box;

     case 2:
       if ((local_node%4 < 2 ) == reflect_flag[1]) {
          (ijk_box[1])--;
          if (ijk_box[1] == -1) {
            if (Type_bc[1][0] == REFLECT) ijk_box[1] = 0;
            else if (Type_bc[1][0] == LAST_NODE || Type_bc[1][0] == LAST_NODE_RESTART) ijk_box[1] = 0;
            else if (Type_bc[1][0] == PERIODIC && Pflag[1]) 
                                          ijk_box[1] = Elements_x[1]-1;
            else     return (-1); /* out of domain -- no contribution */
          }
       }
       else if (ijk_box[1] == Elements_x_box[1]) {
            if (Type_bc[1][1] == REFLECT) ijk_box[1] = Elements_x_box[1]-1;
            else if (Type_bc[1][1] == LAST_NODE || Type_bc[1][1] == LAST_NODE_RESTART) ijk_box[1] = Elements_x_box[1]-1;
            else if (Type_bc[1][1] == PERIODIC && Pflag[1]) ijk_box[1] = 0;
            else     return (-1); /* out of domain -- no contribution */
       }
       iel_box += ijk_box[1] * Elements_x_box[0];

     case 1:
       if ((local_node%2 == 0)  == reflect_flag[0]) {
          (ijk_box[0])--;
          if (ijk_box[0] == -1) {
            if (Type_bc[0][0] == REFLECT) ijk_box[0] = 0;
            else if (Type_bc[0][0] == LAST_NODE || Type_bc[0][0] == LAST_NODE_RESTART) ijk_box[0] = 0;
            else if (Type_bc[0][0] == PERIODIC && Pflag[0]) 
                                          ijk_box[0] = Elements_x[0]-1;
            else     return (-1); /* out of domain -- no contribution */
          }
       }
       else if (ijk_box[0] == Elements_x_box[0]) {
            if (Type_bc[0][1] == REFLECT) ijk_box[0] = Elements_x_box[0]-1;
            else if (Type_bc[0][1] == LAST_NODE || Type_bc[0][1] == LAST_NODE_RESTART) ijk_box[0] = Elements_x_box[0]-1;
            else if (Type_bc[0][1] == PERIODIC && Pflag[0]) ijk_box[0] = 0;
            else     return (-1); /* out of domain -- no contribution */
       }
       iel_box += ijk_box[0];
  }

  return(iel_box);
}
/****************************************************************************/
