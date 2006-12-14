/*
//@HEADER
// ******************************************************************** 
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/
/*
 *  FILE: dft_fill_coarse.c
 *
 *  This file contains the routines to fill coarsened residual equations
 *  into the matrix problem as averages of surrounding nodal values
 */

/*#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"*/
#include "dft_fill_coarse.h"

/********************************************************************/
void load_coarse_node_1dim(int loc_inode, int inode_box, int *ijk_box,int iunk, double **x)
{
    int idim,ijk_tmp[3],jnode_box,loc_jnode,numEntries,
        nodeIndices[2];
    double resid,values[2];

    for (idim=0; idim<Ndim; idim++) ijk_tmp[idim]=0;
    ijk_tmp[Grad_dim]=ijk_box[Grad_dim];
    jnode_box = ijk_box_to_node_box(ijk_tmp);
    loc_jnode = B2L_node[jnode_box]; 

    if (jnode_box <0 ){
        printf("Proc: %d ijk_box: %d %d %d PROBLEMS: jnode_box: %d  ijk_tmp: %d %d %d\n",
        Proc,ijk_box[0],ijk_box[1],ijk_box[2],jnode_box,ijk_tmp[0],ijk_tmp[1],ijk_tmp[2]);
        exit(-1);
    }

    resid = x[iunk][inode_box]-x[iunk][jnode_box];
    dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
    numEntries=2;
    values[0]=1.0; values[1]=-1.0;
    nodeIndices[0]=inode_box; nodeIndices[1]=jnode_box;
    dft_linprobmgr_insertmultinodematrixvalues(LinProbMgr_manager,iunk,loc_inode,
                                               iunk, nodeIndices,values,numEntries);
    return;
}
/********************************************************************/
double load_coarse_node_Ndim(int loc_inode, int inode_box, int iunk, double **x)
{
    double fac_a11,resid,mat_value,fac_coarse;

    if (Unk2Phys[iunk] ==DENSITY || POISSON) fac_a11=1.0; /* temporary factor to keep -1 on diagonal of A11 block */
    else fac_a11=-1.0;

    resid= fac_a11*x[iunk][inode_box];
    mat_value=fac_a11*1.0;
    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_value);

    fac_coarse=0.5*fac_a11;
    locate_neighbor_unks(x,iunk,loc_inode,inode_box,fac_coarse,&resid);
    dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);

    return(resid);
}
/****************************************************************************/
double load_coarse_variable(double **x,int jnode_box,double fac,int iunk,int loc_inode)
{
  double mat_value,resid;
  if (jnode_box<0) resid = -fac*constant_boundary(iunk,jnode_box);
  else {                 
      resid = -fac*x[iunk][jnode_box];
      mat_value=-fac;
      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,jnode_box,mat_value);
  }
  return resid;
}
/****************************************************************************/
void locate_neighbor_unks(double **x,int iunk, int loc_inode, int node_box,double fac,double *resid)
{
  /* identify the neighboring nodes that contribute to a particular coarsened node.  For
     the special case of only two zones (one possible coarsening level), the recursive code
     casts the coarsened node values in terms of independent nodes only */ 

                                   /* the 6 offset patterns for nearest neighbors */
  int offset_idim_pm[18] = {1,0,0,  0,1,0,  0,0,1,  -1,0,0,  0,-1,0,  0,0,-1};
  int *offset_ptr; /* pointer into above */
  int ijk_box[3],loop,jnode_box;
  int     reflect_flag[3];

   node_box_to_ijk_box(node_box,ijk_box);
   for (loop=0;loop<2;loop++){
       if (loop==0) offset_ptr = &offset_idim_pm[3*(-Mesh_coarsen_flag[node_box] - 1)];  /* one node higher */
       else         offset_ptr += 9;                                             /* one node lower */
       jnode_box = offset_to_node_box(ijk_box, offset_ptr, reflect_flag);

       if (jnode_box >=0 && Mesh_coarsen_flag[jnode_box]<0 && Nzone == 2){
           locate_neighbor_unks(x,iunk,loc_inode,jnode_box,fac*0.5,resid);
       }
       else *resid += load_coarse_variable(x,jnode_box,fac,iunk,loc_inode);
   }
   return;
}
/****************************************************************************/
