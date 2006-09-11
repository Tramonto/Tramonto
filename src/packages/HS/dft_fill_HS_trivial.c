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
 *  FILE: dft_fill_HS_trivial.c
 *
 *  This file includes all of the trivial parts of the matrix fill associated with
 *  the hard sphere nonlocal density equations. Trivial fills are identity blocks 
 *  and zero blocks.
 *  
 *
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

/*****************************************************************************/
int fill_hsrhobar_hsrhobar(int inode,int iunk,int resid_only_flag,double **x,double resid)
{
  resid =-x[iunk][inode_box];
  mat_val=-1.0;
  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
  dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);

  return (IDENTITY_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_hsrhobar_poisson(int inode,int iunk,int resid_only_flag,double **x,double resid)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_hsrhobar_diffusion(int inode,int iunk,int resid_only_flag,double **x,double resid)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_hsrhobar_cavwtc(int inode,int iunk,int resid_only_flag,double **x,double resid)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_hsrhobar_usrvar1(int inode,int iunk,int resid_only_flag,double **x,double resid)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_hsrhobar_usrvar2(int inode,int iunk,int resid_only_flag,double **x,double resid)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_hsrhobar_usrvar3(int inode,int iunk,int resid_only_flag,double **x,double resid)
{
   return (ZERO_BLOCK_FLAG);
}

