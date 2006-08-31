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
 *  FILE: dft_fill_control.c
 *
 *  This file contains logic to switch between nodal and physics ordering
 *  whe performing the fill_resid_and_matrix_routines...
 *
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

/****************************************************************************/
void fill_resid_and_matrix_control (double **x, int iter, int resid_only_flag)
{
   int i,iter_tmp;
   if (resid_only_flag) iter_tmp=1;
   else iter_tmp=iter;

   if (MATRIX_FILL_NODAL){
        if (Type_poly != NONE && Type_poly != WTC) fill_resid_and_matrix_P(x,iter,resid_only_flag,NODAL_FLAG);
        else                 fill_resid_and_matrix(x,iter,resid_only_flag,NODAL_FLAG);
   }
   else{
      for (i=0;i<Nunk_per_node;i++){

        if (Type_poly != NONE && Type_poly != WTC) fill_resid_and_matrix_P(x,iter_tmp,resid_only_flag,i);
         else                fill_resid_and_matrix(x,iter_tmp,resid_only_flag,i);

      }
   }
}
/*****************************************************************************************************/
