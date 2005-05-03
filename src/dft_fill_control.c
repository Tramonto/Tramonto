
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
void fill_resid_and_matrix_control (double *x, double *resid,
                            int **bindx_2d, double *fill_time, int fill_flag,
                            int iter, int resid_only_flag)
{
   int i;
   if (MATRIX_FILL_NODAL){
        if (resid_only_flag==FALSE){
           if (fill_flag==MSR_PREPROCESS){
              /* For MSR pre-processing, start nonzeros counter at # of diagonals */ 
              First_time = TRUE;
              Aztec.nonzeros = Aztec.N_update;
              if (Type_poly != -1)
                 fill_resid_and_matrix_P(NULL, NULL, bindx_2d, NULL, fill_flag, iter, resid_only_flag,NODAL_FLAG);
              else
                 fill_resid_and_matrix(NULL, NULL, bindx_2d, NULL, fill_flag, iter, resid_only_flag,NODAL_FLAG);
           }
           else{
              if (Type_poly != -1)
                 fill_resid_and_matrix_P(x,resid,NULL,fill_time,Matrix_fill_flag,iter,resid_only_flag,NODAL_FLAG);
              else
                 fill_resid_and_matrix(x,resid,NULL,fill_time,Matrix_fill_flag,iter,resid_only_flag,NODAL_FLAG);
           }
        }
        else{
           if (Type_poly != -1)
              fill_resid_and_matrix_P(x,resid,NULL,NULL,Matrix_fill_flag,1,resid_only_flag,NODAL_FLAG);
           else
              fill_resid_and_matrix(x,resid,NULL,NULL,Matrix_fill_flag,1,resid_only_flag,NODAL_FLAG);
        }
   }
   else{
      for (i=0;i<Nunk_per_node;i++){

        if (resid_only_flag==FALSE){
           if (fill_flag==MSR_PREPROCESS){
              if (i==0){
                  /* For MSR pre-processing, start nonzeros counter at # of diagonals */ 
                  First_time = TRUE;
                  Aztec.nonzeros = Aztec.N_update;
              }
              if (Type_poly != -1)
                 fill_resid_and_matrix_P(NULL, NULL, bindx_2d, NULL, fill_flag, iter, resid_only_flag,i);
              else
                 fill_resid_and_matrix(NULL, NULL, bindx_2d, NULL, fill_flag, iter, resid_only_flag,i);
           }
           else{
              if (Type_poly != -1)
                 fill_resid_and_matrix_P(x,resid,NULL,fill_time,Matrix_fill_flag,iter,resid_only_flag,i);
              else{
                 fill_resid_and_matrix(x,resid,NULL,fill_time,Matrix_fill_flag,iter,resid_only_flag,i);
              }
           }
        }
        else{
           if (Type_poly != -1)
              fill_resid_and_matrix_P(x,resid,NULL,NULL,Matrix_fill_flag,1,resid_only_flag,i);
           else
              fill_resid_and_matrix(x,resid,NULL,NULL,Matrix_fill_flag,1,resid_only_flag,i);
        }

      }
   }

}
/*****************************************************************************************************/
