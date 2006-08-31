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

#include "dft_linop.h"

struct DFT_AZTEC_LINOP_STRUCT {
  int ** bindx_2d;
  
};
typedef struct DFT_AZTEC_LINOP_STRUCT DFT_AZTEC_LINOP;

void dft_create_linop(DFT_VECSPACE * vecspace, DFT_LINOP ** linop){

  *linop = (DFT_LINOP *) new DFT_LINOP;
  (*linop)->vecspace = vecspace;
  (*linop)->first_time_graph_insert = TRUE;
  (*linop)->nonzeros = 0;
  int N_loc = vecspace->N_loc;
  if (vecspace->vecspace_type==DFT_AZTEC_VECSPACE) {
    int ** bindx_2d = new int*[N_loc+1];
    for (int i=0; i<N_loc+1; i++) bindx_2d = 0;
    if (N_loc>0) bindx_2d[N_loc] = new int[N_loc];
    DFT_AZTEC_LINOP * aztec_linop = new DFT_AZTEC_LINOP;
    aztec_linop->bindx_2d = bindx_2d;
    (*linop)->aux_info = (void *) aztec_linop;
  }
  else {
    exit(1);
  }
  return;
}

void dft_destroy_linop(DFT_LINOP ** linop){

  if ((*linop)->vecspace->vecspace_type==DFT_AZTEC_VECSPACE) {

  }
  
  delete (DFT_LINOP *) (*linop);
  *linop = 0;
  return;
}
int dft_insert_global_graph_row(DFT_LINOP * linop, int i_mat, int * row_indices,
				       int nonzeros_in_row) {
  int N_loc = linop->vecspace->N_loc;
  if (linop->vecspace->vecspace_type==DFT_AZTEC_VECSPACE) {
    DFT_AZTEC_LINOP * aztec_linop = linop->aux_info;
    int ** bindx_2d = aztec_linop->bindx_2d;
    bindx_2d[i_mat] = row_indices;
    bindx_2d[N_loc][i_mat] = nonzeros_in_row;
  }
  linop->nonzeros += nonzeros_in_row;

  return(0);
}

int dft_fill_complete(DFT_LINOP * linop){

  int N_loc = linop->vecspace->N_loc;

  if (linop->vecspace->vecspace_type==DFT_AZTEC_VECSPACE) {
    DFT_AZTEC_LINOP * aztec_linop = linop->aux_info;
    int ** bindx_2d = aztec_linop->bindx_2d;

   /* Allocate MSR matrix with known size */

   min_nonzeros = AZ_gmin_int(Aztec.nonzeros,Aztec.proc_config);
   max_nonzeros = AZ_gmax_int(Aztec.nonzeros,Aztec.proc_config);
   tot_nonzeros = AZ_gsum_int(Aztec.nonzeros,Aztec.proc_config);

   time_preproc = MPI_Wtime()-t1;

   if (Proc ==0 && Iwrite != NO_SCREEN) {
      printf("\n\tMSR Preproc found %d to %d (total %d) nonzeros\n",
                                min_nonzeros,max_nonzeros, tot_nonzeros);
      printf("\tMSR Preproc took %g seconds\n", time_preproc);
   }

   Aztec.bindx  = (int *) array_alloc(1, Aztec.nonzeros+1, sizeof(int));
   Aztec.val = (double *) array_alloc(1, Aztec.nonzeros+1, sizeof(double));
   if (Aztec.val==NULL) {
      printf("%s ERROR: Not enough memory to allocate MSR matrix\n",yo);
      exit(-1);
   }

   /* load diagnal entries of bindx */

   Aztec.bindx[0] = Aztec.N_update + 1;
   for (i=0; i< Aztec.N_update; i++){
     Aztec.bindx[i+1] = Aztec.bindx[i] + bindx_2d[Aztec.N_update][i];
     /*printf("i: %d  Aztec.bindx[i]: %d  bindx_2d[Aztec.N_update][i]: %d\n",
              i,Aztec.bindx[i],bindx_2d[Aztec.N_update][i]);*/
   }

   /* load off-diagnal entries of bindx */

   for (i=0; i<Aztec.N_update; i++)
     for (j=0; j<bindx_2d[Aztec.N_update][i]; j++)
       Aztec.bindx[Aztec.bindx[i] + j] = bindx_2d[i][j];

   AZ_transform(Aztec.proc_config, &(Aztec.external), Aztec.bindx,
                Aztec.val, Aztec.update, &(Aztec.update_index),
                &(Aztec.extern_index), &(Aztec.data_org), Aztec.N_update,
                NULL, NULL, NULL, NULL, AZ_MSR_MATRIX);

    for (int i=0; i<N_loc+; i++) if (bindx_2d[i]>0) delete [] bindx_2d[i];
    delete [] bindx_2d;
  }
  return(0);
}

int dft_insert_local_mat_row(DFT_LINOP * linop, int i_mat, double * values){
  int i;
  /* following needed for ideal gas */

     val[i_mat] = mat_row[i_mat];
     mat_row[i_mat] = 0.0;
     for (j=bindx[i_mat]; j<bindx[i_mat+1]; j++) {
       val[j] = mat_row[bindx[j]];
       mat_row[bindx[j]] = 0.0;
     }
  return(0);
}

void dft_matvec(DFT_LINOP * linop, double *x, double *y){
  return(0);
}
