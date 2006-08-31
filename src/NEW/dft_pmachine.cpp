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

#include "dft_pmachine.h"
#include "az_aztec.h"
#include <cstdlib>
void dft_create_pmachine(int machine_type, MPI_Comm * comm, DFT_PMACHINE ** machine) {

  int * proc_config;
 
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(*comm, &size);
  MPI_Comm_rank(*comm, &rank);
  *machine = (DFT_PMACHINE *) new DFT_PMACHINE;
  (*machine)->machine_type = machine_type;
  (*machine)->size = size;
  (*machine)->rank = rank;

  if (machine_type==DFT_AZTEC_PMACHINE) {
    proc_config = (int *) new int[AZ_PROC_SIZE];
    AZ_set_proc_config(proc_config,*comm);
    (*machine)->aux_info = (void *) proc_config;
  }
  else {
    exit(1);
  }
  return;
}

void dft_destroy_pmachine(DFT_PMACHINE ** machine) {
  
  if ((*machine)->machine_type==DFT_AZTEC_PMACHINE) {
    delete [] (int *) (*machine)->aux_info;
  }
  else {
    exit(1);
  }
  delete (DFT_PMACHINE *) (*machine);
  *machine = 0;
}


void dft_gather_global_vec(DFT_PMACHINE * machine, double *loc_vec, 
			   int *loc_index, int N_loc, double *global_vec){return;}

void dft_gather_global_vec_int(DFT_PMACHINE * machine, int *loc_vec, 
			       int *loc_index,int N_loc, int *global_vec){return;}

double dft_gsum_double(DFT_PMACHINE * machine, double partial_sum){return(0.0);}

double dft_gmax_double(DFT_PMACHINE * machine, double partial_max){return(0.0);}

int dft_gavg_double(DFT_PMACHINE * machine, double in_value){return(0);}

int dft_gmax_int(DFT_PMACHINE * machine, int partial_max){return(0);}

int dft_gmin_int(DFT_PMACHINE * machine,  int partial_max){return(0);}

int dft_gsum_vec_int(DFT_PMACHINE * machine,  int partial_sum, int * sum, int length){return(0);}

void dft_broadcast(DFT_PMACHINE * machine, char * buffer, int size_of_buffer){return;}

void dft_barrier(DFT_PMACHINE * machine){}
