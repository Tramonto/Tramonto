/*@HEADER
// ***********************************************************************
// 
//                Tramonto: Molecular Theories Modeling Code
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// Questions? Contact Laura J.D. Frink (ljfrink@sandia.gov)
// 
// ***********************************************************************
//@HEADER
*/
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
