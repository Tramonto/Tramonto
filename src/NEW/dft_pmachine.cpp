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

struct DFT_AZTEC_COMM_STRUCT {
  int * proc_config;
};
typedef struct DFT_AZTEC_COMM_STRUCT DFT_AZTEC_COMM;

void dft_create_pmachine(int machine_type, MPI_Comm * comm, DFT_PMACHINE ** machine) {

  int * proc_config;
 
  *machine = (DFT_PMACHINE *) malloc(sizeof(DFT_PMACHINE));
  (*machine)->machine_type = machine_type;
  if (machine_type==DFT_AZTEC_PMACHINE) {
    proc_config = (int *) malloc(AZ_PROC_SIZE*sizeof(int));
    AZ_set_proc_config(proc_config,*comm);
    (*machine)->aux_info = (void *) proc_config;
  }
  else {
    abort();
  }
  return;
}

void dft_destroy_pmachine(DFT_PMACHINE ** machine) {
  
  if ((*machine)->machine_type==DFT_AZTEC_PMACHINE) {
    free((char*)(*machine)->aux_info);
  }
  else {
    abort();
  }
}


void dft_gather_global_vec(DFT_PMACHINE * machine, double *loc_vec, 
			   int *loc_index, int N_loc, double *global_vec){}

void dft_gather_global_vec_int(DFT_PMACHINE * machine, int *loc_vec, 
			       int *loc_index,int N_loc, int *global_vec){}

double dft_gsum_double(DFT_PMACHINE * machine, double partial_sum){}

double dft_gmax_double(DFT_PMACHINE * machine, double partial_max){}

int dft_gavg_double(DFT_PMACHINE * machine, double in_value){}

int dft_gmax_int(DFT_PMACHINE * machine, int partial_max){}

int dft_gmin_int(DFT_PMACHINE * machine,  int partial_max){}

int dft_gsum_vec_int(DFT_PMACHINE * machine,  int partial_sum, int * sum, int length){}

double dft_broadcast(DFT_PMACHINE * machine, char * buffer, int size_of_buffer){}

double dft_barrier(DFT_PMACHINE * machine){}
