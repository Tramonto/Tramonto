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
#ifndef DFT_PMACHINE_H
#define DFT_PMACHINE_H
#include <stdio.h>
#include <mpi.h>

/*! \file dft_pmachine.h
\brief Function prototypes required for Tramonto parallel machine functionality

The file must be included by each Tramonto source file that needs access
to parallel machine functionality.  This file also specifies the functions that a 
parallel machine library must implement in order to supply functionality
for Tramonto.

*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*! \def DFT_AZTEC_PMACHINE
  \brief Selector for Aztec implementation of Tramonto's parallel machine.
*/
#define DFT_AZTEC_PMACHINE           1

/*! \def DFT_EPETRA_PMACHINE
  \brief Selector for Epetra's implementation of Tramonto's parallel machine.
*/
#define DFT_EPETRA_PMACHINE          2

/*! \def DFT_LAPACK_PMACHINE
  \brief Selector for a serial implementation of Tramonto's parallel machine.
*/
#define DFT_SERIAL_PMACHINE          3

/*! \struct DFT_PMACHINE_STRUCT dft_pmachine.h

\brief Contains all information about a parallel machine object.

*/
struct DFT_PMACHINE_STRUCT {

  int machine_type;    /*!< Selector for machine type */
  void * aux_info;  /*!< Generic pointer that can be used by the implementing library for addional data */
  
};

/*! \typedef DFT_PMACHINE

\brief Convenient expression for struct DFT_PMACHINE_STRUCT.

*/
typedef struct DFT_PMACHINE_STRUCT DFT_PMACHINE;

/*! \fn void dft_create_pmachine(int machine_type, MPI_Comm * comm, DFT_PMACHINE ** machine)

\brief Create a parallel machine object.

\param comm (In) MPI communicator.
\param machine (Out) The address of the parallel machine.

\pre comm must be properly defined.
\pre op must be an uninitialized pointer.  If it is already initialized, there will be 
a memory leak.

\post The MPI comm pointer will be copied, but the associated  data WILL NOT be copied.  It is the
responsibility of the calling program to maintain the integrity of the associated MPI comm object.

*/

void dft_create_pmachine(int machine_type, MPI_Comm * comm, DFT_PMACHINE ** machine);

/*! \fn void dft_destroy_pmachine(DFT_PMACHINE ** machine)

\brief Destroy a linear operator and delete all associated memory.

\param machine (In/Out) The address of the machine object that will be destroyed.

\pre machine must have been created by a call to dft_create_pmachine.

\post All memory associated with this object will be deleted and the pointer will be set to NULL.  
\post The associated MPI comm object will NOT be deleted.  This is the responsibility of the calling
program.

*/

void dft_destroy_pmachine(DFT_PMACHINE ** machine);

/*! \fn void dft_gather_global_vec(DFT_VEC_SPACE * vec_space, double *loc_vec, 
                                   int *loc_index,int N_loc, double *global_vec)

\brief Collect a distributed vector on root processor.

\param machine (In) A pointer to a valid vector space object.
 \param loc_vec (In) vector containing the solution vector index by the 
        local node number on the current processor.
 \param local_index (In) vector containing global node numbers of each of 
        the local processor nodes.
 \param N_loc (In) Length of loc_vec and loc_index arrays.
 \param global_vec (Output) global vector which has been gathered from all procs.

*/

void dft_gather_global_vec(DFT_PMACHINE * machine, double *loc_vec, 
                           int *loc_index,int N_loc, double *global_vec);

/*! \fn void dft_gather_global_vec_int(DFT_PMACHINE * machine, int *loc_vec, 
                                       int *loc_index,int N_loc, int *global_vec)

\brief Collect a distributed vector on root processor.

\param machine (In) A pointer to a valid vector space object.
 \param loc_vec (In) vector containing the solution vector index by the 
        local node number on the current processor. 
 \param local_index (In) vector containing global node numbers of each of 
        the local processor nodes.
 \param N_loc (In) Length of loc_vec and loc_index arrays.
 \param global_vec (Output) global vector which has been gathered from all procs.

*/

void dft_gather_global_vec_int(DFT_PMACHINE * machine, int *loc_vec, 
                               int *loc_index,int N_loc, int *global_vec);

/*! \fn dft_gsum_double(DFT_PMACHINE * machine, double partial_sum)

\brief Global sum.

\param machine (In) A pointer to a valid parallel machine object.
\param partial_sum (In) Partial sum from local results.

\return Global sum of partial sum values.

*/

double dft_gsum_double(DFT_PMACHINE * machine, double partial_sum);

/*! \fn dft_gmax_double(DFT_PMACHINE * machine, double partial_max)

\brief Global max.

\param machine (In) A pointer to a valid parallel machine object.
\param partial_sum (In) Partial max from local results.

\return Global max of partial max values.

*/
double dft_gmax_double(DFT_PMACHINE * machine, double partial_max);


/*! \fn dft_gavg_double(DFT_PMACHINE * machine, double partial_max)

\brief Global average.

\param machine (In) A pointer to a valid parallel machine object.
\param in_value (In) Partial max from local results.

\return Global max of partial max values.

*/
  int dft_gavg_double(DFT_PMACHINE * machine, double in_value);

/*! \fn dft_gmax_int(DFT_PMACHINE * machine, int partial_max)

\brief Global max.

\param machine (In) A pointer to a valid parallel machine object.
\param partial_sum (In) Partial max from local results.

\return Global max of partial max values.

*/
  int dft_gmax_int(DFT_PMACHINE * machine, int partial_max);

/*! \fn dft_gmin_int(DFT_PMACHINE * machine, int partial_min)

\brief Global max.

\param machine (In) A pointer to a valid parallel machine object.
\param partial_sum (In) Partial min from local results.

\return Global min of partial min values.

*/
  int dft_gmin_int(DFT_PMACHINE * machine,  int partial_max);

/*! \fn dft_gsum_vec_int(DFT_PMACHINE * machine, int * partial_sum, int * sum, in length)

\brief Integer vector global sum.

\param machine (In) A pointer to a valid parallel machine object.
\param partial_sum (In) Partial min from local results.
\param sum (Out) Global Vector sum from partials.
\param length (In) Length of vectors.

*/
  int dft_gsum_vec_int(DFT_PMACHINE * machine,  int partial_sum, int * sum, int length);

/*! \fn void dft_broadcast(DFT_PMACHINE * machine, char * buffer, int size_of_buffer);

\brief Broadcast from root node.

\param machine (In) A pointer to a valid parallel machine object.
\param buffer (In) Array of characters to broadcast from node 0.
\param size_of_buffer (In) Length of buffer in characters.

*/

double dft_broadcast(DFT_PMACHINE * machine, char * buffer, int size_of_buffer);


/*! \fn void dft_broadcast(DFT_PMACHINE * machine, char * buffer, int size_of_buffer);

\brief Barrier forcing all processors to synchronize execution.

\param machine (In) A pointer to a valid parallel machine object.

*/

double dft_barrier(DFT_PMACHINE * machine);




#endif /* DFT_PMACHINE_H */
