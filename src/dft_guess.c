/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 *  FILE: dft_guess.c
 *
 *  This file contains high level logic for setting up an initial 
 *  guess for the solution vector.
 *
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
#include <string.h>
 
void set_initial_guess (int iguess, double** xOwned)
{
 /*
  * Local variable declarations
  */

  char filename[20];
  double t1=0.0;
  double *x_new;

  int iunk,i;
  int i1,start_no_info;

  int idim,iwall,iwall_type;
  double sep,fac=0.0;

  
  /********************** BEGIN EXECUTION ************************************/

  if (Proc==0) { 
    /*  printf("\n%s: Setting an initial guess ... ",yo);*/
      t1 = MPI_Wtime();
  }

 /*
  * Allocate array, and point local x to the start of it
  */

  if (Proc==0 && Iwrite==VERBOSE){
    if (Restart==0 && Imain_loop==0) printf("generating guess from scratch\n");  
    else printf("setting up guess from existing files\n");
  } 

  if (Restart > 0 || Imain_loop > 0){
       start_no_info = FALSE;
       for (i=0;i<NEQ_TYPE;i++) Restart_field[i]=FALSE;

        x_new = (double *) array_alloc(1, Nnodes*Nunk_per_node, sizeof(double));

        if (Proc == 0) {  /* Proc 0 reads in the data file */

           if ( Imain_loop == 0){

              /* START FROM AN OLD FILE - otherwise all of the _old variables
                                    were set in collect_xold (dft_output.c) */

              if (Print_rho_type != PRINT_RHO_0 && Imain_loop>0 ) {
                    if (Imain_loop>0)   i1=Imain_loop-1;
                    else                i1=0;
                    
                    if (Lbinodal && iguess==BINODAL_FLAG)
                          sprintf(filename, "dft_dens2.%0d",i1);
                    else  sprintf(filename, "dft_dens.%0d",i1);
              }
              else{
                    if (Lbinodal && iguess==BINODAL_FLAG)
                         sprintf(filename,"dft_dens2.dat"); 
                    else sprintf(filename,"dft_dens.dat");
              }

              Nodes_old = find_length_of_file(filename);
              if (Restart==5){ /* going to read in 1D file, need to dimension arrays for 2D or 3D system */
                  Nodes_old=Nnodes;  
if (Proc==0 && Iwrite != NO_SCREEN) printf("set Nodes_old to be Nnodes...Nnodes=%d Nodes_old=%d\n",Nnodes,Nodes_old);
              }
              if (Lbinodal && iguess==BINODAL_FLAG){
                 X2_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
              }
              else
                 X_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
              read_in_a_file(iguess,filename); /* Get X_old */
           }
           else{
              for (iunk=0;iunk<Nunk_per_node;iunk++) Restart_field[Unk2Phys[iunk]]=TRUE;
           }

           if (Nodes_old != Nnodes) {

if (Proc==0 && Iwrite != NO_SCREEN) printf("Nodes_old=%d  Nnodes=%d\n",Nodes_old,Nnodes);

             idim = Plane_new_nodes;
             iwall = 0;
             iwall_type = WallType[iwall];


             /* fix up to allow for some mixing of a bulk solution with a previously
                converged solution either with input parameter or as some smart guess, but
                can't be dependent on 2 surfaces in 1D */
/*             fac=Guess_range[0];*/
             fac=1.0;

             if (!start_no_info) shift_the_profile(x_new,fac); 
           }
           else{ 
             if (iguess==Iguess1){
               for (iunk=0; iunk<Nunknowns; iunk++) x_new[iunk] = X_old[iunk];
             }
             else if (iguess==BINODAL_FLAG && Lbinodal)
               for (iunk=0; iunk<Nunknowns; iunk++) x_new[iunk] = X2_old[iunk];
           }
           if (iguess==Iguess1)
              safe_free((void *) &X_old);
           else if (iguess==BINODAL_FLAG && Lbinodal)
              safe_free((void *) &X2_old);
        }

        MPI_Bcast (Restart_field,NEQ_TYPE,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast (&start_no_info,1,MPI_INT,0,MPI_COMM_WORLD);

        if (!start_no_info){
            if (Restart_field[DENSITY]==FALSE) {
                printf("can't restart without density fields yet\n");
                exit(-1);
            }
            communicate_profile(x_new,xOwned);
            check_zero_densities(xOwned);
            if (Lsteady_state && (Restart==3 || Restart_field[DIFFUSION]==FALSE))   setup_chem_pot(xOwned);
            if ((Type_poly==NONE||Type_poly==WTC) &&
                  (Restart==3 || Restart_field[HSRHOBAR]==FALSE)) setup_rho_bar(xOwned);
            if (Ipot_ff_c == COULOMB && (Restart==3 || Restart_field[POISSON]==FALSE)){
                   printf("setting up electrostatic potential guess....\n");
                   setup_elec_pot(xOwned,iguess);
            }
            if ((Type_poly != NONE && Type_poly!=WTC) && Restart_field[CMS_FIELD]==FALSE) setup_polymer_field(xOwned,iguess);   
            if (Restart == 2) chop_profile(xOwned,iguess);
        /*    if (Restart ==6 && Nwall>0) setup_exp_density_with_profile(xOwned);*/
        }
        safe_free((void *) &x_new);


  } /* end of logic for restarting from a file or previous run */

  else start_no_info = TRUE;

  if (start_no_info) {  /* START FROM A SPECIFIED INITIAL GUESS */

 if ((Type_poly==NONE || Type_poly==WTC)) {
    switch(iguess){
      case CONST_RHO:    
            if (Type_poly==WTC){
                 setup_const_density(xOwned,Rho_seg_b,Nseg_tot,0);
            }
            else              setup_const_density(xOwned,Rho_b,Ncomp,0);
            break;

      case CONST_RHO_V:  
            setup_const_density(xOwned,Rho_coex,1,1);
            break;
      case CONST_RHO_L:  
            setup_const_density(xOwned,Rho_coex,1,0);
            break;

      case EXP_RHO:
            setup_exp_density(xOwned,Rho_b,Ncomp,0);
            break;
      case EXP_RHO_V:
            setup_exp_density(xOwned,Rho_coex,1,1);
            break;
      case EXP_RHO_L:
            setup_exp_density(xOwned,Rho_coex,1,0);
            break;

      case STEP_PROFILE:
            setup_stepped_profile(xOwned);

      case CHOP_RHO_L:
            setup_exp_density(xOwned,Rho_coex,1,1);
            chop_profile(xOwned,iguess);
            break;
      case CHOP_RHO_V:
            setup_exp_density(xOwned,Rho_coex,1,0);
            chop_profile(xOwned,iguess);
            break;
      case LINEAR:
            setup_linear_profile(xOwned);
    }  /* end of iguess switch */

    setup_rho_bar(xOwned);
    if (Lsteady_state)   setup_chem_pot(xOwned);
    if (Type_poly==WTC){ 
         setup_BondWTC(xOwned);
         setup_Xi_cavWTC(xOwned);
    }
 }
 else{
    /*setup_polymer(x);*/

   if (Type_poly == CMS_SCFT)
     setup_polymer_simple(xOwned,iguess);
   else
     setup_polymer_rho(xOwned,iguess);
     setup_polymer_field(xOwned,iguess);
     setup_polymer_G(xOwned);
    /* setup_polymer_G_2(xOwned); */

    check_zero_densities(xOwned);
 }

 if (Ipot_ff_c == COULOMB)  setup_elec_pot(xOwned,iguess);

   } /* end of cases for setting initial guess from scratch */

  if (Proc==0 && Iwrite==VERBOSE) printf("\n initial guess took %g secs\n",MPI_Wtime()-t1);
  return;
}
