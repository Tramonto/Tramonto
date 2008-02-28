/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_guess.c
 *
 *  This file contains high level logic for setting up an initial 
 *  guess for the solution vector.
 *
 */

#include "dft_guess.h"
 
void set_initial_guess (int iguess, double** xOwned)
{
  double t1=0.0;
  int i;
  int start_no_info;

  if (Proc==0) { 
    /*  printf("\n%s: Setting an initial guess ... ",yo);*/
      t1 = MPI_Wtime();
  }

  if (Proc==0 && Iwrite==VERBOSE){
      if (Restart==0 && Imain_loop==0) printf("generating guess from scratch\n");  
      else printf("setting up guess from existing files\n");
  } 

  /* Get some or all parts of unknown vector from previously generated restart files.
     Note that these files do not need to have all of the solution vector and can 
     be patched up in some cases. */

  for (i=0;i<NEQ_TYPE;i++) Restart_field[i]=FALSE;
  if (Restart > 0 || Imain_loop > 0){
       start_no_info = FALSE;
       guess_restart_from_files(start_no_info,iguess,xOwned);  
  } 
  else start_no_info = TRUE;

                          /* this loop sets up portions of the solution vector that were not treated using 
                             restart files.  This can be part or all of the solution vector. */
  for (i=0;i<NEQ_TYPE;i++){
     switch(i){
         case DENSITY:
           if (Phys2Nunk[DENSITY]>0 && start_no_info){
               if (Type_poly == CMS || Type_poly ==CMS_SCFT){
                   setup_polymer_rho(xOwned,iguess);
               }
               else{
                   setup_density(xOwned,iguess);
               }
           }
           else if (Phys2Nunk[DENSITY]>0 && Restart_field[DENSITY]==FALSE){
                printf("we don't have the ability to restart without a density field at this time\n");
                exit(-1);
           }
           break;
         case MF_EQ:
           if (Phys2Nunk[MF_EQ]>0 && (start_no_info || Restart_field[MF_EQ]==FALSE ||Restart==3)) setup_mf_attract(xOwned); break;
         case HSRHOBAR:
           if (Phys2Nunk[HSRHOBAR]>0 && (start_no_info || Restart_field[HSRHOBAR]==FALSE ||Restart==3)) setup_rho_bar(xOwned); break;
         case POISSON:
           if (Phys2Nunk[POISSON]>0 && (start_no_info || Restart_field[POISSON]==FALSE ||Restart==3)) setup_elec_pot(xOwned,iguess); break;
         case DIFFUSION: 
           if (Phys2Nunk[DIFFUSION]>0 && (start_no_info || Restart_field[DIFFUSION]==FALSE ||Restart==3)) setup_chem_pot(xOwned); break;
         case CAVWTC:
           if (Phys2Nunk[CAVWTC]>0 && (start_no_info || Restart_field[CAVWTC]==FALSE)) setup_Xi_cavWTC(xOwned); break;
         case BONDWTC:
           if (Phys2Nunk[BONDWTC]>0 && (start_no_info || Restart_field[BONDWTC]==FALSE)) setup_BondWTC(xOwned); break;

         case WJDC_FIELD:
           if (Phys2Nunk[WJDC_FIELD]>0 && (start_no_info || Restart_field[WJDC_FIELD]==FALSE)){
                 setup_polymer_field_wjdc(xOwned); } break;
         case CMS_FIELD:
           if (Phys2Nunk[CMS_FIELD]>0 && (start_no_info || Restart_field[CMS_FIELD]==FALSE)) {
                 setup_polymer_field(xOwned,iguess); } break;
         case G_CHAIN:
           if (Phys2Nunk[G_CHAIN]>0 && (start_no_info || Restart_field[G_CHAIN]==FALSE)){
                if (Type_poly==CMS) setup_polymer_G(xOwned); 
                else if (Type_poly==WJDC) setup_polymer_G_wjdc(xOwned);
           }
           break;
     case YW_DENS:
       break;
         default:
           printf("problem with switch in initial guess\n");
           exit(-1);
           break;
     }
  }
  if (Restart == 2) chop_profile(xOwned,iguess);  /* special case for treating wetting problems */

  check_zero_densities(xOwned);              

/* note need to put the setup_poymer_simple functionality inside the various setup_polymer_CMS routines....*/

  if (Proc==0 && Iwrite==VERBOSE) printf("\n initial guess took %g secs\n",MPI_Wtime()-t1);
  return;
}
