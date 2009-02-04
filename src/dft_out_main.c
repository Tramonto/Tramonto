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
 *  FILE: dft_out_main.c
 *
 *  This file contains routines that top level logic controling output to files.
 *
 */

#include "dft_out_main.h"
/*************************************************************************/
void post_process (double **x,char *output_file3,int *niters,
                   double *time_save, int loop1, int binodal_flag)
{
 /*
  * Local variable declarations
  */



  char *yo = "post_process",*output_file4 = "dft_dens.dat",*output_file5=NULL,
             *output_file6="dft_gofr.dat",*output_flux= "dft_flux.dat",
             *output_file7="dft_dens_site.dat",*output_file8=NULL;
  char filename[20];
  int icomp,iunk;
  double t1,energy;
  double fac_area,fac_vol;
  int i,iwall,idim;
  int out_loop;
  FILE *fp=NULL;
  static int first=TRUE;

  if (Print_rho_type != PRINT_RHO_0){
     if (binodal_flag){
       sprintf(filename, "dft_dens2.%0d", loop1);
       /*sprintf(filename_site, "dft_dens2_site.%0d", loop1);*/
     }
     else{
       sprintf(filename, "dft_dens.%0d", loop1);
/*       sprintf(filename_site, "dft_dens_site.%0d", loop1);*/
     }
     output_file5 = filename;
/*     output_file8 = filename_site;*/
  }

  if (binodal_flag){
     output_file4 = "dft_dens2.dat";
/*     output_file7 = "dft_dens2_site.dat";*/
  }

  if (Proc==0 && Iwrite != NO_SCREEN) 
           printf("\n%s: Doing post_processing and output of results ...\n",yo);
  t1 = MPI_Wtime();
 /*
  * First exchange boundary information as necessary !! 
  */

  if (Proc==0){
      if (binodal_flag)  X2_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
      else                X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
      Vext_old = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
  }

  if (binodal_flag) collect_x_old(x,X2_old);
  else              collect_x_old(x,X_old);
  collect_vext_old();

   if (Proc == 0 && Iwrite != MINIMAL) {
        if (binodal_flag){
           if (Print_rho_type != PRINT_RHO_0) print_profile(output_file5,X2_old);
           else                               print_profile(output_file4,X2_old);
        }
        else{
           if (Print_rho_type != PRINT_RHO_0) print_profile(output_file5,X_old);
           else                               print_profile(output_file4,X_old);
        }
   }
   if (Proc==0 && Lprint_gofr && Nlink==1){
       if (binodal_flag) print_gofr(output_file6,X2_old);
       else print_gofr(output_file6,X_old);
   }

   if (Proc==0) safe_free((void *) &Vext_old);


   /* open dft_output.dat file */
   if (Proc ==0){
      if( (fp = fopen(output_file3,"a"))==NULL) {
	printf("Can't open file %s\n", output_file3);
	exit(1);
      }
   }

  /* begin extra loop here */
   if (first) out_loop = 2;
   else  out_loop = 1;

   for(i=0; i<out_loop; i++) {

   if (Proc == 0) {
      if (Loca.method != -1) {
	if(first) print_cont_type(Loca.cont_type1,fp,0);
	else print_cont_variable(Loca.cont_type1,fp,0);
      }
      if (Loca.method == 3 || Loca.method == 4) {
	if(first) print_cont_type(Loca.cont_type2,fp,1);
	else print_cont_variable(Loca.cont_type2,fp,1);
      }
      if (Nruns > 1) {
	if(first) print_cont_type(0,fp,-1);
	else print_cont_variable(0,fp,-1);
      }
   }

   if (Proc==0){
      if (first){
           fprintf(fp,"niters  time  ");
      }
      else fprintf(fp,"%d  %9.4f  ",*niters,*time_save);
   }

   if(!first) 
     setup_domain_multipliers();

   /* calculate multiplicative factors that result         !!!!!!!!!!!remove this chuck ASAP
      from the presence of reflective boundaries...
      use the area of the 0th wall for an area basis calculation. */

   iwall = 0;
   fac_area = 1.0;
   fac_vol = 1.0;
   for (idim = 0; idim<Ndim; idim++) {
       if (!(Type_bc[idim][0] == REFLECT && Type_bc[idim][1] == REFLECT) && Lcount_reflect){

       if (Type_bc[idim][0] == REFLECT){ 
         
          fac_vol *= 2.0;
          if (WallPos[idim][0] == -0.5*Size_x[idim]) fac_area *= 2.0;
       }
       else if (Type_bc[idim][1] == REFLECT) {

          fac_vol *= 2.0;
          if (WallPos[idim][0] == 0.5*Size_x[idim]) fac_area *= 2.0;
       }
       }
   }

   calc_adsorption(fp,x);
   
   if (Type_coul != NONE) calc_fluid_charge(fp,x); 

/*   if (Ipot_wf_n != LJ12_6_WALL &&  
         Ipot_wf_n != LJ_CHARGED_ATOMS && Ipot_wf_n != LJ_ATOMIC) */
   calc_force(fp,x,fac_area);   
                            /* haven't implemented V_dash 
                               for 12-6 integrated wall yet */

   energy=calc_free_energy(fp,x); 

   if (Type_interface==DIFFUSIVE_INTERFACE && Proc==0) calc_flux(fp,output_flux,X_old);

   if (Proc==0) fprintf(fp,"  \n");

   if(first) first = FALSE;

   } /* end of loop over out_loop */
   
   if (Proc==0) {
      fclose(fp);
   }

   if (Proc==0 && Iwrite !=NO_SCREEN) printf("post processing took %g secs\n",MPI_Wtime()-t1);
}
/******************************************************************************/
