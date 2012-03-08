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
void post_process (double **x,int *niters,
                   double *time_save, int loop1, int binodal_flag, int call_from_flag)
{
 /*
  * Local variable declarations
  */



  char *ForceSR_Filename="dft_force_sum_rule.dat",
       *AdsorptionSR_Filename="dft_ads_sum_rule.dat",
       *MainOutput_Filename="dft_output.dat",
       *GofR_Filename="dft_gofr.dat",
       *output_flux= "dft_flux.dat",
       *outPath=NULL;
  char filename[100],outPath_array[100],Density_file[100],DensityCounter_file[100];
 
  double t1,energy;
  double fac_area,fac_vol;
  int i,iwall,idim,first_local,counter,icomp,pol_num;
  int out_loop,contID0_tmp[2];
  FILE *fp=NULL;
  FILE *fpFSR=NULL;
  FILE *fpASR=NULL;
  static int first=TRUE;
  static int second=TRUE;
  static double omega_s_previous,omega_s_2previous,separation_previous,separation_2previous,derivative_previous;
  static double mu_previous,mu_2previous;
  double derivative,cont_var,surface_sep,derivative_avg,mu;

  if (Proc==0) strcpy(outPath_array,OutputFileDir);

  if (Nruns>1) loop1=Imain_loop;
  counter=loop1;

                        /* set up filename arrays with proper output directory */
  if (Proc==0){
     if (Print_rho_type != PRINT_RHO_0){
        if (binodal_flag==TRUE) sprintf(filename, "dft_dens2.%0d", counter); 
        else                    sprintf(filename, "dft_dens.%0d", counter); 

        strcpy(DensityCounter_file,strcat(strcat(outPath_array,"/"),filename));
        strcpy(outPath_array,OutputFileDir);
     }
     else {
        strcpy(Density_file,strcat(strcat(outPath_array,"/"),"dft_dens.dat"));
        strcpy(outPath_array,OutputFileDir);

        if (binodal_flag==TRUE){
           strcpy(Density_file,strcat(strcat(outPath_array,"/"),"dft_dens2.dat"));
           strcpy(outPath_array,OutputFileDir);
        }
     }
  }

  if (Proc==0 && Iwrite_screen != SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY) {
          printf("\n-------------------------------------------------------------------------------\n");
          printf("\nDoing post_processing and output of results ...\n");
  }
  t1 = MPI_Wtime();
 /*
  * First exchange boundary information as necessary !! 
  */

  if (Proc==0){
      if (binodal_flag==TRUE){
          if (X2_old==NULL) X2_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
       }
      else{
          if (X_old==NULL) X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
      }
       if (Vext_old==NULL) Vext_old = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
  }

  if (binodal_flag==TRUE)  collect_x_old(x,X2_old);
  else  collect_x_old(x,X_old);
  collect_vext_old();


                                /* always write density files */
   if (Proc == 0) {
        if (binodal_flag==TRUE){
           if (Print_rho_type != PRINT_RHO_0) print_profile((char *)DensityCounter_file,X2_old);
           else                               print_profile((char *)Density_file,X2_old);
        }
        else{
           if (Print_rho_type != PRINT_RHO_0) print_profile((char *)DensityCounter_file,X_old);
           else                               print_profile((char *)Density_file,X_old);
        }
   }
   if (Proc==0 && Lprint_gofr && (Nlink==1 || Nlocal_charge>0)){
       if (binodal_flag==TRUE) print_gofr(GofR_Filename,X2_old);
       else print_gofr(GofR_Filename,X_old);
   }

   if (Proc==0) safe_free((void *) &Vext_old);


   /* open dft_output.dat file */
   if (!(Nruns>1 && Loca.method!=-1 && call_from_flag==FROM_MAIN)){
   if (Proc ==0){
      if( (fp = fopen(MainOutput_Filename,"a"))==NULL) {
	printf("Can't open file %s\n", MainOutput_Filename);
      }
      if(Nruns>2 && (Nwall==1 || Nwall==2)){
      if( (fpFSR = fopen(ForceSR_Filename,"a"))==NULL) {
	printf("Can't open file %s\n", ForceSR_Filename);
	exit(1);
      }}
      if (Loca.method!=-1 && Loca.cont_type1==CONT_BETAMU_I){
      if( (fpASR = fopen(AdsorptionSR_Filename,"a"))==NULL) {
	printf("Can't open file %s\n", AdsorptionSR_Filename);
	exit(1);
      } }
   }
   if (first==TRUE) first_local=TRUE;
   else first_local=FALSE;

  /* begin extra loop here */
   if (first==TRUE) out_loop = 2;
   else  out_loop = 1;

   for(i=0; i<out_loop; i++) {

   if (Proc == 0) {
 
              /* first print state variables unless specifically told not to do so */
      if (Print_rho_switch != SWITCH_NO_STATEOUT){
        if ((Loca.method == -1 || 
           (Loca.cont_type1 !=CONT_RHO_I && Loca.cont_type1 !=CONT_BETAMU_I && Loca.cont_type1 != CONT_BETAMU_I_NEW)) &&
           ( (Loca.method !=3 && Loca.method !=4) ||
           (Loca.cont_type2 !=CONT_RHO_I && Loca.cont_type2 !=CONT_BETAMU_I && Loca.cont_type2 != CONT_BETAMU_I_NEW))) {

          for (icomp=0;icomp<2;icomp++) { contID0_tmp[icomp]=Cont_ID[0][icomp]; }
          Cont_ID[0][0]=0;

          if(first==TRUE){
                print_cont_type(CONT_RHO_I,fp,0);
                fprintf(fp," | ");
           }
           else cont_var=print_cont_variable(CONT_RHO_I,fp,0);
           for (icomp=0;icomp<2;icomp++) {
              Cont_ID[0][icomp]=contID0_tmp[icomp];
           }
         }
      }

      if (Loca.method != -1) {
	if(first==TRUE){
            print_cont_type(Loca.cont_type1,fp,0);
            fprintf(fp," | ");
        }
	else cont_var=print_cont_variable(Loca.cont_type1,fp,0);
      }
      if (Loca.method == 3 || Loca.method == 4) {
	if(first==TRUE){
             print_cont_type(Loca.cont_type2,fp,1); 
             fprintf(fp," | ");
        }
	else cont_var=print_cont_variable(Loca.cont_type2,fp,1); 
      }
      if (Nruns > 1) {
	if(first==TRUE){
             print_cont_type(CONT_MESH,fp,-1);
             fprintf(fp," | ");
        }
	else  cont_var=print_cont_variable(CONT_MESH,fp,-1);
      }
 
      if (Lprint_scaleFacWJDC){
         if (first==TRUE){
            for (icomp=0;icomp<Ncomp;icomp++) fprintf(fp,"Scale_fac_WJDC[%d]  ",icomp);
            fprintf(fp," | ");
         }
         else{
            for (pol_num=0; pol_num<Npol_comp;pol_num++) {
               for (icomp=0;icomp<Ncomp;icomp++) if (Nseg_type_pol[pol_num][icomp] > 0) fprintf(fp,"%11.8f ",Scale_fac_WJDC[pol_num][icomp]);
            }
         }
      }
   }

   if (Proc==0){
      if (first==TRUE) fprintf(fp,"niters  time  ");
      else fprintf(fp,"%d  %9.4f  ",*niters,*time_save);
   }

   if (first!=TRUE) setup_domain_multipliers();

   /* calculate multiplicative factors that result      
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
   if (Type_interface != DIFFUSIVE_INTERFACE) calc_force(fp,x,fac_area);   
   energy=calc_free_energy(fp,x); 

   if (((Nruns>2 && (Nwall==1 || Nwall==2)) || (Loca.method!=-1 && Loca.cont_type1==CONT_BETAMU_I)) && fabs(cont_var)>1.e-6){
      if (Loca.method!=-1 && Loca.cont_type1==CONT_BETAMU_I){ 
        mu=cont_var;
        printf("cont_var=%g  Betamu=%g first=%d\n",cont_var,Betamu[0],first);
      }
      else                                surface_sep=cont_var;
      if (first_local){
             /* store the results of the first data point */
            omega_s_previous=energy;
            if (Loca.method!=-1 && Loca.cont_type1==CONT_BETAMU_I) {
                mu_previous=mu;
                if (Proc==0) fprintf(fpASR,"betamu \t numerical derivative d(Omega_s)/dmu \n");
                printf("should have printed header in the file I think....\n");
            }
            else{ 
                separation_previous=surface_sep;
                if (Proc==0) fprintf(fpFSR,"separation \t numerical derivative d(Omega_s)/dH \n");
            }
      }
      else{
         if (first!=TRUE && second==TRUE){
            if (Loca.method!=-1 && Loca.cont_type1==CONT_BETAMU_I) {
               derivative=-(energy-omega_s_previous)/(mu-mu_previous);
               if (Proc==0) fprintf(fpASR,"%11.8f \t %11.8f\n",mu_previous,derivative);
            }
            else{
               derivative=-(energy-omega_s_previous)/(surface_sep-separation_previous);
               if (Proc==0) fprintf(fpFSR,"%11.8f \t %11.8f\n",separation_previous,derivative);
            }
            second=FALSE;
         }
         else{
            if (Loca.method!=-1 && Loca.cont_type1==CONT_BETAMU_I) {
               derivative=-(energy-omega_s_previous)/(mu-mu_previous);
               derivative_avg=0.5*(derivative+derivative_previous);
               if (Proc==0) fprintf(fpASR,"%11.8f \t %11.8f\n",mu_previous,derivative);
            }
            else{
               derivative=-(energy-omega_s_previous)/(surface_sep-separation_previous);
               derivative_avg=0.5*(derivative+derivative_previous);
               if (Proc==0) fprintf(fpFSR,"%11.8f \t %11.8f\n",separation_previous,derivative_avg);
            }
         }
         if (Loca.method!=-1 && Loca.cont_type1==CONT_BETAMU_I) {
           mu_2previous=mu_previous;
           mu_previous=mu;
         }
         else{
            separation_2previous=separation_previous;
            separation_previous=surface_sep;
         }
         omega_s_2previous=omega_s_previous;
         omega_s_previous=energy;
         derivative_previous=derivative;
       }
   }

   if (Type_interface==DIFFUSIVE_INTERFACE && Proc==0){
        calc_flux(fp,output_flux,X_old);
    }

   if (Proc==0) fprintf(fp,"  \n");

   if(first==TRUE) first = FALSE;

   } /* end of loop over out_loop */
   
   if (Proc==0)  fclose(fp); 
   if (Loca.cont_type1==CONT_BETAMU_I){ if (Proc==0)  fclose(fpASR); }
   if (Nruns>2 && (Nwall==1 || Nwall==2)){ if (Proc==0)  fclose(fpFSR); }
   }

  if (Proc==0 && Iwrite_screen != SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY) 
     printf("\n-------------------------------------------------------------------------------\n");

  return;
}
/******************************************************************************/
