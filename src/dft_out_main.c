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
             *output_file6="dft_gofr.dat", *output_flux= "dft_flux.dat",
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

  if (Proc == 0) {
    X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
    Vext_old = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
  }

  collect_x_old(x);
  collect_vext_old();

   if (Proc == 0 && Iwrite != MINIMAL) {
        if (Print_rho_type != PRINT_RHO_0) print_profile(output_file5);
        else                                print_profile(output_file4);
   }
   if (Lprint_gofr && Nlink==1) print_gofr(output_file6);

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
	if(first) print_cont_type(Loca.cont_type1,fp);
	else print_cont_variable(Loca.cont_type1,fp);
      }
      if (Loca.method == 3 || Loca.method == 4) {
	if(first) print_cont_type(Loca.cont_type2,fp);
	else print_cont_variable(Loca.cont_type2,fp);
      }
      if (Nruns > 1) {
	if(first) print_cont_type(0,fp);
	else print_cont_variable(0,fp);
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

   if (Lsteady_state && Proc==0) calc_flux(fp,output_flux,X_old);

   if (Proc==0) fprintf(fp,"  \n");

   if(first) first = FALSE;

   } /* end of loop over out_loop */
   
   if (Proc==0) {
      fclose(fp);
   }

   if (Proc==0 && Iwrite !=NO_SCREEN) printf("post processing took %g secs\n",MPI_Wtime()-t1);
}
/******************************************************************************/
/*print_cont_variable: Here print the value of the variable that
                     is changing in a given run */
void print_cont_variable(int cont_type,FILE *fp)
{
   int i,idim,icomp,iwall,iwall_type,nloop;
   double kappa,kappa_sq,rhosum;


   switch(cont_type){
      case CONT_MESH: 
         if (Print_mesh_switch == SWITCH_SURFACE_SEP && (Nwall == 1 ||
                Nwall==2)){ 
            iwall = 0;
            iwall_type = WallType[iwall];
            idim = Plane_new_nodes;

            if (Nwall == 1) {
              if (Type_bc[idim][0] == REFLECT)
                  fprintf(fp,"%11.8f   ", 2.0*(WallPos[idim][iwall] 
                                              + 0.5*Size_x[0] 
                                              - WallParam[iwall_type]));
           
               else if (Type_bc[idim][1] == REFLECT) 
                  fprintf(fp,"%11.8f   ", 2.0*(0.5*Size_x[0] 
                                              - WallPos[idim][iwall]
                                              -WallParam[iwall_type]));
               else
                for (idim=0; idim<Ndim; idim++)
                    fprintf(fp,"%11.8f   ",WallPos[idim][iwall]);
             }
             else if (Nwall == 2){
                  fprintf(fp,"%11.8f   ",
                    (fabs(WallPos[idim][1] - WallPos[idim][0]) - 2.0*WallParam[iwall_type]));
             }
         }
         else if (Nwall > 1) {
            for (iwall=0; iwall<Nwall; iwall++){
                for (idim=0; idim<Ndim; idim++)
                    fprintf(fp,"%11.8f   ",WallPos[idim][iwall]);
            } 
         }
         else {
            for (idim=0; idim<Ndim; idim++) 
                 fprintf(fp,"%11.8f   ",Size_x[idim]);
         }

         break;

      case CONT_TEMP: 
         if (Ipot_ff_c == 0) fprintf(fp,"%10.7f   ", Temp); 
	 else fprintf(fp,"%7.4f   ",Temp_elec);
         break;

      case CONT_CRFAC:
	 fprintf(fp,"%11.8f   ",Crfac);
         break;
      case CONT_SCALE_RHO:
         fprintf(fp,"%11.8f   ", Scale_fac); 
      case CONT_RHO_0:
      case CONT_RHO_ALL:
      case CONT_LOG_RHO_0:
      case CONT_LOG_RHO_ALL:
         rhosum=0.0;
         nloop=Ncomp;
         for (i=0; i<nloop; i++){
                 fprintf(fp,"%11.8f  ", Rho_b[i]); 
                 rhosum+=Rho_b[i];
         }
         for (i=0;i<nloop;i++) fprintf(fp,"%9.6f  ",Rho_b[i]/rhosum);
         if (Print_rho_switch == SWITCH_RELP && nloop == 1)
              fprintf(fp,"%11.8f   ", P_over_po); 
         else if (Print_rho_switch == SWITCH_ION && Ipot_ff_c == COULOMB) {
             kappa_sq = 0.0;
             for(icomp = 0; icomp<nloop; icomp++)
                kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*
                           Charge_f[icomp]*Charge_f[icomp];
             kappa = sqrt(kappa_sq);
             fprintf(fp,"%11.8f   ", kappa); 
         }
         else if (Print_rho_switch == SWITCH_MU)
           for (i=0; i<nloop; i++) fprintf(fp,"%11.8f   ", Betamu[i]); 
         break;

      case CONT_BETAMU_0:
         fprintf(fp,"%11.8f   ", Betamu[0]); 

      case CONT_SCALE_EPSW:
         fprintf(fp,"%11.8f   ", Scale_fac); 
      case CONT_EPSW_0:
      case CONT_EPSW_ALL:
         if (Mix_type==0) fprintf(fp,"%11.8f   ", Eps_w[0]); 
         else fprintf(fp,"%11.8f   ", Eps_ww[0][0]); 
        
         break;
 
      case CONT_SCALE_EPSWF:
         fprintf(fp,"%11.8f   ", Scale_fac); 
         for (i=0;i<Ncomp;i++) fprintf(fp,"%9.6f  ",Eps_wf[i][0]);
         break;
      case CONT_EPSWF00:
      case CONT_EPSWF_ALL_0:
         fprintf(fp,"%11.8f   ", Eps_wf[0][0]); 
	 /*        fprintf(fp,"%11.8f   ", Eps_wf[2][0]); */
         break;

      case CONT_SCALE_EPSFF:
         fprintf(fp,"%11.8f   ", Scale_fac); 
      case CONT_EPSFF_00:
      case CONT_EPSFF_ALL:
         fprintf(fp,"%11.8f   ", Eps_ff[0][0]); 
     /*    fprintf(fp,"%11.8f  ", Eps_ff[2][0]); */
      /*   fprintf(fp,"%11.8f  ", Eps_ff[2][2]); */
         break;

      case CONT_SCALE_CHG:
         fprintf(fp,"%11.8f   ", Scale_fac); 
         break;
   
      case CONT_SEMIPERM:
         fprintf(fp,"%11.8f   ", Vext_membrane[0][0]); break;

      case CONT_WALLPARAM:
         fprintf(fp,"%11.8f   ", WallParam[1]); break;

   }
   return;

}
/******************************************************************************/
/*print_cont_type: Here print the type of the variable that
                     is changing in a given run 
		     Note the logic here should mirror that in print_cont_variable() above
*/
void print_cont_type(int cont_type,FILE *fp)
{
  int i,idim,icomp,iwall,nloop;

   switch(cont_type){
      case CONT_MESH: 
	if (Print_mesh_switch == SWITCH_SURFACE_SEP && (Nwall==1 || Nwall==2)) {
	  idim = Plane_new_nodes;
	  if(Nwall==1 && Type_bc[idim][0] != REFLECT && Type_bc[idim][1] != REFLECT) {
	    for(idim=0; idim<Ndim; idim++)
                fprintf(fp,"WallPos[%d][0]  ", idim);
	  }
	  else fprintf(fp, "Surf_sep  ");
	}
	else if (Nwall > 1) {  
	   for (iwall=0; iwall<Nwall; iwall++) {
	     for(idim=0; idim<Ndim; idim++)
	       fprintf(fp,"WallPos[%d][%d]  ", idim,iwall);
	   }
	 }
	 else {
	   for(idim=0; idim<Ndim; idim++)
	     fprintf(fp,"Size[idim=%d]  ",idim);
	 }

         break;

      case CONT_TEMP: 
         if (Ipot_ff_c == 0) fprintf(fp,"TEMP  "); 
	 else fprintf(fp,"TEMP_ELEC  "); 
         break;

      case CONT_CRFAC:
	fprintf(fp, "Crfac  "); break;
      case CONT_SCALE_RHO:
         fprintf(fp,"Scale_fac_rho  "); break;
      case CONT_RHO_0:
      case CONT_RHO_ALL:
      case CONT_LOG_RHO_0:
      case CONT_LOG_RHO_ALL:
	nloop=Ncomp;
	for (i=0; i<nloop; i++) fprintf(fp,"Rho_b[%d]  ", i);
	for (i=0; i<nloop; i++) fprintf(fp, "Rho_b[%d]/Rho_sum  ", i);
         if (Print_rho_switch == SWITCH_RELP && Ncomp == 1)
              fprintf(fp,"P_over_Po  "); 
         else if (Print_rho_switch == SWITCH_ION && Ipot_ff_c == COULOMB) 
              for(i=0; i<nloop; i++) fprintf(fp,"KAPPA[%d]   ",i); 
         else if (Print_rho_switch == SWITCH_MU)
              for(i=0; i<nloop; i++) fprintf(fp,"CHEM_POT[%d]  ",i);
         break;
   
      case CONT_BETAMU_0:
         fprintf(fp,"Betamu[0]:  "); break;

      case CONT_SCALE_EPSW:
         fprintf(fp,"Scale_fac_epsw:  "); break;
      case CONT_EPSW_0:
      case CONT_EPSW_ALL:
         if (Mix_type==0) fprintf(fp,"Eps_w[0]  ");
         else              fprintf(fp,"Eps_ww[0][0]  "); 
         break;

      case CONT_SCALE_EPSWF:
         fprintf(fp,"Scale_fac_epswf  "); 
	 for(i=0; i<Ncomp; i++) fprintf(fp, "Eps_wf[%d][0]  ",i);
	 break;
      case CONT_EPSWF00:
      case CONT_EPSWF_ALL_0:
         fprintf(fp,"Eps_wf[0][0]:  "); break;
	 /*         fprintf(fp,"Eps_wf[2][0]:  "); break;*/

      case CONT_SCALE_EPSFF:
         fprintf(fp,"Scale_fac_epsff:  "); break;
      case CONT_EPSFF_00:
      case CONT_EPSFF_ALL:
         fprintf(fp,"Eps_ff[0][0]:  "); break;

      case CONT_SCALE_CHG:
         fprintf(fp,"Scale_fac_chg:  "); break;

   case CONT_SEMIPERM:
     fprintf(fp, "Vext_membrane[0][0]  "); break;

   case CONT_WALLPARAM:
     fprintf(fp, "WallParam[1]  "); break;

   }
   return;

}
