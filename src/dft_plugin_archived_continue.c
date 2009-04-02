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



/* dft_plugin_archived_continue.c:  This routine contains continuation strategies that
   have been implemented for specific purposes. They are archived for future reference.  */

#include <stdio.h>
#include "dft_plugin_archived_continue.h"
/*****************************************************************************/
double get_init_param_archived_plugin(int cont_type,int Loca_contID)
{
  int i,j; 
  double param,sum;

  switch(cont_type){

      case CONT_LOG_RHO_I:
         if (Type_poly==NONE) return log(Rho_b[Cont_ID[Loca_contID][0]]);
         else{
            sum=0.0;
            for (i=0;i<Nseg_tot;i++){
            if (SegAll_to_Poly[i]==Cont_ID[Loca_contID][0]) {
                      sum += Rho_seg_b[i];
            }
            }
            param=log(sum);
            return param;
         }
      break;
 
      case CONT_RHO_CONST_RHOTOT58: 
        return Rho_b[2];
        break;
      case CONT_RHO_CONST_XSOLV: 
        return Rho_b[2];
        break;
      case CONT_RHO_ALL:
        return Rho_b[0];
        break;
      case CONT_EPSW_ALL:
        if (Mix_type==0) return Eps_w[0];
        else             return Eps_ww[0][0];
        break;
      case CONT_EPSWF_ALL:
        return Eps_wf[0][0];
        break; 
      case CONT_EPSWF_SOME:
        return Eps_wf[1][0];
        break;
      case CONT_EPSFF_ALL:
        return Eps_ff[0][0];
        break;
      case CONT_CRFAC:
           return Crfac;
           break;

      default:
        printf("ERROR: Unknown Continuation parameter in archived plugin %d\n",cont_type);
        exit(-1); 
        break;
  }
  return 0.0;
}
/*****************************************************************************/
void assign_param_archived_plugin(int cont_type, int Loca_contID, double param)
/* Note: Post_processing assumes the cont_type flags are the same as those
   used in Tramonto's own continuation */
{
  int i,j,icomp,jcomp,iw,inode;
  double ratio,temp_save,scale_save,eps_wf_save[NCOMP_MAX][NWALL_MAX_TYPE],param_save;
  char     *output_TF,*output_file1, *output_file2;
  
  output_file1 = "dft_out.lis";
  switch(cont_type){
     case CONT_LOG_RHO_I:
           ratio=1./Rho_b[Cont_ID[Loca_contID][0]];
           Rho_b[Cont_ID[Loca_contID][0]]        = exp(param);
           ratio *= Rho_b[Cont_ID[Loca_contID][0]];
           if (Type_poly !=NONE) {
               for (i=0;i<Nseg_tot; i++){
                  if (SegAll_to_Poly[i]==Cont_ID[Loca_contID][0]){
                     Rho_seg_b[i]*=ratio;
                  }
               }
           }
           if (Type_poly == CMS) setup_polymer_cr();
           recalculate_stencils();
           break;


      case CONT_RHO_CONST_RHOTOT58:   
          /*continuation at constant rho_tot=0.58 for a system of an 8-2-8 lipid and solvent*/
             Rho_b[2] = param;
             Rho_b[0] = (16./18.)*(0.58-Rho_b[2]);
             Rho_b[1] = (2./18.)*(0.58-Rho_b[2]);

             if (Type_poly == CMS) setup_polymer_cr();
             recalculate_stencils();
             thermodynamics(output_file1);
             break;

      case CONT_RHO_CONST_XSOLV:
          /* continuation in density for a system of an 8-2-8 lipid and solvent system where the 
             number fraction of solvent is held constant. */
             ratio=1./Rho_b[2]; 
             Rho_b[2]= param;    
             ratio*=Rho_b[2];
             Rho_b[0]*=ratio;
             Rho_b[1]*=ratio; 
             if (Type_poly == CMS) setup_polymer_cr();
             recalculate_stencils();
             thermodynamics(output_file1);
             break;

      case CONT_RHO_ALL: 
          /* continuation of all densities simultaneously */
             ratio=1./Rho_b[0];
             Rho_b[0]=param;
             ratio*=Rho_b[0];
             for (i=1;i<Ncomp;i++) Rho_b[i]*=ratio;

             if (Type_poly != NONE){
                  for (i=0; i<Nseg_tot;i++){
                      Rho_seg_b[i]*=ratio; 
                  }
             } 
             if (Type_poly == CMS) setup_polymer_cr();
             recalculate_stencils();
             thermodynamics(output_file1);
             break; 
                 
      case CONT_EPSW_ALL: 
            /* vary wall-wall interactions for all surfaces simultaneously */
            if (Mix_type==0){
                 ratio = 1.0/Eps_w[0];
                 Eps_w[0]=param;
                 ratio*= Eps_w[0];
                 for (i=0;i<Nwall_type;i++) Eps_w[i] *= ratio;

                 for (i=0; i<Ncomp; i++){ 
                   for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                 }
                 pot_parameters(NULL);
                 for (i=0; i<Ncomp; i++){
                   for (iw=0; iw<Nwall; iw++){
                      ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                      scale_vext_epswf(ratio,i,iw); 
                   }
                 }
             }
             else {
                   ratio = 1.0/Eps_ww[0][0];
                   Eps_ww[0][0]=param;
                   ratio *= Eps_ww[0][0];
                   for (i=0;i<Nwall_type;i++) 
                       for (j=0;j<Nwall_type;j++) Eps_ww[i][j] *=ratio; 
             }
             break;

      case CONT_EPSWF_ALL: 
           /* vary wall-fluid interactions for all surface-fluid interactions simultaneously */
            ratio = 1./Eps_wf[0][0];
            Eps_wf[0][0]=param;
            ratio *= Eps_wf[0][0];
            for (iw=0;iw<Nwall_type;iw++){
               for (i=0; i<Ncomp-1; i++){ 
                  Eps_wf[i][iw]*=ratio;
                  scale_vext_epswf(ratio,i,iw);
               }
            }
            break;

      case CONT_EPSWF_SOME:
          /* vary wall fluid interactions for two species against two surface types simultaneously */
            ratio = param/Eps_wf[1][0];
            Eps_wf[1][0]  = param;
/*            Eps_wf[1][1] *= ratio;*/
            Eps_wf[2][0] *= ratio;
/*            Eps_wf[2][1] *= ratio;*/
            scale_vext_epswf(ratio,1,0); 
/*            scale_vext_epswf(ratio,1,1); */
            scale_vext_epswf(ratio,2,0);
/*            scale_vext_epswf(ratio,2,1); */
            break;


      case CONT_EPSFF_ALL: 
          /* vary all of the fluid-fluid interactions simultaneously */
            ratio = Eps_ff[0][0];
            Eps_ff[0][0]=param;
            ratio *= Eps_ff[0][0];

            if (Mix_type==0) {
                for (i=0; i<Ncomp; i++){ 
                   for (j=0; j<Ncomp; j++) Eps_ff[i][j]*=ratio;
                }

                    /* find new wf interactions and correct vext */
                pot_parameters("dft_out.lis"); 
                for (i=0; i<Ncomp; i++){
                  for (iw=0; iw<Nwall; iw++){
                    ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                    scale_vext_epswf(ratio,i,iw); 
                }

              }
            }
            else{  
                for (i=0; i<Ncomp; i++){ 
                   for (j=0; j<Ncomp; j++) Eps_ff[i][j]*=ratio;
                }
            }
            if (Type_hsdiam == BH_DIAM){ 
                 calc_HS_diams();
                 calc_InvR_params();
            }
            if (Type_poly==CMS && Type_poly==CMS_SCFT) setup_polymer_cr();
            recalculate_stencils();
            thermodynamics("dft_out.lis");
            break;

     case CONT_CRFAC:
         Crfac=param;
         setup_polymer_cr();
         recalculate_stencils();
         break;

      default:
        printf("ERROR_apt: Unknown Continuation parameter in archived plugins. %d\n",cont_type);
        exit(-1); break;
  }
}
/*****************************************************************************/
/*print_cont_type_archived_plugin: Here print the type of the variable that
                     is changing in a given run -- archived special cases    */       
void print_cont_type_archived_plugin(int cont_type,FILE *fp,int Loca_contID)
{
  int i,j,idim,icomp,iwall,nloop,jcomp;

   switch(cont_type){
      case CONT_LOG_RHO_I:
        fprintf(fp,"Rho_b[%d]  ",Cont_ID[Loca_contID][0]);
         /*  alternate print types for the density variable */
        /* for (i=0; i<nloop; i++) fprintf(fp, "Rho_b[%d]/Rho_sum  ", i);
         if (Print_rho_switch == SWITCH_RELP && Ncomp == 1)
              fprintf(fp,"P_over_Po  ");
         else if (Print_rho_switch == SWITCH_ION && Ipot_ff_c == COULOMB)
              for(i=0; i<nloop; i++) fprintf(fp,"KAPPA[%d]   ",i);
         else if (Print_rho_switch == SWITCH_MU)
              for(i=0; i<nloop; i++) fprintf(fp,"CHEM_POT[%d]  ",i);    */
        break;

      case CONT_RHO_CONST_RHOTOT58:
         fprintf(fp,"Rho_b[0]  Rho_b[1]  Rho_b[2]");
         break;

      case CONT_RHO_CONST_XSOLV:
         fprintf(fp,"Rho_b[0]  Rho_b[1]  Rho_b[2]");
         break;

      case CONT_RHO_ALL:
         if (Ncomp <=5){
           for (i=0;i<Ncomp;i++) fprintf(fp,"Rho_b[%d]  ",i);
         } 
         else{
           fprintf(fp,"Rho_b[0]  ");
         }
 
      case CONT_EPSW_ALL:
         if (Mix_type==0){
            if (Nwall_type <=5){
              for (i=0;i<Nwall_type;i++) fprintf(fp,"Eps_w[%d]  ",i);
            } 
            else{
              fprintf(fp,"Eps_w[0](all)  ");
            }
          }
          else{ 
             if (Nwall_type<=2){
                for (i=0;i<Nwall_type;i++) 
                for (j=0;j<Nwall_type;j++) fprintf(fp,"Eps_ww[%d][%d]  ",i,j);
             }
             else fprintf(fp,"Eps_ww[0][0](all)  ");
          }
          break; 

      case CONT_EPSWF_ALL:
          if (Nwall_type*Ncomp<=5){
                for (i=0;i<Ncomp;i++) 
                for (j=0;j<Nwall_type;j++) fprintf(fp,"Eps_wf[%d][%d]  ",i,j);
          }
          else fprintf(fp,"Eps_wf[0][0](all)  ");
          break; 
    
      case CONT_EPSWF_SOME:
          /*fprintf(fp,"Eps_wf[1][0]  Eps_wf[1][1]  Eps_wf[2][0]  Eps_wf[2][1]  ");*/
          fprintf(fp,"Eps_wf[1][0]  Eps_wf[2][0] ");
          break;

      case CONT_EPSFF_ALL:
         if (Ncomp <=2){
              for (i=0;i<Ncomp;i++) 
                 for (j=0;j<Ncomp;j++) fprintf(fp,"Eps_ff[%d][%d]  ",i,j);
         }
         else fprintf(fp,"Eps_ff[0][0](all)  ");

      case CONT_CRFAC:
        fprintf(fp, "Crfac  ");
        break;

      default:
        printf("ERROR_apt: Unknown Continuation parameter in archived plugins. %d\n",cont_type);
        exit(-1); break;
   }
   return;
}
/******************************************************************************/
/*print_cont_variable_archived_plugin: Here print the value of the variable that
                     is changing in a given run  -- archived special cases */
void print_cont_variable_archived_plugin(int cont_type,FILE *fp,int Loca_contID)
{
  int i,j,idim,icomp,iwall,nloop,jcomp;

   switch(cont_type){
      case CONT_LOG_RHO_I:
                  fprintf(fp,"%11.8f  ",Rho_b[Cont_ID[Loca_contID][0]]);
         /* alternate ways to print density */
         /*
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
         */
         break;

      case CONT_RHO_CONST_RHOTOT58:
         fprintf(fp,"%11.7f  %11.7f  %11.7f  ",Rho_b[0],Rho_b[1],Rho_b[2]);
         break;

      case CONT_RHO_CONST_XSOLV:
         fprintf(fp,"%11.7f  %11.7f  %11.7f  ",Rho_b[0],Rho_b[1],Rho_b[2]);
         break;

      case CONT_RHO_ALL:
         if (Ncomp <=5){
           for (i=0;i<Ncomp;i++) fprintf(fp,"%11.7f  ",Rho_b[i]);
         } 
         else{
           fprintf(fp,"%11.7f  ",Rho_b[0]);
         }
 
      case CONT_EPSW_ALL:
         if (Mix_type==0){
            if (Nwall_type <=5){
              for (i=0;i<Nwall_type;i++) fprintf(fp,"%11.7f  ",Eps_w[i]);
            } 
            else{
              fprintf(fp,"%11.7f  ",Eps_w[0]);
            }
          }
          else{ 
             if (Nwall_type<=2){
                for (i=0;i<Nwall_type;i++) 
                for (j=0;j<Nwall_type;j++) fprintf(fp,"%11.7f  ",Eps_ww[i][j]);
             }
             else fprintf(fp,"%11.7f  ",Eps_ww[0][0]);
          }
          break; 

      case CONT_EPSWF_ALL:
          if (Nwall_type*Ncomp<=5){
                for (i=0;i<Ncomp;i++) 
                for (j=0;j<Nwall_type;j++) fprintf(fp,"%11.7f  ",Eps_wf[i][j]);
          }
          else fprintf(fp,"%11.7f  ",Eps_wf[0][0]);
          break; 
    
      case CONT_EPSWF_SOME:
          /*fprintf(fp,"%11.7f  %11.7f  %11.7f  %11.7f  ",Eps_wf[1][0],Eps_wf[1][1],Eps_wf[2][0],Eps_wf[2][1]);*/
          fprintf(fp,"%11.7f  %11.7f  ",Eps_wf[1][0],Eps_wf[2][0]);
          break;

      case CONT_EPSFF_ALL:
         if (Ncomp <=2){
              for (i=0;i<Ncomp;i++) 
                 for (j=0;j<Ncomp;j++) fprintf(fp,"%11.7f  ",Eps_ff[i][j]);
         }
         else fprintf(fp,"%11.f7  ",Eps_ff[0][0]);
  
      case CONT_CRFAC:
         fprintf(fp,"%11.8f   ",Crfac);
         break;

      default:
        printf("ERROR_apt: Unknown Continuation parameter in archived plugins. %d\n",cont_type);
        exit(-1); break;
   }
   return;
}
/******************************************************************************/
