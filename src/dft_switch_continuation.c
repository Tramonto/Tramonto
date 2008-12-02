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



/* dft_switch_continuation.c:  This routine sets up parameters to be used in 
   continuation calculations using the LOCA library.   */

#include "dft_switch_continuation.h"
/*****************************************************************************/
double get_init_param_value(int cont_type)
{
  int i,j; 
  double param,sum;

  switch(cont_type){
      case CONT_MESH: 
       printf("ERROR: Continuation Library cannot do mesh size changes\n");
       exit(-1); break;

      case CONT_TEMP: return Temp; break;

      case CONT_RHO_0:   return Rho_b[0]; break;
      case CONT_RHO_ALL:  return Rho_b[0]; break; 
             if (Type_poly ==NONE){
                for (i=0;i<Ncomp;i++) if (Rho_b[i] != Rho_b[0]){
                   printf("ERROR: need all Rho_b to be the same for CONT_RHO_ALL\n"); 
                   exit(-1);
                }
                param=Rho_b[0];
             }
             else if (Npol_comp ==1){ /* assume we are continuing in one molecular density */
                /*sum=0.;
                for (i=0;i<Nseg_tot;i++) sum += Rho_seg_b[i];
                param=sum;*/
                param=Rho_seg_b[0];  /*segment densities are identical = Rho_b/Nmer */
             }
             else{
                printf("ERROR: continue either with identical initial densities or with a single molecule\n");
                exit(-1);
             }
             return param; break;

      case CONT_LOG_RHO_0: return log(Rho_b[0]); break;

      case CONT_LOG_RHO_ALL:   
             for (i=0;i<Ncomp;i++) if (Rho_b[i] != Rho_b[0]){
                 printf("ERROR: need all Rho_b to be the same for CONT_LOG_RHO_ALL\n"); 
                 exit(-1);
             }
             return log(Rho_b[0]); break;

      case CONT_SCALE_RHO: return Scale_fac; break;

      case CONT_BETAMU_0: 
           if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) return Betamu_chain[0];
           else                 return Betamu[0];
           break;

      case CONT_BETAMU_1: 
           if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) return Betamu_chain[1];
           else                 return Betamu[1];
           break;

      case CONT_EPSW_0:   if (Mix_type==0) return Eps_w[0];
                          else             return Eps_ww[0][0];
                          break;

      case CONT_EPSW_ALL: 
             if (Mix_type ==0 ){
                for (i=0;i<Nwall_type;i++) if (Eps_w[i] != Eps_w[0]){
                    printf("ERROR: need all Eps_w to be the same for CONT_EPS_W_ALLL\n"); 
                    exit(-1);
                }
                return Eps_w[0];
             }
             else{
                for (i=0;i<Nwall_type;i++){
                   for (j=0; j<Nwall_type;j++)  if (Eps_ww[i][j] != Eps_ww[0][0]){
                      printf("ERROR: need all Eps_ww to be the same for CONT_EPS_W_ALLL Mix_type=1\n"); 
                   }
                }
                return Eps_ww[0][0];
             }
             break;
      case CONT_SCALE_EPSW: return Scale_fac; break;

      case CONT_EPSWF00: return Eps_wf[0][0]; /*return Eps_wf[2][0];*/

      case CONT_EPSWF_ALL_0: 
/*           for (i=0;i<Ncomp-1;i++) if (Eps_wf[i][0] != Eps_wf[0][0]) {
                 printf("ERROR: all Eps_wf must be equal for CONT_EPSWF_ALL_0\n");
                 exit(-1);
           }
           return Eps_wf[0][0]; break;*/
           return Eps_wf[2][0]; /*return Eps_wf[2][0];*/
      case CONT_SCALE_EPSWF: return Scale_fac; break;

      case CONT_EPSFF_00:   
           return Eps_ff[0][2]; 
           break;
      case CONT_EPSFF_ALL:
/*           for (i=0; i<Ncomp; i++) if (Eps_ff[i][i] != Eps_ff[0][0]) {
                 printf("ERROR: need all Eps_ff[i][i] to be equal for  (CONT_EPSFF_ALL)\n");
                 exit(-1);
           }
           return Eps_ff[0][0]; break;*/
           return Eps_ff[0][1]; break;
      case CONT_SCALE_EPSFF: return Scale_fac; break;

      case CONT_SCALE_CHG:  return Scale_fac; break;

      case CONT_SEMIPERM: return Vext_membrane[0][0]; break;

      case CONT_WALLPARAM: return WallParam[WallType[1]]; break;

      case CONT_CRFAC:
              return Crfac; break;

      default:
        printf("ERROR: Unknown Continuation parameter %d\n",cont_type);
        exit(-1); break;
  }
  return 0.0;
}
/*****************************************************************************/
void assign_parameter_tramonto(int cont_type, double param)
/* Note: Post_processing assumes the cont_type flags are the same as those
   used in Tramonto's own continuation */
{
  int i,j,icomp,jcomp,iw,inode;
  double ratio,temp_save,scale_save,eps_wf_save[NCOMP_MAX][NWALL_MAX_TYPE],param_save;
  char     *output_TF,*output_file1, *output_file2;
  
  output_file1 = "dft_out.lis";
  switch(cont_type){
     case CONT_MESH: 
       printf("ERROR: Continuation Library cannot do mesh size changes\n");
       exit(-1); break;

      case CONT_TEMP: 
                      temp_save = Temp;
                      Temp      = param;
                      ratio = Temp/temp_save;
                      if (Ipot_ff_c == COULOMB ) Temp_elec *=ratio;
                      if (Mix_type==0){
                         if (Ipot_ff_n == LJ12_6){
                            for (icomp=0; icomp<Ncomp; icomp++) Eps_ff[icomp][icomp] /= ratio;
                         }
                         for (i=0; i<Nwall_type;i++){
                             if(Ipot_wf_n[i] != VEXT_HARD) Eps_w[i] /= ratio;
                         }
                         pot_parameters(NULL);
                      }
                      else if (Mix_type==1){
                           for (icomp=0; icomp<Ncomp; icomp++){
                             for(jcomp=0; jcomp<Ncomp; jcomp++) Eps_ff[icomp][jcomp] /= ratio;
                             for (i=0; i<Nwall_type;i++) Eps_wf[icomp][i] /= ratio;
                           }
                           for (i=0; i<Nwall_type;i++) {
                              for (j=0;j<Nwall_type;j++) Eps_ww[i][j] /= ratio;
                           }
                      }
                      if (Type_func != NONE){ calc_HS_diams(); 
                                              calc_InvR_params();}
                      if (Type_poly ==CMS) setup_polymer_cr();
                      recalculate_stencils();
                      if (Nwall>0) scale_vext_temp(ratio);
                      break;

      case CONT_RHO_0:   
                       Rho_b[0]=param;
/*                       Rho_b[1]=Rho_b[0]/8;*/
/*                      ratio=1./Rho_b[0];  vary rho tot at const x_s
                      Rho_b[2]= param;    */
/*                      ratio*=Rho_b[2];
                      Rho_b[0]*=ratio;
                      Rho_b[1]*=ratio;  vary rho tot at const x_s*/

/*                         Rho_b[0] = (16./18.)*(0.3771-Rho_b[2]);
                         Rho_b[1] = (2./18.)*(0.3771-Rho_b[2]);*/

                          /*continuation at constant rho_tot=0.58 */
                      /*   Rho_b[0] = (16./18.)*(0.58-Rho_b[2]);
                         Rho_b[1] = (2./18.)*(0.58-Rho_b[2]);*/

                          /*continuation at constant rho_tot=0.825 */
                         /*Rho_b[0] = (16./18.)*(0.825-Rho_b[2]);
                         Rho_b[1] = (2./18.)*(0.825-Rho_b[2]);*/

                         if (Type_poly == CMS) setup_polymer_cr();
                         recalculate_stencils();
                         break;
      case CONT_RHO_ALL: 
               /*Rho_b[2]=param;
		  Rho_b[0]=param;
		  Rho_b[1]=param/8;*/

                     /* ratio=1./Rho_b[2];*/  /*vary rho tot at const x_s*/
                    /*  Rho_b[2]= param;    
                      ratio*=Rho_b[2];
                      Rho_b[0]*=ratio;
                      Rho_b[1]*=ratio; */  /*vary rho tot at const x_s*/

             if (Type_poly==NONE)    for (i=0; i<Ncomp;i++)  Rho_b[i]= param;   
             else if (Npol_comp ==1){
                  for (i=0; i<Ncomp;i++)  Rho_b[i]= 0.;   
                  for (i=0; i<Nseg_tot;i++){
                      Rho_b[Unk2Comp[i]] += param;    
                      Rho_seg_b[i]=param; 
                  }
             } 
             if (Type_poly == CMS) setup_polymer_cr();
             recalculate_stencils();
             break; 
                 
      case CONT_LOG_RHO_0: Rho_b[0]        = exp(param);    break;
      case CONT_LOG_RHO_ALL: for (i=0;i<Ncomp;i++) Rho_b[i]= exp(param);    break;
      case CONT_SCALE_RHO: scale_save = Scale_fac;
                           Scale_fac = param;
                           ratio = Scale_fac/scale_save;
                           for (i=0; i<Ncomp-1; i++) Rho_b[i] *= ratio;
                           Rho_b[Ncomp-1]=0.678;
                           for (i=0; i<Ncomp-1; i++) Rho_b[Ncomp-1] -= Rho_b[i];
                           recalculate_stencils();
                           break;

      case CONT_BETAMU_0: if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)  Betamu_chain[0]=param;
                          else                 Betamu[0]=param;
                          break;

      case CONT_BETAMU_1: if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)  Betamu_chain[1]=param;
                          else                 Betamu[1]=param;
                          break;

      case CONT_EPSW_0:  
                      if (Mix_type ==0) {
                         Eps_w[0] = param;
                         for (i=0; i<Ncomp; i++){ 
                           for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                         }
                         pot_parameters(NULL);
                         for (i=0; i<Ncomp; i++){
                            ratio = Eps_wf[i][0]/eps_wf_save[i][WallType[0]];
                            scale_vext_epswf(ratio,i,0); 
                         }
                      }
                      else Eps_ww[0][0]=param;
                      break;

      case CONT_EPSW_ALL: 
                      if (Mix_type==0){
                          for (i=0;i<Nwall_type;i++) Eps_w[i] = param;
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
                            for (i=0;i<Nwall_type;i++) 
                                for (j=0;j<Nwall_type;j++) Eps_ww[i][j] = param;
                      }
                      break;

      case CONT_SCALE_EPSW: 
                      scale_save = Scale_fac;
                      Scale_fac = param;
                      ratio = Scale_fac/scale_save;
                      if (Mix_type==0){
                          for (iw=0; iw <Nwall_type;iw++) Eps_w[iw] *= ratio;
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
                        else{
                            for (i=0;i<Nwall_type;i++) 
                                for (j=0;j<Nwall_type;j++) Eps_ww[i][j] *=ratio;
                        }
                        break;

      case CONT_EPSWF00: 
                         ratio = param/Eps_wf[0][0];
                         Eps_wf[0][0]    = param;
                         scale_vext_epswf(ratio,0,0); break;

                          /* component 2 wall 0 */
/*                         ratio = param/Eps_wf[2][0];
                         Eps_wf[2][0]    = param;
                         scale_vext_epswf(ratio,2,0); break;*/

/*                         ratio = param/Eps_wf[1][0];
                         Eps_wf[1][0]    = param;
                         Eps_wf[1][1]    = param;
                         Eps_wf[2][0]    = param;
                         Eps_wf[2][1]    = param;
                         scale_vext_epswf(ratio,1,0); break;
                         scale_vext_epswf(ratio,1,1); break;
                         scale_vext_epswf(ratio,2,0); break;
                         scale_vext_epswf(ratio,2,1); break;*/

      case CONT_EPSWF_ALL_0: 

                         ratio = param/Eps_wf[2][0];
                         Eps_wf[2][0]    = param;
                         scale_vext_epswf(ratio,2,0); break;

/*                         for (i=0; i<Ncomp-1; i++){ 
                            for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                         }
                         for (i=0; i<Ncomp-1; i++){
                              Eps_wf[i][0] = param;
                              ratio = Eps_wf[i][0]/eps_wf_save[i][WallType[0]];
                              scale_vext_epswf(ratio,i,0);
                         }
                         break;*/
      
      case CONT_SCALE_EPSWF:
                          scale_save = Scale_fac;
                          Scale_fac = param;
                          ratio = Scale_fac/scale_save;
                          for (i=0; i<Ncomp; i++){ 
                            for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                          }
                          for (i=0; i<Ncomp; i++){
                             for (iw=0; iw<Nwall_type; iw++){
                                Eps_wf[i][iw] *= ratio;
                             }
                          }
                          for (i=0; i<Ncomp; i++){
                             for (iw=0; iw<Nwall; iw++){
                                ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                                scale_vext_epswf(ratio,i,iw);
                             }
                          }
                          break;

      case CONT_EPSFF_00: Eps_ff[0][0]=param;  

/*now do a special case where we change two of them at once */
/*                         Eps_ff[2][0]=param;
                         Eps_ff[0][2]=param;*/
                  /*       Eps_ff[1][2]=param;
                         Eps_ff[2][1]=param;*/

                      /*   Eps_ff[2][2]=param;
                         Eps_ff[1][2]=param;
                         Eps_ff[2][1]=param;*/

                       /*  Eps_ff[2][0]=param;*/
                  /*       if (Mix_type==0) {
                             for (iw=0; iw<Nwall_type; iw++) eps_wf_save[0][iw]=Eps_wf[0][iw];
                             pot_parameters("dft_out.lis"); 
                             for (iw=0; iw<Nwall; iw++){
                                 ratio = Eps_wf[0][WallType[iw]]/eps_wf_save[0][WallType[iw]];
                                 scale_vext_epswf(ratio,0,iw); 
                             }
                         }*/
                         if (Type_func != NONE){ calc_HS_diams(); 
                                                 calc_InvR_params();}
                         if (Type_poly == CMS) setup_polymer_cr();
                         recalculate_stencils();
                         break;

      case CONT_EPSFF_ALL: 
                       Eps_ff[0][1]=param;
                       Eps_ff[0][2]=param;
                       Eps_ff[1][0]=param;
                       Eps_ff[2][0]=param;
/*                         for (i=0; i<Ncomp; i++) 
                            for (j=0; j<Ncomp; j++) if (fabs(Eps_ff[i][j])>1.e-15) Eps_ff[i][j] = param;
                      
                         if (Mix_type==0) {
                             for (i=0; i<Ncomp; i++){ 
                              for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                             }
                             pot_parameters("dft_out.lis"); 
                             for (i=0; i<Ncomp; i++){
                               for (iw=0; iw<Nwall; iw++){
                             
                                ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                                scale_vext_epswf(ratio,i,iw); 
                             }
                           }
                         }
                         if (Type_poly==CMS && Type_poly==CMS_SCFT) setup_polymer_cr();*/
                         recalculate_stencils();
                         break;

      case CONT_SCALE_EPSFF:
                         scale_save = Scale_fac;
                         Scale_fac = param;
                         ratio = Scale_fac/scale_save;
                         for (i=0; i<Ncomp; i++) Eps_ff[i][i] *= ratio;
                         if (Mix_type==0) {
                             for (i=0; i<Ncomp; i++){ 
                              for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
                             }
                             pot_parameters("dft_out.lis"); 
                             for (i=0; i<Ncomp; i++){
                               for (iw=0; iw<Nwall; iw++){
                             
                                ratio = Eps_wf[i][WallType[iw]]/eps_wf_save[i][WallType[iw]];
                                scale_vext_epswf(ratio,i,iw); 
                             }
                           }
                         }
                         if (Type_poly==CMS) setup_polymer_cr();
                         recalculate_stencils();
                         break;

      case CONT_SCALE_CHG: scale_save=Scale_fac;
                           Scale_fac=param;
                           ratio = Scale_fac/scale_save;
                           scale_elec_param(ratio); 
                           break;

      case CONT_SEMIPERM:
           param_save=Vext_membrane[0][0];
           for (inode=0;inode<Nnodes_per_proc;inode++){
                  if (fabs(Vext[inode][0]-param_save)<1.e-10) Vext[inode][0]=param;
           }
           Vext_membrane[0][0]=param;
           break;

      case CONT_WALLPARAM:
           WallParam[1]=param;
           WallPos[0][1] = 0.5*Size_x[0]-WallParam[WallType[1]]; 
           WallPos[0][0] = 0.5*Size_x[0]-2.*WallParam[WallType[1]]-WallParam[WallType[0]]; 
           free_mesh_arrays();
           boundary_free();
           if (Type_interface==DIFFUSIVE_INTERFACE && Ndim==1) safe_free((void *) &Area_IC);
           safe_free((void *) &Comm_node_proc);
           safe_free((void *) &Comm_unk_proc);
           safe_free((void *) &Comm_offset_node);
           safe_free((void *) &Comm_offset_unk);
           output_file2 = "dft_vext.dat";
           output_TF = "dft_zeroTF.dat";
           set_up_mesh(output_file1,output_file2);
           boundary_setup(output_file1);
           if (Iwrite==VERBOSE) {
                 print_vext(Vext,output_file2);
                 print_zeroTF(Zero_density_TF,output_TF);
           }

           break;

     case CONT_CRFAC:
         Crfac=param;
         setup_polymer_cr();
         recalculate_stencils();
         break;

      default:
        printf("ERROR_apt: Unknown Continuation parameter %d\n",cont_type);
        exit(-1); break;
  }

  /* for most cases...recalculate thermo based on new parameter.  However if
     calculating bulk_coexistence or varying Betamu do not call thermo */
  if (Loca.cont_type1 != CONT_BETAMU_0 && Loca.cont_type1 != CONT_BETAMU_1 &&
     !(Loca.method==4 && Loca.cont_type2 == CONT_BETAMU_0) &&
     !(Loca.method==4 && Loca.cont_type2 == CONT_BETAMU_1)){
          thermodynamics(output_file1);
  }
}
/*****************************************************************************/
