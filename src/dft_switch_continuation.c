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

#include <stdio.h>
#include "dft_switch_continuation.h"
/*****************************************************************************/
double get_init_param_value(int cont_type,int Loca_contID)
{
  int i,j; 
  double param,sum;

  switch(cont_type){
      case CONT_MESH: 
       printf("ERROR: Continuation Library cannot do mesh size changes\n");
       exit(-1); break;

      case CONT_TEMP: return Temp; break;

      case CONT_RHO_I:   
             if (Type_poly==NONE) return Rho_b[Cont_ID[Loca_contID][0]]; 
             else{
                sum=0.0;
                for (i=0;i<Nseg_tot;i++){
                  if (SegAll_to_Poly[i]==Cont_ID[Loca_contID][0]) {
                      sum += Rho_seg_b[i];
                  }
                }
                param=sum;
                return param;
             }
             break;   

      case CONT_BETAMU_I: 
           if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
                return Betamu_chain[Cont_ID[Loca_contID][0]];
           }
           else                 return Betamu[Cont_ID[Loca_contID][0]];
           break;    

      case CONT_EPSW_I:   
          if (Mix_type==0) return Eps_w[Cont_ID[Loca_contID][0]];
          else             return Eps_ww[Cont_ID[Loca_contID][0]][Cont_ID[Loca_contID][1]];
          break;

      case CONT_EPSWF_IJ: 
          return Eps_wf[Cont_ID[Loca_contID][0]][Cont_ID[Loca_contID][1]]; 
          break;

      case CONT_EPSFF_IJ:   
           return Eps_ff[Cont_ID[Loca_contID][0]][Cont_ID[Loca_contID][1]]; 
           break;

      case CONT_ELECPARAM_I:  
           return Elec_param_w[Cont_ID[Loca_contID][0]]; 
           break;

      case CONT_ELECPARAM_ALL:  
           return Elec_param_w[0];
           break;

      case CONT_SEMIPERM_IJ: 
           return Vext_membrane[Cont_ID[Loca_contID][0]][Cont_ID[Loca_contID][1]]; 
           break;

      default:
        if (cont_type > 99 && cont_type < 199) {
           param = get_init_param_archived_plugin(cont_type,Loca_contID); 
           return param;
        }
        else if (cont_type >= 200 && cont_type<299){
           param = get_init_param_user_plugin(cont_type,Loca_contID); 
           return param;
        }
        else{
           printf("ERROR: Unknown Continuation parameter %d\n",cont_type);
           exit(-1); 
        }
        break;
  }
  return 0.0;
}
/*****************************************************************************/
void assign_parameter_tramonto(int cont_type, double param,int Loca_contID)
/* Note: Post_processing assumes the cont_type flags are the same as those
   used in Tramonto's own continuation */
{
  int i,j,icomp,jcomp,iw,iwall_type,inode;
  double ratio,temp_save,scale_save,eps_wf_save[NCOMP_MAX][NWALL_MAX_TYPE],param_save,rho_chain;
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

      case CONT_RHO_I:   
            if (Type_poly==NONE) Rho_b[Cont_ID[Loca_contID][0]]=param; 
            else{
              rho_chain=0.0;
              for (i=0;i<Nseg_tot; i++){
                  if (SegAll_to_Poly[i]==Cont_ID[Loca_contID][0]) rho_chain+= Rho_seg_b[i];
              }
              ratio=param/rho_chain;

              for (i=0;i<Ntype_mer;i++) Rho_b[i]=0.0;
              for (i=0;i<Nseg_tot; i++){
                 if (SegAll_to_Poly[i]==Cont_ID[Loca_contID][0]) Rho_seg_b[i]*=ratio;
                 Rho_b[Unk2Comp[i]]+=Rho_seg_b[i];
              }
            }

            if (Type_poly == CMS) setup_polymer_cr();
            recalculate_stencils();
            break;

      case CONT_BETAMU_I: 
          if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
                Betamu_chain[Cont_ID[Loca_contID][0]]=param;
          }
          else  Betamu[Cont_ID[Loca_contID][0]]=param;
          break;

      case CONT_EPSW_I:  
          iwall_type=Cont_ID[Loca_contID][0];
          if (Mix_type ==0) {
              Eps_w[iwall_type] = param;
              for (i=0; i<Ncomp; i++){ 
                 for (iw=0; iw<Nwall_type; iw++) eps_wf_save[i][iw]=Eps_wf[i][iw];
              }
              pot_parameters(NULL);
              for (i=0; i<Ncomp; i++){
                 ratio = Eps_wf[i][iwall_type]/eps_wf_save[i][iwall_type];
                 scale_vext_epswf(ratio,i,0); 
              }
          }
          else Eps_ww[iwall_type][iwall_type]=param;
          break;

      case CONT_EPSWF_IJ: 
          icomp=Cont_ID[Loca_contID][0];
          iwall_type=Cont_ID[Loca_contID][1];

          ratio = param/Eps_wf[icomp][iwall_type];
          Eps_wf[icomp][iwall_type]    = param;
          scale_vext_epswf(ratio,icomp,iwall_type); break;

      case CONT_EPSFF_IJ: 
          icomp=Cont_ID[Loca_contID][0];
          jcomp=Cont_ID[Loca_contID][1];
          Eps_ff[icomp][jcomp]=param;  
          if (icomp != jcomp){
             Eps_ff[jcomp][icomp]=param;
          }
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
          if (Type_hsdiam == BH_DIAM){ 
              calc_HS_diams(); 
              calc_InvR_params();
          }
         if (Type_poly == CMS) setup_polymer_cr();
         recalculate_stencils();
         break;

      case CONT_ELECPARAM_I: 
           Elec_param_w[Cont_ID[Loca_contID][0]]=param;
           break;

      case CONT_ELECPARAM_ALL: 
           ratio=1./Elec_param_w[0];
           Elec_param_w[0]=param;
           ratio *= Elec_param_w[0];
           scale_elec_param(ratio); 
           break;

      case CONT_SEMIPERM_IJ:
		  iwall_type=Cont_ID[Loca_contID][0];
		  icomp=Cont_ID[Loca_contID][1];
		  param_save=Vext_membrane[iwall_type][icomp];
			  for (inode=0;inode<Nnodes_per_proc;inode++){
                  if (fabs(Vext[inode][icomp]-param_save)<1.e-10) Vext[inode][icomp]=param;
			  }
			  Vext_membrane[iwall_type][icomp]=param;
		  break;
		  
      default:
        if (cont_type > 99 && cont_type < 199) {
           assign_param_archived_plugin(cont_type,Loca_contID,param); 
        }
        else if (cont_type >= 200 && cont_type<299){
           assign_param_user_plugin(cont_type,Loca_contID,param); 
        }
        else{
           printf("ERROR: Unknown Continuation parameter %d\n",cont_type);
           exit(-1); 
        }
        break;
  }

  /* for most cases...recalculate thermo based on new parameter.  However if
     calculating bulk_coexistence or varying Betamu do not call thermo */
  if (Loca.cont_type1 != CONT_BETAMU_I && !(Loca.method==4 && Loca.cont_type2 == CONT_BETAMU_I) ){
          thermodynamics(output_file1);
   }
}
/*****************************************************************************/
/*print_cont_type: Here print the type of the variable that
                     is changing in a given run
                     Note the logic here should mirror that in print_cont_variable() above
*/
void print_cont_type(int cont_type,FILE *fp,int Loca_contID)
{
  int i,idim,icomp,iwall,nloop,jcomp;

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
         if (Type_attr != NONE){
         for (icomp=0;icomp<Ncomp;icomp++)
            for (jcomp=0;jcomp<Ncomp;jcomp++)
                fprintf(fp,"EPS_ff[%d][%d]  ",icomp,jcomp);
         }
         break;

      case CONT_RHO_I:
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

      case CONT_BETAMU_I:
         if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)
            fprintf(fp,"Betamu_chain[%d]:  ",Cont_ID[Loca_contID][0]);
         else
            fprintf(fp,"Betamu[%d]:  ",Cont_ID[Loca_contID][0]);
         break;

      case CONT_EPSW_I:
         if (Mix_type==0) fprintf(fp,"Eps_w[%d]  ",Cont_ID[Loca_contID][0]);
         else             fprintf(fp,"Eps_ww[%d][%d]  ",Cont_ID[Loca_contID][0],Cont_ID[Loca_contID][1]);
         break;

      case CONT_EPSWF_IJ:
         fprintf(fp,"Eps_wf[%d][%d]:  ",Cont_ID[Loca_contID][0],Cont_ID[Loca_contID][1]); 
         break;

      case CONT_EPSFF_IJ:
         fprintf(fp,"Eps_ff[%d][%d]:  ",Cont_ID[Loca_contID][0],Cont_ID[Loca_contID][1]); 
         break;

      case CONT_ELECPARAM_I:
         fprintf(fp,"Elec_param_w[%d]:  ",Cont_ID[Loca_contID][0]); 
         break;

      case CONT_ELECPARAM_ALL:
         fprintf(fp,"Elec_param_w[0]:  "); 
         break;

     case CONT_SEMIPERM_IJ:
         fprintf(fp, "Vext_membrane[%d][%d]  ",Cont_ID[Loca_contID][0],Cont_ID[Loca_contID][1]); 
         break;

      default:
        if (cont_type > 99 && cont_type < 199) {
           print_cont_type_archived_plugin(cont_type,fp,Loca_contID);
        }
        else if (cont_type >= 200 && cont_type<299){
           print_cont_type_user_plugin(cont_type,fp,Loca_contID); 
        }
        else{
           printf("ERROR: Unknown Continuation parameter %d\n",cont_type);
           exit(-1); 
        }
        break;
   }
   return;
}
/*****************************************************************************/
/*print_cont_variable: Here print the value of the variable that
                     is changing in a given run */
void print_cont_variable(int cont_type,FILE *fp,int Loca_contID)
{                 
   int i,idim,icomp,iwall,iwall_type,nloop,jcomp;
   double kappa,kappa_sq,rhosum,rho_chain;
         

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
         if (Type_attr != NONE){
           for (icomp=0;icomp<Ncomp;icomp++){
              for (jcomp=0;jcomp<Ncomp;jcomp++)
                  fprintf(fp,"%10.7f  ",Eps_ff[icomp][jcomp]);
           }
         }
         break;

      case CONT_RHO_I:
         if (Type_poly==NONE) fprintf(fp,"%11.8f  ",Rho_b[Cont_ID[Loca_contID][0]]);
         else{
              rho_chain=0.0;
              for (i=0;i<Nseg_tot; i++){
                  if (SegAll_to_Poly[i]==Cont_ID[Loca_contID][0]) rho_chain+= Rho_seg_b[i];
              }
              fprintf(fp,"%11.8f  ",rho_chain);
         }

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

      case CONT_BETAMU_I:
         if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)
                  fprintf(fp,"%11.8f   ", Betamu_chain[Cont_ID[Loca_contID][0]]);
         else     fprintf(fp,"%11.8f   ", Betamu[Cont_ID[Loca_contID][0]]);
         break;

      case CONT_EPSW_I:
         if (Mix_type==0) fprintf(fp,"%11.8f   ", Eps_w[Cont_ID[Loca_contID][0]]);
         else fprintf(fp,"%11.8f   ", Eps_ww[Cont_ID[Loca_contID][0]][Cont_ID[Loca_contID][1]]);

         break;

      case CONT_EPSWF_IJ:
         fprintf(fp,"%11.8f   ", Eps_wf[Cont_ID[Loca_contID][0]][Cont_ID[Loca_contID][1]]);
         break;

      case CONT_EPSFF_IJ:
         fprintf(fp,"%11.8f   ", Eps_ff[Cont_ID[Loca_contID][0]][Cont_ID[Loca_contID][1]]);
         break;

      case CONT_ELECPARAM_I:
         fprintf(fp,"%11.8f   ", Elec_param_w[Cont_ID[Loca_contID][0]]);
         break;

      case CONT_ELECPARAM_ALL:
         fprintf(fp,"%11.8f   ", Elec_param_w[0]);
         break;

      case CONT_SEMIPERM_IJ:
         fprintf(fp,"%11.8f   ", Vext_membrane[Cont_ID[Loca_contID][0]][Cont_ID[Loca_contID][1]]); 
         break;

      default:
        if (cont_type > 99 && cont_type < 199) {
           print_cont_variable_archived_plugin(cont_type,fp,Loca_contID);
        }
        else if (cont_type >= 200 && cont_type<299){
           print_cont_variable_user_plugin(cont_type,fp,Loca_contID);
        }
        else{
           printf("ERROR: Unknown Continuation parameter %d\n",cont_type);
           exit(-1);
        }
        break;

   }
   return;

}
/*****************************************************************************/
