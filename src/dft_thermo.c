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

#include "dft_thermo.h"

/***************************************************************************************/
void  thermodynamics(char *output_file1,int iwrite)
{
   char *yo = "thermodynamics";
   double scale_fac_tmp[NCOMP_MAX][NCOMP_MAX];
   int pol_num,icomp;
   if (Proc==0 && iwrite!=NO_SCREEN){
          printf("\n-------------------------------------------------------------------------------\n");
          printf("%s: Doing Thermo precalculations\n",yo);

   }

    /* first call any functions needed to preprocess global bulk variables associated with the functionals chosen for this run */

    if (L_HSperturbation){
                                                                    /* set up segment densities */
       if (Type_poly == WTC || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) WTC_thermo_precalc(output_file1);   
                                                                    /* set up bulk WJDC_field and G values */
       if (Type_func != NONE)  HS_thermo_precalc(output_file1); 

       if (Type_attr != NONE ) ATT_thermo_precalc();
       if (Type_poly == WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) {
            /* the Scale_fac must be set to zero on the first call since the scaling factor is not used in computation
               of the bulk properties. */
            for (icomp=0; icomp<Ncomp;icomp++){
              for (pol_num=0; pol_num<Npol_comp;pol_num++) {
                 scale_fac_tmp[pol_num][icomp]=Scale_fac_WJDC[pol_num][icomp];
                 Scale_fac_WJDC[pol_num][icomp]=0.0;
              }
            }
            WJDC_thermo_precalc(output_file1);
            /* Return the Scale_fac array to its previous values
               of the bulk properties. */
            for (icomp=0; icomp<Ncomp;icomp++)
              for (pol_num=0; pol_num<Npol_comp;pol_num++) { Scale_fac_WJDC[pol_num][icomp]=scale_fac_tmp[pol_num][icomp]; } 
       }
       /* if (Type_coul == DELTAC)  nothing to do here... */
    }
    /*else { } no CMS precalculations implemented here yet */

    /* if a two phase coexistence calculation is desired run it now.
       Note that we may want to run the DFT calculation at coexistence (for a
       wetting study for example) or we may simply want to know where we are 
       relative to two phase coexistence.  To reimplement this capability we need
       a redesign.  The old implementation for two species can be found in the source_archives
       directory. */

    /*if (Lcoexistence){*/
       /* compute and print coexistence properties */
/*    }
    if (LDFT_at_coex){*/
       /* adjust bulk densities before final pressure and chempot calculations*/
   /* }*/

    /* compute bulk pressure and chemical potentials for the DFT calculation */
    /* must calculate chemical potentials first because we use mu in the WTC calculation of the pressure */
    calc_chempot(output_file1,iwrite);
    calc_pressure(output_file1,iwrite);
    /*adjust bulk terms for scaling */
    /*if  (Physics_scaling && Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) WJDC_thermo_precalc(output_file1); */
    if  (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) WJDC_thermo_precalc(output_file1); 

    return;
}
/***************************************************************************************/
/* calc_pressure: this routine contains the logic for assembly of the pressure */
void calc_pressure(char *output_file1,int iwrite)
{
   double betap_hs_DFT,betap_att,betap_elec,betap_chain;
   /*double betap_hs_PY;*/
   double betap_att_LBB, betap_att_RTF;
   double betap_elec_LBB, betap_elec_RTF;
   FILE *fp;

   betap_att = 0.0;
   betap_att_LBB = 0.0;
     betap_att_RTF = 0.0;

   if( (fp = fopen(output_file1,"a+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
   }
 
   if (Type_interface!=UNIFORM_INTERFACE){      
      if (L_HSperturbation){
				/* IDEAL contributions */
          if (Lseg_densities){
             Betap_LBB=pressure_ideal_gas(Rho_seg_LBB);
             Betap_RTF=pressure_ideal_gas(Rho_seg_RTF);
          }
          else{
             Betap_LBB=pressure_ideal_gas(Rho_b_LBB);
             Betap_RTF=pressure_ideal_gas(Rho_b_RTF);
          }

				/* HS FMT contributions */
          if (Type_func!= NONE){
               Betap_LBB += pressure_FMT_hs(Rhobar_b_LBB,Dphi_Drhobar_LBB);
               Betap_RTF += pressure_FMT_hs(Rhobar_b_RTF,Dphi_Drhobar_RTF);
          }
				/* MF ATT contributions */
          if (Type_attr != NONE){
	    betap_att_LBB = pressure_att(Rho_b_LBB);
	    Betap_LBB += betap_att_LBB;
	    betap_att_RTF = pressure_att(Rho_b_RTF); 
	    Betap_RTF += betap_att_RTF;
          }
				/* electrostatics contributions */
          if (Type_coul == DELTAC_RPM || Type_coul == DELTAC_GENERAL){
            betap_elec_LBB = pressure_elec_MSA(Rho_b_LBB);
            Betap_LBB += betap_elec_LBB;
            betap_elec_RTF = pressure_elec_MSA(Rho_b_RTF);
            Betap_RTF += betap_elec_RTF;
          }
				/* WTC contributions */
	  /* note these aren't additive,instead we recalculate the HS, ideal terms here */
	  /* must then correct contributions from attractions, Coulomb */
	  /* note chem. potential term gives twice the att. pressure, so subtract att. pressure term here */
          if (Type_poly == WTC || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
	    betap_chain = pressure_WTC(Rho_seg_LBB,Xi_cav_LBB);
                   if (Proc==0 && iwrite != NO_SCREEN) printf("\t LBB chain pressure is %9.6f\n",betap_chain);
            Betap_LBB += betap_chain;
	    betap_chain = pressure_WTC(Rho_seg_RTF,Xi_cav_RTF);
                   if (Proc==0 && iwrite != NO_SCREEN) printf("\t RTF chain pressure is %9.6f\n",betap_chain);
            Betap_RTF += betap_chain;
          }
       }
       else{ 
           /* CMS pressure with diffusion would go here */
       }
       if (Proc==0){
            if (iwrite != NO_SCREEN) {
                print_to_screen(Betap_LBB,"Betap_LBB");
                print_to_screen(Betap_RTF,"Betap_RTF");
            }
            print_to_file(fp,Betap_LBB,"Betap_LBB",2);
            print_to_file(fp,Betap_RTF,"Betap_RTF",2);
       }    
   }
   else if (Type_interface==UNIFORM_INTERFACE){ 
      if (L_HSperturbation){
				/* IDEAL contributions */
          if (Lseg_densities)  Betap=pressure_ideal_gas(Rho_seg_b);
          else                 Betap=pressure_ideal_gas(Rho_b);
                   if (Proc==0 && iwrite != NO_SCREEN) printf("\t ideal gas pressure is %9.6f\n",Betap);

				/* HS FMT contributions */
          if (Type_func != NONE) {
               betap_hs_DFT = pressure_FMT_hs(Rhobar_b,Dphi_Drhobar_b);
                   if (Proc==0 && iwrite != NO_SCREEN) printf("\tDFT HS pressure is %9.6f\n",betap_hs_DFT);
               /*betap_hs_PY = pressure_PY_hs(Rho_b);
                   if (Proc==0 && iwrite != NO_SCREEN) printf("\tPY HS pressure is %9.6f\n",betap_hs_PY);*/
               Betap += betap_hs_DFT;
          }
				/* MF ATT contributions */
          if (Type_attr != NONE){
               betap_att = pressure_att(Rho_b);
                   if (Proc==0 && iwrite != NO_SCREEN) printf("\t att pressure is %9.6f\n",betap_att);
	       Betap += betap_att;
          }
				/* electrostatics contributions */
          if (Type_coul == DELTAC_RPM || Type_coul==DELTAC_GENERAL){
               betap_elec = pressure_elec_MSA(Rho_b);
                   if (Proc==0 && Iwrite != NO_SCREEN) printf("\t elec_MSA pressure is %9.6f\n",betap_elec);
               Betap += betap_elec;
          }
				/* WTC contributions */
	  /* note these aren't additive,instead we recalculate the HS, ideal terms here */
	  /* must then correct contributions from attractions, Coulomb */
          if (Type_poly == WTC || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
	    betap_chain = pressure_WTC(Rho_seg_b,Xi_cav_b);
                   if (Proc==0 && iwrite != NO_SCREEN) printf("\t chain pressure is %9.6f\n",betap_chain);
            Betap += betap_chain;
          }
         if (Proc==0){
              if (iwrite != NO_SCREEN) {
                  print_to_screen(Betap,"Betap");
              }
              print_to_file(fp,Betap,"Betap",2);
         }    
      }
      else{
         /* put a calculation of CMS pressure here */
      }
   }
   fclose(fp);
   return;
}
/***************************************************************************************/
/* calc_chempot: this routine contains the logic for assembly of the chemical potentials */
void calc_chempot(char *output_file1,int iwrite)
{

   int icomp,iseg,ipol;
   FILE *fp;

   if( (fp = fopen(output_file1,"a+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
   }

   if (Type_interface !=UNIFORM_INTERFACE){          /* CASE WITH DIFFUSION */
      if (L_HSperturbation){
          if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){  /* this is different than all the others because we compute
                                    chain chemical potentials.  Note that a segment chemical potential
                                    should also be implemented. */
             chempot_chain_wjdc(Rho_seg_LBB,Betamu_chain_LBB,Field_WJDC_LBB,G_WJDC_LBB);
             chempot_chain_wjdc(Rho_seg_RTF,Betamu_chain_RTF,Field_WJDC_RTF,G_WJDC_RTF);
             if (Type_interface==PHASE_INTERFACE){
                for (icomp=0;icomp<Npol_comp;icomp++) Betamu_chain[icomp]=0.5*(Betamu_chain_LBB[icomp]+Betamu_chain_RTF[icomp]);
             }
          }
/*          else{*/

				/* IDEAL contributions */
          chempot_ideal_gas(Rho_b_LBB,Betamu_LBB);
          chempot_ideal_gas(Rho_b_RTF,Betamu_RTF);

				/* HS FMT contributions */
          if (Type_func!= NONE){
               chempot_FMT_hs(Dphi_Drhobar_LBB);
               for (icomp=0; icomp<Ncomp; icomp++) Betamu_LBB[icomp] += Betamu_hs_ex[icomp];    
               chempot_FMT_hs(Dphi_Drhobar_RTF);
               for (icomp=0; icomp<Ncomp; icomp++) Betamu_RTF[icomp] += Betamu_hs_ex[icomp];    
          }
				/* MF ATT contributions */
          if (Type_attr != NONE){
               chempot_att(Rho_b_LBB);
               for (icomp=0; icomp<Ncomp; icomp++) Betamu_LBB[icomp] += Betamu_att[icomp];
               chempot_att(Rho_b_RTF);
               for (icomp=0; icomp<Ncomp; icomp++) Betamu_RTF[icomp] += Betamu_att[icomp];
          }
				/* electrostatics contributions */
          if (Type_coul != NONE){
             for (icomp=0;icomp<Ncomp;icomp++){ 
               Betamu_LBB[icomp] += Charge_f[icomp]*(Elec_pot_LBB);
               Betamu_RTF[icomp] += Charge_f[icomp]*(Elec_pot_RTF);
             }
             if (Type_coul==DELTAC_RPM || Type_coul==DELTAC_GENERAL) {
                 printf("can't do inhomogeneous boundaries with constant Delta c (electrostatic) corrections !!\n");
                 if (Type_coul==DELTAC_RPM)          chempot_ELEC_MSA_RPM(Rho_b_LBB);
                 else if (Type_coul==DELTAC_GENERAL) chempot_ELEC_MSA_GENERAL(Rho_b_LBB);
                 for (icomp=0;icomp<Ncomp;icomp++){ 
                    Betamu_LBB[icomp]+=Deltac_b[icomp];
                    printf("chemical potential contribution on left is icomp=%d mu_deltac=%g\n",icomp,Deltac_b[icomp]);
                 }

                 if (Type_coul==DELTAC_RPM)          chempot_ELEC_MSA_RPM(Rho_b_RTF);
                 else if (Type_coul==DELTAC_GENERAL) chempot_ELEC_MSA_GENERAL(Rho_b_RTF);
                 for (icomp=0;icomp<Ncomp;icomp++){
                       Betamu_RTF[icomp]+=Deltac_b[icomp];
                       printf("chemical potential contribution on right is icomp=%d mu_deltac=%g\n",icomp,Deltac_b[icomp]);
                 }
                 printf("with the bulk electrostatics correction models (Type_coul=%d),  setting the chemical \n",Type_coul);
                 printf("potential between the limits on the left and right sides is not well defined.\n");
                 printf("Exiting code now.\n");
                 exit(-1);
                 
             }
          }
          if (Type_interface==PHASE_INTERFACE) {
             for (icomp=0;icomp<Ncomp;icomp++) Betamu[icomp]=0.5*(Betamu_RTF[icomp]+Betamu_LBB[icomp]);
          }
				/* WTC contributions */
          if (Type_poly ==WTC || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
                                 /* note that this should come last because the segment 
                                    chemical potentials are built using the component chemical potentials
                                    for other physics types */
               chempot_WTC(Rho_seg_LBB,Betamu_LBB,Xi_cav_LBB);
               for (iseg=0;iseg<Nseg_tot;iseg++){ 
                     Betamu_seg_LBB[iseg]=Betamu_seg[iseg];
                     Betamu_wtc_LBB[iseg]=Betamu_wtc_LBB[iseg];
               }
               for (ipol=0;ipol<Npol_comp;ipol++){
                  Betamu_chain_LBB[ipol]=0.0;
                  for (iseg=0;iseg<Nmer[ipol];iseg++){
                     Betamu_chain_LBB[ipol]+= Betamu_seg_LBB[SegChain2SegAll[ipol][iseg]];
                  }
               }
               chempot_WTC(Rho_seg_RTF,Betamu_RTF,Xi_cav_RTF);
               for (iseg=0;iseg<Nseg_tot;iseg++){ 
                     Betamu_seg_RTF[iseg]=Betamu_seg[iseg];
                     Betamu_wtc_RTF[iseg]=Betamu_wtc_RTF[iseg];
               }
               for (ipol=0;ipol<Npol_comp;ipol++){
                  Betamu_chain_RTF[ipol]=0.0;
                  for (iseg=0;iseg<Nmer[ipol];iseg++){
                     Betamu_chain_RTF[ipol]+= Betamu_seg_RTF[SegChain2SegAll[ipol][iseg]];
                  }
               }
          }
          if (Proc==0){
               if (Lseg_densities){
                  if (Type_poly==WTC){
                  if (iwrite != NO_SCREEN) {
                      for (iseg=0;iseg<Nseg_tot;iseg++) print_to_screen_comp(iseg,Betamu_seg_LBB[iseg],"Betamu_seg_LBB");
                      for (iseg=0;iseg<Nseg_tot;iseg++) print_to_screen_comp(iseg,Betamu_seg_RTF[iseg],"Betamu_seg_RTF");
                  }
                  for (iseg=0;iseg<Nseg_tot;iseg++) print_to_file_comp(fp,iseg,Betamu_seg_LBB[iseg],"Betamu_seg_LBB",2);
                  for (iseg=0;iseg<Nseg_tot;iseg++) print_to_file_comp(fp,iseg,Betamu_seg_RTF[iseg],"Betamu_seg_RTF",2);
                  }
                  if (iwrite != NO_SCREEN){
                      for (ipol=0;ipol<Npol_comp;ipol++) print_to_screen_comp(ipol,Betamu_chain_LBB[ipol],"Betamu_chain_LBB");
                      for (ipol=0;ipol<Npol_comp;ipol++) print_to_screen_comp(ipol,Betamu_chain_RTF[ipol],"Betamu_chain_RTF");
                  }
                  for (ipol=0;iseg<Npol_comp;ipol++) print_to_file_comp(fp,ipol,Betamu_chain_LBB[ipol],"Betamu_chain_LBB",2);
                  for (ipol=0;iseg<Npol_comp;ipol++) print_to_file_comp(fp,ipol,Betamu_chain_RTF[ipol],"Betamu_chain_RTF",2);
               }
               else{
                  if (iwrite != NO_SCREEN) {
                      for (icomp=0;icomp<Ncomp;icomp++) print_to_screen_comp(icomp,Betamu_LBB[icomp],"Betamu_LBB");
                      for (icomp=0;icomp<Ncomp;icomp++) print_to_screen_comp(icomp,Betamu_RTF[icomp],"Betamu_RTF");
                  }
                  for (icomp=0;icomp<Ncomp;icomp++) print_to_file_comp(fp,icomp,Betamu_LBB[icomp],"Betamu_LBB",2);
                  for (icomp=0;icomp<Ncomp;icomp++) print_to_file_comp(fp,icomp,Betamu_RTF[icomp],"Betamu_RTF",2);
               }
          }    
       /*}*/
       }
       else{ 
           /* CMS chemical potentials with diffusion would go here */
       }
   }
   else if(Type_interface==UNIFORM_INTERFACE){          	
      if (L_HSperturbation){

          if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){  /* this is different than all the others because we compute
                                    chain chemical potentials.  Note that a segment chemical potential
                                    should also be implemented. */
             chempot_chain_wjdc(Rho_seg_b,Betamu_chain,Field_WJDC_b,G_WJDC_b);
          }
				/* IDEAL contributions */
          chempot_ideal_gas(Rho_b,Betamu);
          for (icomp=0;icomp<Ncomp;icomp++) Betamu_id[icomp]=Betamu[icomp];

				/* HS FMT contributions */
          if (Type_func != NONE) {
               chempot_FMT_hs(Dphi_Drhobar_b);
               for (icomp=0; icomp<Ncomp; icomp++){ Betamu[icomp] += Betamu_hs_ex[icomp]; }
          }
				/* MF ATT contributions */
          if (Type_attr != NONE){
               chempot_att(Rho_b);
               for (icomp=0; icomp<Ncomp; icomp++) Betamu[icomp] += Betamu_att[icomp];
          }
				/* electrostatics contributions */
          if (Type_coul != NONE){
              /* note that the current assumption is that electrostatic potential is zero in the bulk */
              /* also note that Betamu_deltaC = Deltac_b (for the case of the RPM/MSA electrolyte) */
              if (Type_coul==DELTAC_RPM) {
                 chempot_ELEC_MSA_RPM(Rho_b);
                 for (icomp=0;icomp<Ncomp;icomp++) Betamu[icomp]+=Deltac_b[icomp];
              }
              else if (Type_coul==DELTAC_GENERAL){
                 chempot_ELEC_MSA_GENERAL(Rho_b);
                 for (icomp=0;icomp<Ncomp;icomp++) Betamu[icomp]+=Deltac_b[icomp];
              }
          }
				/* WTC contributions */
          if (Type_poly==WTC || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){   
                                 /* note that this should come last because the segment 
                                    chemical potentials are built using the component chemical potentials
                                    for other physics types */
              if ((Physics_scaling && (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)) || Type_poly==WTC) chempot_WTC(Rho_seg_b,Betamu,Xi_cav_b);
               if (Type_poly==WTC){
                  for (ipol=0;ipol<Npol_comp;ipol++){
                     Betamu_chain[ipol]=0.0;
                     for (iseg=0;iseg<Nmer[ipol];iseg++){
                        Betamu_chain[ipol]+= Betamu_seg[SegChain2SegAll[ipol][iseg]];
                     }
                  }
               }
          }
/************ Optimizing Scale_fac Terms !!! ************/
/*          for (ipol=0;ipol<Npol_comp;ipol++){
              for (icomp=0; icomp<Ncomp; icomp++){
                if (Nseg_type_pol[ipol][icomp] > 0) Scale_fac_WJDC[ipol][icomp]=Betamu_chain[ipol]/Nseg_type_pol[ipol][icomp];
                if (fabs(Scale_fac_WJDC[ipol][icomp])>4.0) Scale_fac_WJDC[ipol][icomp]=4.*Scale_fac_WJDC[ipol][icomp]/fabs(Scale_fac_WJDC[ipol][icomp]);
                printf("modified Scale_fac_WJDC[ipol=%d][icomp=%d]=%g\n",ipol,icomp,Scale_fac_WJDC[ipol][icomp]);
              }
           }*/
/************ End Optimizing Scale_fac Terms !!! ************/
          if (Proc==0){
               if (Lseg_densities){
                  if (iwrite != NO_SCREEN) for (iseg=0;iseg<Nseg_tot;iseg++) print_to_screen_comp(iseg,Betamu_seg[iseg],"Betamu_seg");
                  for (iseg=0;iseg<Nseg_tot;iseg++) print_to_file_comp(fp,iseg,Betamu_seg[iseg],"Betamu_seg",2);
               }
               if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
                  if (iwrite != NO_SCREEN)printf("\n");
                  if (iwrite != NO_SCREEN) for (ipol=0;ipol<Npol_comp;ipol++) print_to_screen_comp(ipol,Betamu_chain[ipol],"Betamu_chain");
                  for (ipol=0;ipol<Npol_comp;ipol++) print_to_file_comp(fp,ipol,Betamu_chain[ipol],"Betamu_chain",2);

                 
               }
               else{
                  if (iwrite != NO_SCREEN) for (icomp=0;icomp<Ncomp;icomp++) print_to_screen_comp(icomp,Betamu[icomp],"Betamu");
                  for (icomp=0;icomp<Ncomp;icomp++) print_to_file_comp(fp,icomp,Betamu[icomp],"Betamu",2);
               }
          }    
      }
      else{
         /* put a calculation of CMS chemical potentials here */
      }
   }
   fclose(fp);
   return;
}
/***************************************************************************************/
