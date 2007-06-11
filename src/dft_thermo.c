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
void  thermodynamics(char *output_file1)
{
   char *yo = "thermodynamics";

   if (Proc==0 && Iwrite!=NO_SCREEN){
          printf("\n-------------------------------------------------------------------------------\n");
          printf("%s: Doing Thermo precalculations\n",yo);
   }

    /* first call any functions needed to preprocess global bulk variables associated with the functionals chosen for this run */
    if (L_HSperturbation){
       if (Type_poly == WTC) WTC_thermo_precalc(output_file1);   /* do this first to set up segment densities */
       if (Type_func != NONE) HS_thermo_precalc(output_file1); 
       if (Type_attr != NONE ) ATT_thermo_precalc();
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
    calc_chempot(output_file1);
    calc_pressure(output_file1);

    return;
}
/***************************************************************************************/
/* calc_pressure: this routine contains the logic for assembly of the pressure */
void calc_pressure(char *output_file1)
{
   double betap_hs_DFT,betap_hs_PY,betap_hs_bulk,betap_att;
   double betap_att_LBB, betap_att_RTF,betap_hs_bulk_LBB,betap_hs_bulk_RTF;
   int icomp;
   FILE *fp;

   betap_att = 0.0;
   betap_att_LBB = 0.0;
     betap_att_RTF = 0.0;

   if( (fp = fopen(output_file1,"a+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
   }
   fprintf(fp,"\n!!!!!!!!!!!!! output from dft_thermo.c !!!!!!!!!!!!!!!!!!\n");
 
   if (Lsteady_state){          /* CASE WITH DIFFUSION */
      if (L_HSperturbation){
				/* IDEAL contributions */
          Betap_LBB=pressure_ideal_gas(Rho_b_LBB);
          Betap_RTF=pressure_ideal_gas(Rho_b_RTF);
          

				/* HS FMT contributions */
          if (Type_func!= NONE){
               Betap_LBB += pressure_FMT_hs(Rho_b_LBB,&betap_hs_bulk_LBB);
               Betap_RTF += pressure_FMT_hs(Rho_b_RTF,&betap_hs_bulk_RTF);
          }
				/* MF ATT contributions */
          if (Type_attr != NONE){
	    betap_att_LBB = pressure_att(Rho_b_LBB);
	    Betap_LBB += betap_att_LBB;
	    betap_att_RTF = pressure_att(Rho_b_RTF); 
	    Betap_RTF += betap_att_RTF;
          }
				/* electrostatics contributions */
          if (Type_coul != NONE){
         /* please put call to electrostatic pressure function here */
          }
				/* WTC contributions */
	  /* note these aren't additive,instead we recalculate the HS, ideal terms here */
	  /* must then add in contributions from attractions, Coulomb */
          if (Type_poly ==WTC){
	    Betap_LBB = pressure_WTC(Rho_seg_LBB,betap_hs_bulk_LBB);
	    Betap_LBB += betap_att_LBB;
	    Betap_RTF = pressure_WTC(Rho_seg_RTF,betap_hs_bulk_RTF);
	    Betap_RTF += betap_att_RTF;
          }
       }
       else{ 
           /* CMS pressure with diffusion would go here */
       }
       if (Proc==0){
            if (Iwrite != NO_SCREEN) {
                print_to_screen(Betap_LBB,"Betap_LBB");
                print_to_screen(Betap_RTF,"Betap_RTF");
            }
            print_to_file(fp,Betap_LBB,"Betap_LBB",TRUE);
            print_to_file(fp,Betap_RTF,"Betap_RTF",TRUE);
       }    
   }
   else{          		/* CASE WITH NO DIFFUSION */
      if (L_HSperturbation){
				/* IDEAL contributions */
          Betap=pressure_ideal_gas(Rho_b);

				/* HS FMT contributions */
          if (Type_func != NONE) {
               betap_hs_DFT = pressure_FMT_hs(Rho_b,&betap_hs_bulk);
                   if (Proc==0 && Iwrite != NO_SCREEN) printf("\tDFT HS pressure is %9.6f\n",betap_hs_DFT);
               betap_hs_PY = pressure_PY_hs(Rho_b);
                   if (Proc==0 && Iwrite != NO_SCREEN) printf("\tPY HS pressure is %9.6f\n",betap_hs_PY);
               Betap += betap_hs_DFT;
               for (icomp=0;icomp<Ncomp;icomp++) Betap-=Rho_b[icomp]; /* calculation of betap_hs includes ideal terms */
          }
				/* MF ATT contributions */
          if (Type_attr != NONE){
               betap_att = pressure_att(Rho_b);
	       Betap += betap_att;
          }
				/* electrostatics contributions */
          if (Type_coul != NONE){
         /* please put call to electrostatic pressure function here */
          }
				/* WTC contributions */
	  /* note these aren't additive,instead we recalculate the HS, ideal terms here */
	  /* must then add in contributions from attractions, Coulomb */
          if (Type_poly ==WTC){
	    Betap = pressure_WTC(Rho_seg_b,betap_hs_bulk);
	    Betap += betap_att;
	    /* if (Proc==0 && Iwrite != NO_SCREEN) printf("\tWTC pressure is %9.6f\n",);*/
          }
         if (Proc==0){
              if (Iwrite != NO_SCREEN) {
                  print_to_screen(Betap,"Betap");
              }
              print_to_file(fp,Betap,"Betap",TRUE);
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
void calc_chempot(char *output_file1)
{

   double betamu_hs[NCOMP_MAX];
   int icomp,iseg;
   FILE *fp;

   if( (fp = fopen(output_file1,"a+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
   }
 
   if (Lsteady_state){          /* CASE WITH DIFFUSION */
      if (L_HSperturbation){

				/* IDEAL contributions */
          chempot_ideal_gas(Rho_b_LBB,Betamu_LBB);
          chempot_ideal_gas(Rho_b_RTF,Betamu_RTF);

				/* HS FMT contributions */
          if (Type_func!= NONE){
               chempot_FMT_hs(Rho_b_LBB);
               for (icomp=0; icomp<Ncomp; icomp++) Betamu_LBB[icomp] += Betamu_hs_ex[icomp];    
               chempot_FMT_hs(Rho_b_RTF);
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
             if (Type_coul==DELTAC) {
                 chempot_ELEC_MSA(Rho_b_LBB);
                 for (icomp=0;icomp<Ncomp;icomp++) Betamu_LBB[icomp]+=Deltac_b[icomp];
                 chempot_ELEC_MSA(Rho_b_RTF);
                 for (icomp=0;icomp<Ncomp;icomp++) Betamu_RTF[icomp]+=Deltac_b[icomp];
                 
             }
          }
				/* WTC contributions */
          if (Type_poly ==WTC){
                                 /* note that this should come last because the segment 
                                    chemical potentials are built using the component chemical potentials
                                    for other physics types */
               chempot_WTC(Rho_seg_LBB,Betamu_LBB);
               for (iseg=0;iseg<Nseg_tot;iseg++){ 
                     Betamu_seg_LBB[iseg]=Betamu_seg[iseg];
                     Betamu_wtc_LBB[iseg]=Betamu_wtc_LBB[iseg];
               }
               chempot_WTC(Rho_seg_RTF,Betamu_RTF);
               for (iseg=0;iseg<Nseg_tot;iseg++){ 
                     Betamu_seg_RTF[iseg]=Betamu_seg[iseg];
                     Betamu_wtc_RTF[iseg]=Betamu_wtc_RTF[iseg];
               }
          }
          if (Proc==0){
               if (Type_poly==WTC){
                  if (Iwrite != NO_SCREEN) {
                      for (iseg=0;iseg<Nseg_tot;iseg++) print_to_screen_comp(iseg,Betamu_seg_LBB[iseg],"Betamu_seg_LBB");
                      for (iseg=0;iseg<Nseg_tot;iseg++) print_to_screen_comp(iseg,Betamu_seg_RTF[iseg],"Betamu_seg_RTF");
                  }
                  for (iseg=0;iseg<Nseg_tot;iseg++) print_to_file_comp(fp,iseg,Betamu_seg_LBB[iseg],"Betamu_seg_LBB",TRUE);
                  for (iseg=0;iseg<Nseg_tot;iseg++) print_to_file_comp(fp,iseg,Betamu_seg_RTF[iseg],"Betamu_seg_RTF",TRUE);
               }
               else{
                  if (Iwrite != NO_SCREEN) {
                      for (icomp=0;icomp<Ncomp;icomp++) print_to_screen_comp(icomp,Betamu_LBB[icomp],"Betamu_LBB");
                      for (icomp=0;icomp<Ncomp;icomp++) print_to_screen_comp(icomp,Betamu_RTF[icomp],"Betamu_RTF");
                  }
                  for (icomp=0;icomp<Ncomp;icomp++) print_to_file_comp(fp,icomp,Betamu_LBB[icomp],"Betamu_LBB",TRUE);
                  for (icomp=0;icomp<Ncomp;icomp++) print_to_file_comp(fp,icomp,Betamu_RTF[icomp],"Betamu_RTF",TRUE);
               }
          }    
       }
       else{ 
           /* CMS chemical potentials with diffusion would go here */
       }
   }
   else{          		/* CASE WITH NO DIFFUSION */
      if (L_HSperturbation){
				/* IDEAL contributions */
          chempot_ideal_gas(Rho_b,Betamu);
          for (icomp=0;icomp<Ncomp;icomp++) Betamu_id[icomp]=Betamu[icomp];

				/* HS FMT contributions */
          if (Type_func != NONE) {
               chempot_FMT_hs(Rho_b);
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
              if (Type_coul==DELTAC) {
                 chempot_ELEC_MSA(Rho_b);
                 for (icomp=0;icomp<Ncomp;icomp++) Betamu[icomp]+=Deltac_b[icomp];
              }
          }
				/* WTC contributions */
          if (Type_poly==WTC){   
                                 /* note that this should come last because the segment 
                                    chemical potentials are built using the component chemical potentials
                                    for other physics types */
               chempot_WTC(Rho_seg_b,Betamu);
          }
          if (Proc==0){
               if (Type_poly==WTC){
                  if (Iwrite != NO_SCREEN) for (iseg=0;iseg<Nseg_tot;iseg++) print_to_screen_comp(iseg,Betamu_seg[iseg],"Betamu_seg");
                  for (iseg=0;iseg<Nseg_tot;iseg++) print_to_file_comp(fp,iseg,Betamu_seg[iseg],"Betamu_seg",TRUE);
               }
               else{
                  if (Iwrite != NO_SCREEN) for (icomp=0;icomp<Ncomp;icomp++) print_to_screen_comp(icomp,Betamu[icomp],"Betamu");
                  for (icomp=0;icomp<Ncomp;icomp++) print_to_file_comp(fp,icomp,Betamu[icomp],"Betamu",TRUE);
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
