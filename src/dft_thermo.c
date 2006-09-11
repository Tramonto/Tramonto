/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

/* ---------------------------------------------------------
Calculate the relevant thermodynamic properties of the bulk fluid.
We use the Percus-Yevick compressibilty equation of state.

Input:  bulk densities of each component: Rho_b[icomp]

Output:	1) hard sphere pressure:  p sigma_ff[1]^3 / kT
        2) excess hard sphere chemical potentials: mu_hs_ex/kT
        3)       hard sphere chemical potentials: mu_hs/kT
------------------------------------------------------------*/
#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

double   coexistence();
void print_thermo(char *,double,double *);

void  thermodynamics( char *output_file1, int print_flag)
{
   int icomp,jcomp,kcomp,i,ipol,iseg,ibond,jseg,kseg;
   double betap_hs=0., betap_hs_tmp=0.,betamu_hs[NCOMP_MAX],betap_att,p_coex=0.0,betap_TC=0.0;
   double y,dydxi2,dydxi3,sig2,sig3;
   char *yo = "thermodynamics";

 if (!Sten_Type[POLYMER_CR]){
   if (Proc==0 && print_flag &&Iwrite!=NO_SCREEN){
          printf("\n-------------------------------------------------------------------------------\n");
          printf("%s: Doing Thermo precalculations\n",yo);
   }

   if (Type_poly==WTC){
      for (iseg=0;iseg<Nseg_tot;iseg++) Betamu_wtc[iseg]=0.0;

      /* adjust theory for overlapping bonded spheres */
      for (icomp=0;icomp<Ncomp;icomp++){
        for (jcomp=0;jcomp<Ncomp;jcomp++){
          if (Sigma_ff[icomp][icomp]>=Sigma_ff[jcomp][jcomp]){
            if (Bond_ff[icomp][jcomp] >= 0.5*(Sigma_ff[icomp][icomp]+Sigma_ff[jcomp][jcomp])){
               Fac_overlap[icomp][jcomp]=1.0;
            }
            else if (Bond_ff[icomp][jcomp] <= 0.5*fabs(Sigma_ff[icomp][icomp]-Sigma_ff[jcomp][jcomp])){
               Fac_overlap[icomp][jcomp]=0.0;
            }
            else{
               Fac_overlap[icomp][jcomp]=(Bond_ff[icomp][jcomp]/Sigma_ff[jcomp][jcomp]
                                          -0.5*(Sigma_ff[icomp][icomp]-Sigma_ff[jcomp][jcomp])/Sigma_ff[jcomp][jcomp]);
            }
          }
        }    
      }
      Fac_overlap[1][1]=Fac_overlap[0][1]*Fac_overlap[0][1];

      for (icomp=0;icomp<Ncomp;icomp++){
        Fac_overlap_hs[icomp]=0.0;
        for (jcomp=0;jcomp<Ncomp;jcomp++){
          if (Sigma_ff[icomp][icomp]<Sigma_ff[jcomp][jcomp]) Fac_overlap[icomp][jcomp]=Fac_overlap[jcomp][icomp];
          if (Fac_overlap[icomp][jcomp]> Fac_overlap_hs[icomp]) Fac_overlap_hs[icomp]=Fac_overlap[icomp][jcomp];
        }
        }
   }
   else{
      for (icomp=0;icomp<Ncomp;icomp++) Fac_overlap_hs[icomp]=1.0;
   }

   /* use this to turn off all of the overlap nonsense */
Fac_overlap_hs[0]=1.;
Fac_overlap_hs[1]=1.;
    Fac_overlap[0][0]=1.0;
    Fac_overlap[1][0]=1.0;
    Fac_overlap[0][1]=1.0;
    Fac_overlap[1][1]=1.0;


   /* Find bulk coexistence  for a very special case of only one atomistic component with
      mean field attractions turned on !! */

/*   if (Ncomp == 1 && Ipot_ff_n==2 && Iliq_vap != -2){ 
      p_coex = coexistence();
      if (Proc==0 &&Iwrite!=NO_SCREEN) {
         if (print_flag) {
           printf("For kT/epsilon=:%9.6f  the coexistence densities are:\n",1.0/Eps_ff[0][0]);
           printf("\t\t\tliquid: %9.6f  vapor: %9.6f\n",Rho_coex[0],Rho_coex[1]);
           printf("and the pressure at coexistence is:%9.6f\n",p_coex);
           printf("\n\n your density is %9.6f\n",Rho_b[0]);
         }
         else {
           printf("\n\tThermodynamics called for For kT/epsilon=:%9.6f and bulk density %9.6f\n",
                   1.0/Eps_ff[0][0], Rho_b[0]);
         }
      }

      if (Iliq_vap == 1 ){
         if (Proc==0 && print_flag&&Iwrite!=NO_SCREEN) printf("\n RESETTING RHO_BULK = RHO_VAPOR\n");
         Rho_b[0] = Rho_coex[1];
      }
      else if (Iliq_vap == 2 ||Iliq_vap==3){
         if (Proc==0 && print_flag&&Iwrite!=NO_SCREEN) printf("\n RESETTING RHO_BULK = RHO_LIQUID\n");
         Rho_b[0] = Rho_coex[0];
      }
      if (Rho_b[0] < Rho_coex[0] && Rho_b[0] > Rho_coex[1] && 
          Iliq_vap == -1 && Loca.method == -1){
           printf("PROBLEMS: The density you've selected is\n");
           printf("          in the coexistence region\n ");
           printf("          setting density to rho_vapor    \n ");
           Rho_b[0] = Rho_coex[1];
      }
   }*/

  /* for lack of a better place, precalc these quantities here */
 
  if (Ipot_ff_n != IDEAL_GAS) {

     Inv_4pi = 1.0 / (4.0*PI);
     for (icomp=0; icomp<Ncomp; icomp++) {
     if (Sigma_ff[icomp][icomp] < 1.e-6){
        printf("ERROR : icomp=%d Sigma_ff[%d][%d]=%9.6f --- check particle sizes!\n",
                    icomp,icomp,icomp,Sigma_ff[icomp][icomp]);
         exit(-1);
     }
        Inv_rad[icomp] = 2.0 / HS_diam[icomp];
        Inv_4pir[icomp] = Inv_4pi * Inv_rad[icomp];
        Inv_4pirsq[icomp] = Inv_4pir[icomp] * Inv_rad[icomp];
     }
  }

  if (Type_poly==WTC){
     for (iseg=0;iseg<Nseg_tot;iseg++){
         Rho_seg_b[iseg]=Rho_b[Unk2Comp[iseg]]/(double)Nmer_t_total[Unk2Comp[iseg]];
printf("iseg=%d  icomp=%d Nmer_t_total=%d Rho_b=%9.6f Rho_seg_b=%9.6f\n",
  iseg,Unk2Comp[iseg],Nmer_t_total[Unk2Comp[iseg]],Rho_b[Unk2Comp[iseg]],Rho_seg_b[iseg]);
         if (Lsteady_state){
            Rho_seg_LBB[iseg]=Rho_b_LBB[Unk2Comp[iseg]]/Nmer_t_total[Unk2Comp[iseg]];
            Rho_seg_RTF[iseg]=Rho_b_RTF[Unk2Comp[iseg]]/Nmer_t_total[Unk2Comp[iseg]];
         }
     }
     compute_bulk_nonlocal_wtc_properties(output_file1);
  }


  if(Type_func != NONE){
      compute_bulk_nonlocal_hs_properties(output_file1);
  }


  /* Now find pressure and chemical potential contributions for the fluid of interest */

  /* start with ideal gas contributions - need this for all fluids. */ 

   Betap=calc_ideal_gas(Rho_b,Betamu);
   for (icomp=0;icomp<Ncomp;icomp++) Betamu_id[icomp]=Betamu[icomp];
   if (Lsteady_state){
      Betap_LBB=calc_ideal_gas(Rho_b_LBB,Betamu_LBB);
      Betap_RTF=calc_ideal_gas(Rho_b_RTF,Betamu_RTF);
      if (Type_coul !=NONE){
             Betamu_LBB[icomp] += Charge_f[icomp]*(Elec_pot_LBB);
             Betamu_RTF[icomp] += Charge_f[icomp]*(Elec_pot_RTF);
      }
   }

   /* now add in the hard sphere contributions */
  
   if (Type_func!= NONE){


      if (Lsteady_state){
         betap_hs = calc_hs_properties(betamu_hs,Rho_b_LBB);
         for (icomp=0; icomp<Ncomp; icomp++) Betamu_LBB[icomp] += betamu_hs[icomp];    

         betap_hs = calc_hs_properties(betamu_hs,Rho_b_RTF);
         for (icomp=0; icomp<Ncomp; icomp++) Betamu_RTF[icomp] += betamu_hs[icomp];
      }
      else{
         betap_hs = calc_hs_properties(betamu_hs,Rho_b);
if (Proc==0) printf("old pressure calculates %9.6f\n",betap_hs);
         betap_hs_tmp = calc_hs_properties_new(betamu_hs,Rho_b);
if (Proc==0) printf("new pressure calculates %9.6f\n",betap_hs_tmp);
         Betap += betap_hs;
         for (icomp=0; icomp<Ncomp; icomp++){ 
             Betamu[icomp] += betamu_hs[icomp];
             Betap-=Rho_b[icomp]; /* note that the calculation of betap_hs apparantly includes ideal terms .... check this */
         }
     }
   }

printf("start mean field attractions calculation\n");
     /* now add in the mean field attraction contributions */
   if (Type_attr != NONE){
      if (Lsteady_state){
         betap_att = calc_att_properties(Betamu_att,Rho_b_LBB);
         for (icomp=0; icomp<Ncomp; icomp++) 
               Betamu_LBB[icomp] += Betamu_att[icomp];

         betap_att = calc_att_properties(Betamu_att,Rho_b_RTF);
         for (icomp=0; icomp<Ncomp; icomp++) 
               Betamu_RTF[icomp] += Betamu_att[icomp];
      }
      else{
         betap_att = calc_att_properties(Betamu_att,Rho_b);
         Betap += betap_att;
         for (icomp=0; icomp<Ncomp; icomp++) Betamu[icomp] += Betamu_att[icomp];
      }
   }

printf("start mean coulomb calculation\n");
     /* now add in an applied electric field if one exists */
   if (Lsteady_state && Type_coul != NONE){
      for (icomp=0; icomp<Ncomp; icomp++){
         Betamu_LBB[icomp] += Charge_f[icomp]*Elec_pot_LBB;
         Betamu_RTF[icomp] += Charge_f[icomp]*Elec_pot_RTF;
      }
   }

printf("start WTC calculations\n");
    /* now add in the WTC polymer contributions */
   if (Type_poly==WTC){
      if (Lsteady_state){
          printf("have not implemented the thermodynamics for diffusion plus Wertheim-Tripathi-Chapman functionals yet....please debug.\n");
          exit(-1);
      }
      else{

         /* first compute terms that are based on bond pairs starting from segment iseg*/
         /* note that in the bulk the BondWTC is identical to the site density of the jth segment */
         for (iseg=0; iseg<Nseg_tot;iseg++){
            icomp=Unk2Comp[iseg];
            Betamu_seg[iseg]=Betamu[icomp]-log(Nmer_t_total[icomp]); /* correct the ideal gas term for segments */
            for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
                jseg=Bonds_SegAll[iseg][ibond];
                jcomp=Unk2Comp[jseg];
                y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                  Betamu_seg[iseg] += 0.5*(1.0-Fac_overlap[icomp][jcomp]*log(y)-log(Rho_seg_b[jseg])  
                                      - Rho_seg_b[jseg]/Rho_seg_b[iseg]);
                  Betamu_wtc[iseg] += 0.5*(1.0-Fac_overlap[icomp][jcomp]*log(y)-log(Rho_seg_b[jseg])     
                                       -Rho_seg_b[jseg]/Rho_seg_b[iseg]);
            }
         }
         /* now add in the term that involves _all_ the bond pairs in the system for each segment via that cavity
            correlation function.  Note that this term is nonzero even for bond pairs on different polymer
            components than the one where the iseg segment is found ! */
         for (kseg=0; kseg<Nseg_tot;kseg++){
           kcomp = Unk2Comp[kseg];
           sig2=Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp];
           sig3=Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp];
           for (iseg=0; iseg<Nseg_tot;iseg++){
              icomp=Unk2Comp[iseg];
              for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
                jseg=Bonds_SegAll[iseg][ibond];
                jcomp=Unk2Comp[jseg];
                y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                dydxi2 = dy_dxi2_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                dydxi3 = dy_dxi3_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                Betamu_seg[kseg] -= Fac_overlap[icomp][jcomp]*(PI/12.0)*(Rho_seg_b[iseg]/y)*(dydxi2*sig2+dydxi3*sig3);
                Betamu_wtc[kseg] -= Fac_overlap[icomp][jcomp]*(PI/12.0)*(Rho_seg_b[iseg]/y)*(dydxi2*sig2+dydxi3*sig3);
              }
           }
         }
      }
   }

   if (Ipot_ff_c == COULOMB && Sten_Type[THETA_CHARGE])
      calc_charge_correlations_b();

   if (Ncomp == 1 && Ipot_ff_n==2 && Iliq_vap != -2) P_over_po=Betap/p_coex;
   if (Proc==0 && print_flag && Iwrite!=NO_SCREEN){
      print_thermo(output_file1,betap_hs,betamu_hs);
      printf("-------------------------------------------------------------------------------\n");
   }
 
 }
 else{  /* POLYMER properties */
   
   /* the stencil has already been calculated at this point. attractions
      need to be added to c(r) earlier in the code */
   if (Proc==0 && print_flag && Iwrite==VERBOSE) printf("\n%s: No polymer thermo precalculations\n",yo);

   if (Ipot_ff_c == COULOMB && Sten_Type[THETA_CHARGE])
      calc_charge_correlations_b();

 } /* end of   if (!Sten_Type[POLYMER_CR])   */
   
}
/****************************************************************************/
/* pot_parameters: calculate the cross terms (sigmaij,epsilonij,cutoffij)
   for this potential */
void pot_parameters(char *output_file1)
{
 int i,j,iw,jw, printproc=FALSE;
 double cut=0.0;
 FILE *fp2=NULL;
 if (Proc==0) printproc = TRUE;
 else printproc=FALSE;
 if (printproc) fp2 = fopen(output_file1,"a+");

 if (printproc) fprintf(fp2,"****** SUMMARY OF POTENTIAL PARAMETERS ************\n");

 for (i=0; i<Ncomp; i++){
     if (Ipot_ff_n==IDEAL_GAS) Sigma_ff[i][i]=0.0;
     if (Type_poly ==NONE) Bond_ff[i][i]=0.0;

     for (j=0; j<Ncomp; j++){

        if (i != j){
           Sigma_ff[i][j] = 0.5*(Sigma_ff[i][i]+Sigma_ff[j][j]);
           Eps_ff[i][j] = sqrt(Eps_ff[i][i]*Eps_ff[j][j]);
           Cut_ff[i][j] = 0.5*(Cut_ff[i][i]+Cut_ff[j][j]);
           Bond_ff[i][j] = 0.5*(Bond_ff[i][i]+Bond_ff[j][j]);
        }

        if (printproc) {
           fprintf(fp2,"\ti:%d  j:%d",i,j);
           fprintf(fp2,"   Sigma_ff: %9.6f  Cut_ff: %9.6f  Eps_ff: %9.6f  Bond_ff: %9.6f \n",
                   Sigma_ff[i][j],Cut_ff[i][j],Eps_ff[i][j],Bond_ff[i][j]);
        }
     }
 } 

 for (iw=0; iw<Nwall_type; iw++){
   for (jw=0; jw<Nwall_type; jw++){
        Sigma_ww[iw][jw] = 0.5*(Sigma_w[iw]+Sigma_w[jw]);
        Eps_ww[iw][jw] = sqrt(Eps_w[iw]*Eps_w[jw]);
        if (iw != jw){
           Cut_ww[iw][jw] = 0.5*(Cut_ww[iw][iw]+Cut_ww[jw][jw]);
        }
        if (printproc){
          fprintf(fp2,"\tiwall_type:%d  jwall_type:%d",iw,jw);
          fprintf(fp2,"   Sigma_ww: %9.6f  Cut_ww: %9.6f  Eps_ww: %9.6f\n",
                           Sigma_ww[iw][jw],Cut_ww[iw][jw],Eps_ww[iw][jw]);
        }
   }
 }

 for (i=0; i<Ncomp; i++){

     for (iw=0; iw<Nwall_type; iw++){
        Sigma_wf[i][iw] = 0.5*(Sigma_ff[i][i] + Sigma_w[iw]);
        Eps_wf[i][iw] = sqrt(Eps_ff[i][i]*Eps_w[iw]);
        Cut_wf[i][iw] = 0.5*(Cut_ff[i][i]+Cut_ww[iw][iw]);

        if (printproc) {
          fprintf(fp2,"\ti:%d  iwall_type:%d",i,iw);
          fprintf(fp2,"   Sigma_wf: %9.6f  Cut_wf: %9.6f  Eps_wf: %9.6f\n",
           Sigma_wf[i][iw],Cut_wf[i][iw],Eps_wf[i][iw]);
        }
     }
 } 


 if (printproc) {
     fprintf(fp2,"***************************************************\n");
     fclose(fp2);
 }
}
/*******************************************************************************/
/* coexistence: This routine finds the coexisting densities of a pure fluid
                given a temperature */
double coexistence()
{
  double resid[2],Jac[2][2],Jac_inv[2][2];
  double betamu_hs_l[1],betamu_hs_g[1];
  double betamu_att_l[1],betamu_att_g[1];
  double rho_coex_old[2],p_coex;
  int i,j,niter;

  /* for one component only */

  /* Initial guesses for the densitites */
   Rho_coex[0] = 0.7;
   Rho_coex[1] = 0.01;

   /* calculate residuals at these densities */

   resid[0] = ( Rho_coex[0]+calc_hs_properties(betamu_hs_l,&(Rho_coex[0])) 
                      + calc_att_properties(betamu_att_l,&(Rho_coex[0])) ) -
              ( Rho_coex[1]+calc_hs_properties(betamu_hs_g,&(Rho_coex[1])) 
                      + calc_att_properties(betamu_att_g,&(Rho_coex[1])) );

   resid[1] = ( log(Rho_coex[0])+betamu_hs_l[0] + betamu_att_l[0]) - 
              ( log(Rho_coex[1])+betamu_hs_g[0] + betamu_att_g[0]);

   niter = 0;
if (niter==0) printf("Rho_coex=%9.6f  %9.6f\n",Rho_coex[0],Rho_coex[1]);
if (niter==0) printf("betamu_hs_l=%9.6f  betamu_hs_g=%9.6f  betamu_att_l=%9.6f  betamu_att_g=%9.6f\n",
                      betamu_hs_l[0],betamu_hs_g[0],betamu_att_l[0],betamu_att_g[0]);
if (niter==0) printf("resid= %9.6f  %9.6f\n", resid[0],resid[1]);


   /* perform Newton iterations to find the root */

   while (fabs(resid[0]) > 1e-10 || fabs(resid[1]) > 1e-10){

        rho_coex_old[0] = Rho_coex[0];
        rho_coex_old[1] = Rho_coex[1];

        Jac[0][0] = dp_drho_hs (&(Rho_coex[0])) + dp_drho_att(&(Rho_coex[0]));
        Jac[0][1] = dp_drho_hs (&(Rho_coex[1])) + dp_drho_att(&(Rho_coex[1]));
        Jac[1][0] = dmu_drho_hs(&(Rho_coex[0])) + dmu_drho_att(&(Rho_coex[0]));
        Jac[1][1] = dmu_drho_hs(&(Rho_coex[1])) + dmu_drho_att(&(Rho_coex[1]));

if (niter==0) printf("Jac: %9.6f  %9.6f  %9.6f  %9.6f\n",Jac[0][0],Jac[0][1],Jac[1][0],Jac[1][1]);


        /* matrix inversion is simple for this two by two */

        Jac_inv[1][0] = 1.0/(Jac[0][1]-Jac[0][0]*Jac[1][1]/Jac[1][0]);
        Jac_inv[1][1] = 1.0/(Jac[1][1]-Jac[1][0]*Jac[0][1]/Jac[0][0]);
        Jac_inv[0][0] = -(Jac[1][1]/Jac[1][0])*Jac_inv[1][0];
        Jac_inv[0][1] = -(Jac[0][1]/Jac[0][0])*Jac_inv[1][1];
if (niter==0) printf("Jac: %9.6f  %9.6f  %9.6f  %9.6f\n",Jac_inv[0][0],Jac_inv[0][1],Jac_inv[1][0],Jac_inv[1][1]);

        for (i=0; i< 2; i++) {
          for (j=0; j< 2; j++) {
             Rho_coex[i] += Jac_inv[i][j]*resid[j];
          }
        }

        if (Rho_coex[1] < 0.0) Rho_coex[1] = 0.0001 ;
        if (Rho_coex[0] > 1.0) Rho_coex[0] = rho_coex_old[0] ;


        resid[0] = ( Rho_coex[0]+calc_hs_properties(betamu_hs_l,&(Rho_coex[0])) 
                      + calc_att_properties(betamu_att_l,&(Rho_coex[0])) ) -
                   ( Rho_coex[1]+calc_hs_properties(betamu_hs_g,&(Rho_coex[1])) 
                      + calc_att_properties(betamu_att_g,&(Rho_coex[1])) );

        resid[1] = ( log(Rho_coex[0])+betamu_hs_l[0] + betamu_att_l[0]) - 
                   ( log(Rho_coex[1])+betamu_hs_g[0] + betamu_att_g[0]);
        niter++;
        if (niter > 100) {
           printf("problems locating coexistence....more than 100 iterations");
           exit(-1);
        }
   }
   p_coex = calc_hs_properties(betamu_hs_l,&(Rho_coex[0])) 
                      + calc_att_properties(betamu_att_l,&(Rho_coex[0])) ;
   return(p_coex);
}
/******************************************************************************
print_thermo_properties:  do all the file printing here.*/
void print_thermo(char *output_file1, double betap_hs, 
                      double *betamu_hs)
{

  FILE *fp2;
  int icomp,jcomp,ipol,iseg;
   double betap_id;

  fp2 = fopen(output_file1,"a+");
  fprintf(fp2,"\n!!!!!!!!!!!!! output from dft_thermo.c !!!!!!!!!!!!!!!!!!\n");

      fprintf(fp2,"\n     ********** component sphere based contributions ***********\n");
      fprintf(fp2,"           **** ideal gas contributions ****\n");
      betap_id=0.0;
      for (icomp=0;icomp<Ncomp;icomp++)
         betap_id += Rho_b[icomp];

      fprintf(fp2,"\tideal gas pressure : betap_id  = %e\n",betap_id);
      for (icomp=0;icomp<Ncomp;icomp++)
      fprintf(fp2,"\tideal component chemical potentials:betap_mu_id[%d]=%e\n",
                                                 icomp,log(Rho_b[icomp]));
   if (Ipot_ff_n != IDEAL_GAS) {
      fprintf(fp2,"\n           **** hard sphere contributions ****\n");
      fprintf(fp2,"\thard sphere pressure : betap_hs  = %e\n",betap_hs);
      fprintf(fp2,"\texcess hard sphere chemical potentials  (units of kT)\n" );
      for (icomp=0; icomp<Ncomp; ++icomp) 
           fprintf(fp2,"\t\tBetamu_hs_ex[%d] = %e\n",icomp, Betamu_hs_ex[icomp]);
   }

   if (Ipot_ff_n == LJ12_6){
      fprintf(fp2,"\n           **** attractive contributions ****\n");
      for (icomp=0; icomp<Ncomp; ++icomp) 
        fprintf(fp2,"\tBetamu_att[%d] = %e\n\n",icomp, Betamu_att[icomp]); 

      for (icomp=0; icomp<Ncomp; ++icomp) 
       for (jcomp=0; jcomp<Ncomp; ++jcomp) 
        fprintf(fp2,"\tAvdw[%d][%d] = %e\n\n",icomp,jcomp,Avdw[icomp][jcomp]); 
   }
   if (Ipot_ff_c > 0 && Sten_Type[THETA_CHARGE]){
      fprintf(fp2,"\n           **** charge contributions ****\n");
      for (icomp=0; icomp<Ncomp; ++icomp) 
        fprintf(fp2,"\tDeltac_b[%d] = %e\n\n",icomp, Deltac_b[icomp]);                  
   }
      fprintf(fp2,"\n           ******* TOTAL OF ABOVE !!!*******\n");
  fprintf(fp2,"\tpressure : psigma^3/kT : \t%9.6f\n",Betap);
  fprintf(fp2,"\tchemical potentials : mu/kT : \n");
  for (icomp=0; icomp<Ncomp; ++icomp) 
    fprintf(fp2,"\t\t\t Betamu[%d] = %9.6f\n",icomp, Betamu[icomp]);
      fprintf(fp2,"           *********************************\n");

   if (Type_poly==WTC){
      fprintf(fp2,"\n     ***********************************************************\n");
      fprintf(fp2,"     *** NOTE THE ABOVE IS NOT APPROPRIATE FOR BONDED SYSTEMS ***\n");
      fprintf(fp2,"     *** INSTEAD WE REQUIRE SEGMENT BASED CHEMICAL POTENTIALS ***\n");
      fprintf(fp2,"     *** THAT INCLUDE THE CONTRIBUTIONS OF THE BONDS !!       ***\n");
      fprintf(fp2,"     ***********************************************************\n");
      for (iseg=0; iseg<Nseg_tot;iseg++){
          fprintf(fp2,"\t\t Betamu_wtc[iseg=%d] = %9.6f\n", iseg, Betamu_wtc[iseg]);
      }
      fprintf(fp2,"\n    ******* TOTAL SEGMENT BASED CHEMICAL POTENTIALS !!!*******\n");
      for (iseg=0; iseg<Nseg_tot;iseg++){
          fprintf(fp2,"\t\t Betamu_seg[iseg=%d] = %9.6f\n", iseg, Betamu_seg[iseg]);
      }
      fprintf(fp2,"           *********************************\n");
   }
   if (Type_poly==WTC){
      fprintf(fp2,"correcting for overlapping particles in WTC bonded systems\n");
      fprintf(fp2,"\t\t cavity correlation corrections\n");
      for (icomp=0;icomp<Ncomp;icomp++)
         for (jcomp=0;jcomp<Ncomp;jcomp++)
           fprintf(fp2,"icomp=%d jcomp=%d Fac_overlap=%9.6f\n",icomp,jcomp,Fac_overlap[icomp][jcomp]);
   }
   fprintf(fp2,"\t\t checking hard sphere corrections\n");
   for (icomp=0;icomp<Ncomp;icomp++)
        fprintf(fp2,"icomp=%d Fac_overlap_hs=%9.6f\n",icomp,Fac_overlap_hs[icomp]);
   fprintf(fp2,"           *********************************\n");

   fprintf(fp2,"\n!!!!!!!!!!!!! end of dft_thermo.c output !!!!!!!!!!!!!!!!!!\n");
   fclose(fp2);

   return;
}
/***********************end of thermo file **********************************/
