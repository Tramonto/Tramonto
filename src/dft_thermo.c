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
double calc_ideal_gas(double *,double *);
double calc_hs_properties(double *,double *);
double calc_att_properties(double *, double *);
void calc_poly_TC_properties(double *,double *);
void   calc_charge_correlations_b();
double int_stencil_bulk(int,int,int);
double   coexistence();
double dp_drho_hs(double *);
double dp_drho_att(double *);
void print_thermo(char *,double,double *);
void compute_bulk_nonlocal_properties(char *);
void sum_rhobar(double,int, double *);

void  thermodynamics( char *output_file1, int print_flag)
{
   int icomp,jcomp,kcomp,i,ipol,iseg,ibond,jseg,kseg;
   double betap_hs=0, betamu_hs[NCOMP_MAX],betap_att,p_coex=0.0,betap_TC=0.0;
   double y,dydxi2,dydxi3,sig2,sig3,sum;
   char *yo = "thermodynamics";

 if (!Sten_Type[POLYMER_CR]){
   if (Proc==0 && print_flag &&Iwrite!=NO_SCREEN){
          printf("\n-------------------------------------------------------------------------------\n");
          printf("%s: Doing Thermo precalculations\n",yo);
   }
   if (Ipot_ff_n == LJ12_6) {
      Avdw = (double **) array_alloc(2, Ncomp, Ncomp, sizeof(double));
   }
   Betamu_att = (double *) array_alloc(1, Ncomp, sizeof(double));
   for (i=0;i<Ncomp;i++) Betamu_att[i]=0.0;
   betap_att=0.0;

 /*  if (Type_poly_TC){
         for (icomp=0; icomp<Ncomp; icomp++)
             for (i=0;i<Nbonds;i++)  Betamu_ex_bondTC[icomp][i]=0.0;
   }*/
   if (Type_poly_TC){
      for (iseg=0;iseg<Nseg_tot;iseg++) Betamu_wtc[iseg]=0.0;
   }

   /* Find bulk coexistence  for a very special case of only one atomistic component with
      mean field attractions turned on !! */

   if (Ncomp == 1 && Ipot_ff_n==2 && Iliq_vap != -2){ 
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
   }

  /* for lack of a better place, precalc these quantities here */
 
  if (Ipot_ff_n != IDEAL_GAS) {

     Inv_4pi = 1.0 / (4.0*PI);
     for (icomp=0; icomp<Ncomp; icomp++) {
     if (Sigma_ff[icomp][icomp] < 1.e-6){
        printf("ERROR : icomp=%d Sigma_ff[%d][%d]=%9.6f --- check particle sizes!\n",
                    icomp,icomp,icomp,Sigma_ff[icomp][icomp]);
         exit(-1);
     }
        Inv_rad[icomp] = 2.0 / Sigma_ff[icomp][icomp];
        Inv_4pir[icomp] = Inv_4pi * Inv_rad[icomp];
        Inv_4pirsq[icomp] = Inv_4pir[icomp] * Inv_rad[icomp];
     }
  }

  if (Type_poly_TC){
     for (iseg=0;iseg<Nseg_tot;iseg++){
         Rho_seg_b[iseg]=Rho_b[Unk2Comp[iseg]]/Nmer_t_total[Unk2Comp[iseg]];
         if (Lsteady_state){
            Rho_seg_LBB[iseg]=Rho_b_LBB[Unk2Comp[iseg]]/Nmer_t_total[Unk2Comp[iseg]];
            Rho_seg_RTF[iseg]=Rho_b_RTF[Unk2Comp[iseg]]/Nmer_t_total[Unk2Comp[iseg]];
         }
     }
  }

   compute_bulk_nonlocal_properties(output_file1);

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
  
   if (Type_func==ROSENFELD || Type_func==ROSENFELD2){
      if (Lsteady_state){
         betap_hs = calc_hs_properties(betamu_hs,Rho_b_LBB);
         for (icomp=0; icomp<Ncomp; icomp++) Betamu_LBB[icomp] += betamu_hs[icomp];    

         betap_hs = calc_hs_properties(betamu_hs,Rho_b_RTF);
         for (icomp=0; icomp<Ncomp; icomp++) Betamu_RTF[icomp] += betamu_hs[icomp];
      }
      else{
         betap_hs = calc_hs_properties(betamu_hs,Rho_b);
         Betap += betap_hs;
         for (icomp=0; icomp<Ncomp; icomp++) Betamu[icomp] += betamu_hs[icomp];
     }
   }

     /* now add in the mean field attraction contributions */
   if (Type_attr == LJ_WCA_CS){
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

     /* now add in an applied electric field if one exists */
   if (Lsteady_state && Type_coul != NONE){
      for (icomp=0; icomp<Ncomp; icomp++){
         Betamu_LBB[icomp] += Charge_f[icomp]*Elec_pot_LBB;
         Betamu_RTF[icomp] += Charge_f[icomp]*Elec_pot_RTF;
      }
   }

    /* now add in the WTC polymer contributions */
   if (Type_poly_TC){
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
/*                y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);*/
                y = y_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                  Betamu_seg[iseg] += 0.5*(1.0-log(y)-log(Rho_seg_b[jseg])  
                                      - Rho_seg_b[jseg]/Rho_seg_b[iseg]);
                  Betamu_wtc[iseg] += 0.5*(1.0-log(y)-log(Rho_seg_b[jseg])     
                                       -Rho_seg_b[jseg]/Rho_seg_b[iseg]);
            }
         }
         /* now add in the term that involves _all_ the bond pairs in the system for each segment via that cavity
            correlation function.  Note that this term is nonzero even for bond pairs on different polymer
            components than the one where the iseg segment is found ! */
         for (kseg=0; kseg<Nseg_tot;kseg++){
 sum=0.0;
           kcomp = Unk2Comp[kseg];
           sig2=Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp];
           sig3=Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp];
           for (iseg=0; iseg<Nseg_tot;iseg++){
              icomp=Unk2Comp[iseg];
              for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
                jseg=Bonds_SegAll[iseg][ibond];
                jcomp=Unk2Comp[jseg];
/*                y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                dydxi2 = dy_dxi2_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                dydxi3 = dy_dxi3_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);*/
                y = y_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                dydxi2 = dy_dxi2_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                dydxi3 = dy_dxi3_cav(Bond_ff[icomp][icomp],Bond_ff[jcomp][jcomp],Xi_cav_b[2],Xi_cav_b[3]);
                Betamu_seg[kseg] -= (PI/12.0)*(Rho_seg_b[iseg]/y)*(dydxi2*sig2+dydxi3*sig3);
                Betamu_wtc[kseg] -= (PI/12.0)*(Rho_seg_b[iseg]/y)*(dydxi2*sig2+dydxi3*sig3);
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
   if (Proc==0 && print_flag && Iwrite!=NO_SCREEN) printf("\n%s: No polymer thermo precalculations\n",yo);

   if (Ipot_ff_c == COULOMB && Sten_Type[THETA_CHARGE])
      calc_charge_correlations_b();

 } /* end of   if (!Sten_Type[POLYMER_CR])   */
   
}
/****************************************************************************/
   double y_cav(double sigma_1,double sigma_2,double xi_2, double xi_3)
{           
   double one_m_xi3_sq,one_m_xi3_cb,sig_sum,sig_m,y;
                   
   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3); 
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3);
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;
              
   y = 1./(1.-xi_3) +(3.*sig_m/sig_sum)*(xi_2/one_m_xi3_sq)+
       2.*(sig_m/sig_sum)*(sig_m/sig_sum)*(xi_2*xi_2/one_m_xi3_cb);
                 
   return y;     
}
/********************************************************************************************/
   double dy_dxi2_cav(double sigma_1,double sigma_2,double xi_2, double xi_3)
{                
   double one_m_xi3_sq,one_m_xi3_cb,sig_sum,sig_m,dy_dxi2;

   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3); 
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3); 
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;
                 
  dy_dxi2 = (3.*sig_m/sig_sum)*(1.0/one_m_xi3_sq)+
       4.*(sig_m/sig_sum)*(sig_m/sig_sum)*(xi_2/one_m_xi3_cb);
                 
   return dy_dxi2;
}
/********************************************************************************************/
   double dy_dxi3_cav(double sigma_1,double sigma_2,double xi_2, double xi_3)
{
   double one_m_xi3_sq,one_m_xi3_cb,sig_sum,sig_m,dy_dxi3,one_m_xi3_4th;

   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3);
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3);
   one_m_xi3_4th= (1.-xi_3)*one_m_xi3_cb;
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;

  dy_dxi3 = 1./one_m_xi3_sq +(6.*sig_m/sig_sum)*(xi_2/one_m_xi3_cb)+
       6.*(sig_m/sig_sum)*(sig_m/sig_sum)*(xi_2*xi_2/one_m_xi3_4th);

   return dy_dxi3;
}
/******************************************************************************/
/* calc_ideal_gas: This routine computes ideal gas properties at a density of interest */
double calc_ideal_gas(double *rho,double *betamu)
{
   double betap=0.0;
   int i;
   
   for (i=0;i<Ncomp;i++){
       betamu[i] = log(rho[i]);
                  /*           - 3.0*log(Sigma_ff[icomp][icomp]) -
                               1.5*log(Mass[icomp]*Temp); */
       betap+=rho[i];
   }
   return(betap);
}
/********************************************************************************
calc_hs_properties:  This routine calculates the pressure and excess chemical 
                     potential for hard spheres at the density of interest 
                     using the PY equations.            */

double calc_hs_properties(double *betamu_hs,double *rho)
{
   int icomp;
   double pi6, hs_diam_cubed, xsi0, xsi1, xsi2, xsi3, y1, y2, y3,
          betap_hs;

   xsi0 = 0.0;
   xsi1 = 0.0;
   xsi2 = 0.0;
   xsi3 = 0.0;
   pi6 = PI/6.0;                 /* shorthand  for pi/6                 */

   /*  Determine the effective hard sphere diameters 
    For now we will set these to unity, but in the future we
    can define a temperature dependent diameter. Doing so
    offers a way to provide a better mean field equation of
    state. In essence, for  an attractive (e.g LJ) fluid (mixture)
    we can use the effective hard sphere diameter to off set the 
    shortcomings of the PY + mean field approximation. See the
    work by Telo da Gama et al.. 
   */

   /* calculate the hard sphere diamtere ... this can be 
      turned into a T-dependent diameter */

   for (icomp=0; icomp<Ncomp; ++icomp) Hs_diam[icomp] = 1.0;  

   /*  calculate the constants xsi and introduce some shorthand */

   for (icomp=0; icomp<Ncomp; ++icomp) {
      hs_diam_cubed = POW_DOUBLE_INT(Hs_diam[icomp],3);
      xsi0 +=pi6 * rho[icomp] * hs_diam_cubed;
      xsi1 +=pi6 * rho[icomp] * hs_diam_cubed * Sigma_ff[icomp][icomp];
      xsi2 +=pi6 * rho[icomp] * hs_diam_cubed * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
      xsi3 +=pi6 * rho[icomp] * hs_diam_cubed * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);
   }
   y1 = 1.0 - xsi3;
   y2 = y1 * y1;
   y3 = y1 * y1 * y1;

   /* the hard sphere pressure in units of kT and Sigma_ff[1]^3 */

   betap_hs = (1.0/pi6) * (xsi0/y1 + 3.0 * xsi1 * xsi2/y2 +
                                     3.0 * POW_DOUBLE_INT(xsi2,3)/y3  );

   /* the excess hard sphere chemical potential in units of kT */

   for (icomp=0; icomp<Ncomp; ++icomp) 
      Betamu_hs_ex[icomp] = - log(y1) + 
              pi6 * betap_hs * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3) +
              3.0 * xsi2 * Sigma_ff[icomp][icomp]/y1 +
              3.0 * xsi1 * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2)/y1 +
              4.5 * POW_DOUBLE_INT((xsi2 * Sigma_ff[icomp][icomp]),2)/y2 ;

   /* 
    * add the ideal gas term to give
    * the hard sphere chemical potential in units of kT 
    */
   for (icomp=0; icomp<Ncomp; ++icomp) {
      betamu_hs[icomp]  =  Betamu_hs_ex[icomp]; 
   }

   return (betap_hs);
}
/*************************************************************
calc_att_properties: In this routine calculate the strict mean field
                     attraction contribution to the pressure and
                     chemical potential */
double calc_att_properties(double *betamu_att, double *rho)
{
  int icomp,jcomp;
  double betap_att;

  betap_att = 0.0; 
  for (icomp=0; icomp<Ncomp; icomp++) betamu_att[icomp] = 0.0;

  for (icomp=0; icomp<Ncomp; icomp++) {
     for (jcomp=0; jcomp<Ncomp;jcomp++){
       Avdw[icomp][jcomp] = int_stencil_bulk(U_ATTRACT,icomp,jcomp);
       betamu_att[icomp] += rho[jcomp]*Avdw[icomp][jcomp];
       betap_att += 0.5*Avdw[icomp][jcomp]*rho[icomp]*rho[jcomp];
     }
  }
  return(betap_att);
}
/********************************************************************
calc_charge_correlations_b: Here we set up the bulk cross correlations
      between the hard sphere and coulomb parts of the potential*/
void calc_charge_correlations_b()
{
   int icomp,jcomp;

   Deltac_b = (double *) array_alloc (1, Ncomp, sizeof(double));
   for (icomp=0; icomp<Ncomp; icomp++) Deltac_b[icomp] = 0.0;

   for (icomp=0; icomp<Ncomp; icomp++) 
      for (jcomp=0; jcomp<Ncomp;jcomp++)
          Deltac_b[icomp]+= Rho_b[jcomp]*
                            int_stencil_bulk(THETA_CHARGE,icomp,jcomp);
   return;
}
/***********************************************************************
int_stencil_bulk: this routine sums the appropriate stencil to get
                  the bulk contributions to various terms in the E-L
                  equation. */
double int_stencil_bulk(int sten_type,int icomp,int jcomp)
{
  int izone, isten;
  double sum, weight, *sten_weight;
  struct Stencil_Struct *sten;

  sum = 0.0;
  izone = 0;
  sten = &(Stencil[sten_type][izone][icomp+Ncomp*jcomp]);
  sten_weight = sten->Weight;

  for (isten = 0; isten<sten->Length; isten++){
     weight = sten_weight[isten];
     sum += weight;
  }
  return(sum);
}
/*******************************************************************************/
/* deltaC_MSA:  given r12, calculate the attractive part of a cut and
           shifted 12-6 LJ potential. */

double deltaC_MSA(double r,int i, int j)
{
  double deltac,B,kappa,kappa_sq;
  int icomp;

  kappa_sq = 0.0;
  for(icomp = 0; icomp<Ncomp; icomp++)
     kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*
                  Charge_f[icomp]*Charge_f[icomp];
  kappa = sqrt(kappa_sq);
  B = (kappa + 1.0 - sqrt(1.0+2.0*kappa))/kappa;

/*  printf("\t r: %9.6f icomp: %d  jcomp: %d  kappa: %9.6f  B: %9.6f  Sigma_ff: %9.6f\n",
          r,i,j,kappa,B,Sigma_ff[i][j]);*/

  if (r == 0.0) printf("trouble with deltaC term .... r=0");
  if (r <= Sigma_ff[i][j] && r>0) {

     deltac = -Charge_f[i]*Charge_f[j]/Temp_elec*
              (  2*B/Sigma_ff[i][j] - 1.0/r
               - POW_DOUBLE_INT(B/Sigma_ff[i][j],2)*r );
  }
  else deltac = 0.0;

  return deltac;
}
/*******************************************************************************/
/* deltaC_MSA_int:  given range of integrattion, r, calculate the definite
           integral of deltac_MSA over all space */

double deltaC_MSA_int(double r,int i, int j)
{
  double deltac_int,B,kappa,kappa_sq;
  int icomp;

  kappa_sq = 0.0;
  for(icomp = 0; icomp<Ncomp; icomp++)
     kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*
                  Charge_f[icomp]*Charge_f[icomp];
  kappa = sqrt(kappa_sq);
  B = (kappa + 1.0 - sqrt(1.0+2.0*kappa))/kappa;

  deltac_int = -(4*PI*Charge_f[i]*Charge_f[j]/Temp_elec)*
                r*r*
               (  2*B*r/(3.0*Sigma_ff[i][j]) - 0.5
               - 0.25*POW_DOUBLE_INT(B/Sigma_ff[i][j],2)*r*r );

  return deltac_int;
}
/*******************************************************************************/
/* compute_bulk_nonlocal_properties: compute some additional bulk properties we
   need to carry around the calculation. These are based on the input densities
   and particle sizes. */
void compute_bulk_nonlocal_properties(char *output_file1)
{
  int i,loc_inode,loc_i,inode_box,inode,ijk[3],icomp,idim,iunk,printproc;
  int ibond,iseg,jseg,pol_number,type_jseg,nloop,iloop;
  double vol,area,x_dist;
  FILE *fp2=NULL;
  if (Proc==0) printproc = TRUE;
  else printproc=FALSE;
  if (printproc) fp2 = fopen(output_file1,"a+");

   /* compute bulk nonlocal densities needed for Rosenfeld terms */
  if (Type_func != NONE){   
     for (iunk=0; iunk<Nrho_bar; iunk++){
       Rhobar_b[iunk] = 0.0;
       Rhobar_b_LBB[iunk] = 0.0;
       Rhobar_b_RTF[iunk] = 0.0;
     }

    if (Type_poly_TC) nloop=Nseg_tot;
    else             nloop=Ncomp;

     for (iloop=0; iloop<nloop; iloop++){
       if (Type_poly_TC) icomp=Unk2Comp[iloop];
       else             icomp=iloop;

       if (Type_poly_TC){
           sum_rhobar(Rho_seg_b[iloop],icomp,Rhobar_b);
       }
       else{
          if (Lsteady_state){
                 sum_rhobar(Rho_b_LBB[icomp],icomp,Rhobar_b_LBB);
                 sum_rhobar(Rho_b_RTF[icomp],icomp,Rhobar_b_RTF);
          }
          else if (Nwall ==0 && Iliq_vap==3){
                 sum_rhobar(Rho_coex[1],icomp,Rhobar_b_LBB);
                 sum_rhobar(Rho_coex[0],icomp,Rhobar_b_RTF);
          }
          else sum_rhobar(Rho_b[icomp],icomp,Rhobar_b);
   
       }
     }
     if (printproc){
        fprintf(fp2,"Rhobar_bulk, LBB, and RTF variables for Rosenfeld HS functionals:\n");
        fprintf(fp2,"Note that vector terms are strictly zero in the bulk!\n");
        fprintf(fp2,"\t i  Rhobar_b[i]  Rhobar_b_LBB[i]  Rhobar_b_RTF[i]\n");
        for (i=0;i<4;i++) fprintf(fp2,"\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
                 Rhobar_b[i], Rhobar_b_LBB[i], Rhobar_b_RTF[i]);
        if (Iwrite==VERBOSE){
           printf("Rhobar_bulk, LBB, and RTF variables for Rosenfeld HS functionals:\n");
           printf("Note that vector terms are strictly zero in the bulk!\n");
           printf("\t i  Rhobar_b[i]  Rhobar_b_LBB[i]  Rhobar_b_RTF[i]\n");
           for (i=0;i<4;i++) printf("\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
                 Rhobar_b[i], Rhobar_b_LBB[i], Rhobar_b_RTF[i]);
        }
     }
  }

  /* compute bulk nonlocal densities needed for Wertheim-Tripathi-Chapman functionals */
  if (Type_poly_TC){  
     for (i=0; i<4; i++){
        Xi_cav_b[i]=0.0;
        Xi_cav_LBB[i]=0.0;
        Xi_cav_RTF[i]=0.0;
     }
     for (icomp=0;icomp<Ncomp;icomp++){
        for (i=0;i<4;i++){
           Xi_cav_b[i]+=(PI/6.0)*Rho_b[icomp]*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],i);
           if (Lsteady_state){
             Xi_cav_LBB[i]+=(PI/6.0)*Rho_b_LBB[icomp]*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],i);
             Xi_cav_RTF[i]+=(PI/6.0)*Rho_b_RTF[icomp]*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],i);
           }
        }
     }
     if (printproc){
        fprintf(fp2,"Xi_cavity_bulk, LBB, and RTF variables for WTC polymer run are:\n");
        fprintf(fp2,"\t i  Xi_cav_b[i]  Xi_cav_LBB[i]  Xi_cav_RTF[i]\n");
        for (i=0;i<4;i++) fprintf(fp2,"\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
                 Xi_cav_b[i], Xi_cav_LBB[i], Xi_cav_RTF[i]);
        if (Iwrite==VERBOSE){
           printf("Xi_cav_bulk, LBB, and RTF variables for WTC polymer run are:\n");
           printf("\t i  Xi_cav_b[i]  Xi_cav_LBB[i]  Xi_cav_RTF[i]\n");
            for (i=0;i<4;i++) printf("\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
                 Xi_cav_b[i], Xi_cav_LBB[i], Xi_cav_RTF[i]);
        }
     }

     for (ibond=0; ibond<Nbonds; ibond++){
       BondWTC_b[ibond]=NO_BOND_PAIR;
       BondWTC_LBB[ibond]=NO_BOND_PAIR;
       BondWTC_RTF[ibond]=NO_BOND_PAIR;
    }

    for (ibond=0; ibond<Nbonds; ibond++){
        iseg=Unk_to_Seg[ibond];
        pol_number=Unk_to_Poly[ibond];
        jseg=Bonds[pol_number][iseg][Unk_to_Bond[ibond]];
        type_jseg=Type_mer[pol_number][jseg];
        BondWTC_b[ibond]=Rho_seg_b[jseg];
        if (Lsteady_state){
           BondWTC_LBB[ibond]=Rho_seg_LBB[jseg];
           BondWTC_RTF[ibond]=Rho_seg_RTF[jseg];
        }
    }
    if (printproc){
        fprintf(fp2,"BondWTC_bulk, LBB, and RTF variables for WTC polymer run are:\n");
        fprintf(fp2,"\t i  BondWTC_b[i]  BondWTC_LBB[i]  BondWTC_RTF[i]\n");
        for (i=0;i<Nbonds;i++) fprintf(fp2,"\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
                 BondWTC_b[i], BondWTC_LBB[i], BondWTC_RTF[i]);
        if (Iwrite==VERBOSE){
           printf("BondWTC_bulk, LBB, and RTF variables for WTC polymer run are:\n");
           printf("\t i  BondWTC_b[i]  BondWTC_LBB[i]  BondWTC_RTF[i]\n");
           for (i=0;i<Nbonds;i++) printf("\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
                 BondWTC_b[i], BondWTC_LBB[i], BondWTC_RTF[i]);
        }
       fclose(fp2);
     }

  } /* end of Type_poly_TC rhobars (bulk)*/


}
/*********************************************************************************************/
void sum_rhobar(double rho,int icomp, double *rhobar)
{
   double vol,area;
   int idim;

   vol = PI*Sigma_ff[icomp][icomp]*
           Sigma_ff[icomp][icomp]*Sigma_ff[icomp][icomp]/6.0;
   area = PI*Sigma_ff[icomp][icomp]*Sigma_ff[icomp][icomp];

   rhobar[0] += vol*rho;
   rhobar[1] += area*rho;
   rhobar[2] += area*rho*Inv_4pir[icomp];
   rhobar[3] += area*rho*Inv_4pirsq[icomp];
   for (idim=0; idim<Ndim; idim++){
       rhobar[4+idim] = 0.0;
       rhobar[4+Ndim+idim] = 0.0;
   }
   return;
}
/*******************************************************************************/
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
     if (Type_poly <0 && !Type_poly_TC) Bond_ff[i][i]=0.0;

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
/*************************************************************************
dp_drho_hs: the derivative of the hard sphere pressure with respect to rho*/
double dp_drho_hs(double *rho)
{
   int icomp;
   double Hs_diam[NCOMP_MAX],hs_diam_cubed,pi6;
   double xsi0,xsi1,xsi2,xsi3,y1,y2,y3,dy1_drho,dy2_drho,dy3_drho;
   double dp_drho;

   xsi0 = xsi1 = xsi2 = xsi3 = 0.0;
   pi6 = PI/6.0;                 /* shorthand  for pi/6                 */

   for (icomp=0; icomp<Ncomp; ++icomp) Hs_diam[icomp] = 1.0;  

   /*  calculate the constants xsi and introduce some shorthand */

   for (icomp=0; icomp<Ncomp; ++icomp) {
      hs_diam_cubed = POW_DOUBLE_INT(Hs_diam[icomp],3);
      xsi0 +=pi6 * rho[icomp] * hs_diam_cubed;
      xsi1 +=pi6 * rho[icomp] * hs_diam_cubed * Sigma_ff[icomp][icomp];
      xsi2 +=pi6 * rho[icomp] * hs_diam_cubed 
                 * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
      xsi3 +=pi6 * rho[icomp] * hs_diam_cubed 
                 * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);
   }
   y1 = 1.0 - xsi3;
   y2 = y1 * y1;
   y3 = y1 * y1 * y1;

   dy1_drho = xsi3/rho[0];
   dy2_drho = 2*y1*dy1_drho;
   dy3_drho = 3*y1*y1*dy1_drho;

   dp_drho = (1.0/pi6) * (  xsi0/(rho[0]*y1) - xsi0*dy1_drho/(y1*y1) +
                           (3.0*xsi1*xsi2/y2)*( 2.0/rho[0] - dy2_drho/y2) +
                           (3.0*xsi2*xsi2*xsi2/y3) *(3.0/(rho[0]) - dy3_drho/y3) );
   return (dp_drho);
}
/*************************************************************************
dp_drho_att: the derivative of the attractive part of the pressure 
             with respect to rho*/
double dp_drho_att(double *rho)
{
 double dp_drho;

 dp_drho = Avdw[0][0]*rho[0];
 return (dp_drho);
}
/*************************************************************************
dmu_drho_hs: the derivative of the hard sphere chemical potential 
             with respect to rho*/
double dmu_drho_hs(double *rho)
{
   int icomp;
   double Hs_diam[NCOMP_MAX],hs_diam_cubed,pi6;
   double xsi0,xsi1,xsi2,xsi3,y1,y2,y3,dy1_drho,dy2_drho,dy3_drho;
   double dmu_drho;

   xsi0 = xsi1 = xsi2 = xsi3 = 0.0;
   pi6 = PI/6.0;                 /* shorthand  for pi/6                 */

   for (icomp=0; icomp<Ncomp; ++icomp) Hs_diam[icomp] = 1.0;  

   /*  calculate the constants xsi and introduce some shorthand */

   for (icomp=0; icomp<Ncomp; ++icomp) {
      hs_diam_cubed = POW_DOUBLE_INT(Hs_diam[icomp],3);
      xsi0 +=pi6 * rho[icomp] * hs_diam_cubed;
      xsi1 +=pi6 * rho[icomp] * hs_diam_cubed * Sigma_ff[icomp][icomp];
      xsi2 +=pi6 * rho[icomp] * hs_diam_cubed 
                 * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
      xsi3 +=pi6 * rho[icomp] * hs_diam_cubed 
                 * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);
   }
   y1 = 1.0 - xsi3;
   y2 = y1 * y1;
   y3 = y1 * y1 * y1;

   dy1_drho = xsi3/rho[0];
   dy2_drho = 2*y1*dy1_drho;
   dy3_drho = 3*y1*y1*dy1_drho;

   dmu_drho = - dy1_drho/y1 + 
        pi6 * dp_drho_hs(rho) * POW_DOUBLE_INT(Sigma_ff[0][0],3) +
       (3.0 * xsi2 * Sigma_ff[0][0]/y1)
                                 *( 1.0/rho[0] - dy1_drho/y1) +
       (3.0 * xsi1 * POW_DOUBLE_INT(Sigma_ff[0][0],2)/y1)
                                 *( 1.0/rho[0] - dy1_drho/y1) +
       (9.0 * POW_DOUBLE_INT((xsi2 * Sigma_ff[0][0]),2)/y2)
                                 *( 1.0/rho[0] - 0.5*dy2_drho/y2) +
       (1.0/rho[0]);

   return(dmu_drho);
}
/*************************************************************************
dmu_drho_att: the derivative of the attractive part of the chemical potential 
               with respect to rho*/
double dmu_drho_att(double *rho)
{
   double dmu_drho;

   dmu_drho = Avdw[0][0];
   return (dmu_drho);
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

   if (Type_poly_TC){
      fprintf(fp2,"\n     ***********************************************************\n");
      fprintf(fp2,"     *** NOTE THE ABOVE IS NOT APPROPRIATE FOR BONDED SYSTEMS ***\n");
      fprintf(fp2,"     *** INSTEAD WE REQUIRE SEGMENT BASED CHEMICAL POTENTIALS ***\n");
      fprintf(fp2,"     *** THAT INCLUDE THE CONTRIBUTIONS OF THE BONDS !!       ***\n");
      fprintf(fp2,"     ***********************************************************\n");
      for (iseg=0; iseg<Nseg_tot;iseg++){
          fprintf(fp2,"\t\t Betamu_seg[iseg=%d] = %9.6f\n", iseg, Betamu_seg[iseg]);
      }
   }

   fprintf(fp2,"\n!!!!!!!!!!!!! end of dft_thermo.c output !!!!!!!!!!!!!!!!!!\n");
   fclose(fp2);

   return;
}
/***********************end of thermo file **********************************/
