/*
//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/* ---------------------------------------------------------
Routines to Calculate the thermodynamic properties of hard sphere fluids.
Two routines are included for hard sphere properties -
(1) is based on the Percus-Yevick compressibilty equation of state.
the other (2) uses the same notation as the DFT Euler-Lagrange fill later 
in the code.  These should give the same results.
------------------------------------------------------------*/
#include "dft_thermo_hs.h"

/*#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

void rhobar_icomp(double,int, double *);
void dphi_drb_bulk(double *,double *);*/



/********************************************************************************
/*HS_thermo_precalc: logic to control preprocessing of needed parameters when 
                   we are using FMT functionals */
void HS_thermo_precalc(char *output_file1)
{
   int icomp;

   /* make sure overlap factors are turned off if not doing WTC functionals */
   if (Type_poly != WTC) for (icomp=0;icomp<Ncomp;icomp++) Fac_overlap_hs[icomp]=1.0;

   calc_InvR_params();

  compute_bulk_FMT_properties(output_file1);


   return;
}
/********************************************************************************
calc_InvR_params: This routine computes factors of 1/R needed in FMT functionals */
void calc_InvR_params()
{
   int icomp;

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

   return;
}
/********************************************************************************
calc_HS_diams: This routine computes diameter to use in hard sphere functionals */
void calc_HS_diams()
{
  int icomp,i;
  double del,r,sum;
  for (icomp=0;icomp<Ncomp;icomp++){
     if (Type_attr==LJ_BH_MF){ /* Barker-Henderson Treatment of hard core diameters */
        del=Sigma_ff[icomp][icomp]/1000.;
        r=sum=0.0;
        for (i=0;i<1000;i++){
          if (i==0 || i==999) sum+=0.5*integrand_BH(r,icomp)*del;
          else sum+=integrand_BH(r,icomp)*del;
          r+=del;
        }
        HS_diam[icomp]=sum;
     }
     else HS_diam[icomp]=Sigma_ff[icomp][icomp];

     if (Iwrite==VERBOSE && Proc==0) printf("BH hard cores: icomp=%d  Sigma_ff=%9.6f  eps/kT=%9.6f HS_diam=%9.6f\n",
             icomp,Sigma_ff[icomp][icomp],Eps_ff[icomp][icomp],HS_diam[icomp]);
  }
  return;
}
/********************************************************************************/
double integrand_BH(double r,int icomp)
{
  double integrand,rcut;
  rcut=1.e6;  /* set very large to eliminate the cut and shift for the BH diameters */
  if (r>.01) integrand = 1.-exp(pairPot_switch(r,Sigma_ff[icomp][icomp],Eps_ff[icomp][icomp],rcut));
  else integrand = 1.0;
  return(integrand);
}
/********************************************************************************
pressure_FMT_hs  This routine calculates the pressure contribution for 
                 hard sphere fundamental measures theory based DFTs */

double pressure_FMT_hs(double *rho)
{
   int i;
   double betap_hs,n[4+2*NDIM_MAX],rho_bar[4+2*NDIM_MAX];
 
   for (i=0;i<Nrho_bar_s;i++) rho_bar[i]=Rhobar_b[i];
   for (i=0;i<2*Ndim;i++) rho_bar[Nrho_bar_s+i]=0.0; 
   solutionVec_to_nOrdering(rho_bar,n);

   betap_hs = -phispt_switch(n);
   for (i=0;i<4;i++) {
      betap_hs += Dphi_Drhobar_b[i]*n[i];
   }
   betap_hs +=n[0];

   return (betap_hs);
}
/********************************************************************************
chempot_FMT_hs:  This routine calculates the excess chemical potential for 
                 hard sphere fundamental measures theory based DFTs */

void chempot_FMT_hs(double *rho)
{
   int icomp,i;
   double sten_sum[4],n[4+2*NDIM_MAX];
  
   for (icomp=0; icomp<Ncomp;icomp++){
      sten_sum[0]=1.;
      sten_sum[1]= HS_diam[icomp]/2.;
      sten_sum[2]=(4.*PI)*POW_DOUBLE_INT(HS_diam[icomp]/2.,2);
      sten_sum[3]=(4.*PI/3.)*POW_DOUBLE_INT(HS_diam[icomp]/2.,3);
  
      Betamu_hs_ex[icomp]=0.0; 
      for (i=0;i<4;i++){
        Betamu_hs_ex[icomp] += Dphi_Drhobar_b[i]*sten_sum[i]*Fac_overlap_hs[icomp];
      }
   }
   return;
}
/********************************************************************************
pressure_PY_hs:  This routine calculates the pressure 
                 hard spheres at the density of interest using the PY equations.            */

double pressure_PY_hs(double *rho)
{
   int icomp;
   double pi6, hs_diam_cubed, xsi0, xsi1, xsi2, xsi3, y1, y2, y3,
          betap_hs;

   xsi0=xsi1=xsi2=xsi3=0.0;
   pi6 = PI/6.0;                 /* shorthand  for pi/6                 */

   /*  calculate the constants xsi and introduce some shorthand */

   for (icomp=0; icomp<Ncomp; ++icomp) {
      hs_diam_cubed = POW_DOUBLE_INT(HS_diam[icomp],3);
      xsi0 += Fac_overlap_hs[icomp]*pi6 * rho[icomp] * hs_diam_cubed;
      xsi1 += Fac_overlap_hs[icomp]*pi6 * rho[icomp] * hs_diam_cubed * Sigma_ff[icomp][icomp];
      xsi2 += Fac_overlap_hs[icomp]*pi6 * rho[icomp] * hs_diam_cubed * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
      xsi3 += Fac_overlap_hs[icomp]*pi6 * rho[icomp] * hs_diam_cubed * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);
   }
   y1 = 1.0 - xsi3;
   y2 = y1 * y1;
   y3 = y1 * y1 * y1;

   /* the hard sphere pressure in units of kT and Sigma_ff[1]^3 */

   betap_hs = (1.0/pi6) * (xsi0/y1 + 3.0 * xsi1 * xsi2/y2 +
                                     3.0 * POW_DOUBLE_INT(xsi2,3)/y3  );

   return (betap_hs);
}
/********************************************************************************
chempot_PY_hs:  This routine calculates the excess chemical potential for 
                hard spheres at the density of interest using the PY equations.            */

void chempot_PY_hs(double *rho)
{
   int icomp;
   double pi6, hs_diam_cubed, xsi0, xsi1, xsi2, xsi3, y1, y2, y3;

   xsi0=xsi1=xsi2=xsi3=0.0;
   pi6 = PI/6.0;                 /* shorthand  for pi/6                 */

   /*  calculate the constants xsi and introduce some shorthand */

   for (icomp=0; icomp<Ncomp; ++icomp) {
      hs_diam_cubed = POW_DOUBLE_INT(HS_diam[icomp],3);
      xsi0 += Fac_overlap_hs[icomp]*pi6 * rho[icomp] * hs_diam_cubed;
      xsi1 += Fac_overlap_hs[icomp]*pi6 * rho[icomp] * hs_diam_cubed * Sigma_ff[icomp][icomp];
      xsi2 += Fac_overlap_hs[icomp]*pi6 * rho[icomp] * hs_diam_cubed * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2);
      xsi3 += Fac_overlap_hs[icomp]*pi6 * rho[icomp] * hs_diam_cubed * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3);
   }
   y1 = 1.0 - xsi3;
   y2 = y1 * y1;
   y3 = y1 * y1 * y1;

   /* the excess hard sphere chemical potential in units of kT */

   for (icomp=0; icomp<Ncomp; ++icomp) 
      Betamu_hs_ex[icomp] = - log(y1) +
              pi6 * pressure_PY_hs(rho) * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],3) +
              3.0 * xsi2 * Sigma_ff[icomp][icomp]/y1 +
              3.0 * xsi1 * POW_DOUBLE_INT(Sigma_ff[icomp][icomp],2)/y1 +
              4.5 * POW_DOUBLE_INT((xsi2 * Sigma_ff[icomp][icomp]),2)/y2 ;

   return;
}
/*******************************************************************************/
/* compute_bulk_FMT_properties: compute some additional bulk values of
   nonlocal densities.  It is convenient to compute these once up front.  */
void compute_bulk_FMT_properties(char *output_file1)
{
  int i,loc_inode,loc_i,inode_box,inode,ijk[3],icomp,idim,iunk,printproc;
  int ibond,iseg,jseg,pol_number,type_jseg,nloop,iloop;
  double vol,area,x_dist;
  FILE *fp2=NULL;
  if (Proc==0) printproc = TRUE;
  else printproc=FALSE;
  if (printproc) {
    if( (fp2 = fopen(output_file1,"a+"))==NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
    }
  }
   /* compute bulk nonlocal densities needed for Rosenfeld terms */
  for (iunk=0; iunk<Nrho_bar; iunk++){
     Rhobar_b[iunk] = 0.0;
     Rhobar_b_LBB[iunk] = 0.0;
     Rhobar_b_RTF[iunk] = 0.0;
     Dphi_Drhobar_b[iunk]=0.0;
     Dphi_Drhobar_LBB[iunk]=0.0;
     Dphi_Drhobar_RTF[iunk]=0.0;
  }

  if (Type_poly==WTC) nloop=Nseg_tot;
  else             nloop=Ncomp;

  for (iloop=0; iloop<nloop; iloop++){
       if (Type_poly==WTC) icomp=Unk2Comp[iloop];
       else             icomp=iloop;

       if (Type_poly==WTC){
           rhobar_icomp(Rho_seg_b[iloop],icomp,Rhobar_b);
       }
       else{
          if (Lsteady_state){
                 rhobar_icomp(Rho_b_LBB[icomp],icomp,Rhobar_b_LBB);
                 rhobar_icomp(Rho_b_RTF[icomp],icomp,Rhobar_b_RTF);
          }
          else if (Nwall ==0 && Iliq_vap==3){
                 rhobar_icomp(Rho_coex[1],icomp,Rhobar_b_LBB);
                 rhobar_icomp(Rho_coex[0],icomp,Rhobar_b_RTF);
          }
          else rhobar_icomp(Rho_b[icomp],icomp,Rhobar_b);
   
       }
  }
  dphi_drb_bulk(Rhobar_b,Dphi_Drhobar_b);
  if (Lsteady_state || (Nwall==0&&Iliq_vap==3)){
     dphi_drb_bulk(Rhobar_b_LBB,Dphi_Drhobar_LBB);
     dphi_drb_bulk(Rhobar_b_RTF,Dphi_Drhobar_RTF);
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
        fclose(fp2);
  }
  return;
}
/*********************************************************************************************/
void rhobar_icomp(double rho,int icomp, double *rhobar)
{
   double vol,area,hs_diam;
   int idim;

   hs_diam = HS_diam[icomp];

   vol = PI*hs_diam*hs_diam*hs_diam/6.0;
   area = PI*hs_diam*hs_diam;

   rhobar[0] += Fac_overlap_hs[icomp]*vol*rho;
   rhobar[1] += Fac_overlap_hs[icomp]*area*rho;
   rhobar[2] += Fac_overlap_hs[icomp]*area*rho*Inv_4pir[icomp];
   rhobar[3] += Fac_overlap_hs[icomp]*area*rho*Inv_4pirsq[icomp];
   for (idim=0; idim<Ndim; idim++){
       rhobar[4+idim] = 0.0;
       rhobar[4+Ndim+idim] = 0.0;
   }
   return;
}
/*****************************************************************************/
void dphi_drb_bulk(double *rhobar,double *dphi_drb)
{
 double n[4+2*NDIM_MAX], inv_n3[5];
 int idim;

  n[3]=rhobar[0]; n[2]=rhobar[1];
  n[1]=rhobar[2]; n[0]=rhobar[3];
  for (idim=0;idim<Ndim;idim++){
        n[Nrho_bar_s+Ndim+idim]=0.;
        n[Nrho_bar_s+idim]=0.;
  }

  inv_n3[0] = (1.0 - n[3]);
  inv_n3[1] = 1.0 / inv_n3[0];
  inv_n3[2] = inv_n3[1]*inv_n3[1];
  inv_n3[3] = inv_n3[2]*inv_n3[1];
  inv_n3[4] = inv_n3[3]*inv_n3[1];

  FMT1stDerivBulk_switch(n,inv_n3,dphi_drb);

  return;
}
/*************************************************************************
dp_drho_hs_PY: the derivative of the hard sphere pressure with respect to rho
               in the PY approximation */
double dp_drho_hs_PY(double *rho)
{
   int icomp;
   double hs_diam_cubed,pi6;
   double xsi0,xsi1,xsi2,xsi3,y1,y2,y3,dy1_drho,dy2_drho,dy3_drho;
   double dp_drho;

   xsi0 = xsi1 = xsi2 = xsi3 = 0.0;
   pi6 = PI/6.0;                 /* shorthand  for pi/6                 */

   /*  calculate the constants xsi and introduce some shorthand */

   for (icomp=0; icomp<Ncomp; ++icomp) {
      hs_diam_cubed = POW_DOUBLE_INT(HS_diam[icomp],3);
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
dmu_drho_hs_PY: the derivative of the hard sphere chemical potential 
             with respect to rho - PY hard spheres*/
double dmu_drho_hs_PY(double *rho)
{
   int icomp;
   double hs_diam_cubed,pi6;
   double xsi0,xsi1,xsi2,xsi3,y1,y2,y3,dy1_drho,dy2_drho,dy3_drho;
   double dmu_drho;

   xsi0 = xsi1 = xsi2 = xsi3 = 0.0;
   pi6 = PI/6.0;                 /* shorthand  for pi/6                 */

   /*  calculate the constants xsi and introduce some shorthand */

   for (icomp=0; icomp<Ncomp; ++icomp) {
      hs_diam_cubed = POW_DOUBLE_INT(HS_diam[icomp],3);
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
        pi6 * dp_drho_hs_PY(rho) * POW_DOUBLE_INT(Sigma_ff[0][0],3) +
       (3.0 * xsi2 * Sigma_ff[0][0]/y1)
                                 *( 1.0/rho[0] - dy1_drho/y1) +
       (3.0 * xsi1 * POW_DOUBLE_INT(Sigma_ff[0][0],2)/y1)
                                 *( 1.0/rho[0] - dy1_drho/y1) +
       (9.0 * POW_DOUBLE_INT((xsi2 * Sigma_ff[0][0]),2)/y2)
                                 *( 1.0/rho[0] - 0.5*dy2_drho/y2) +
       (1.0/rho[0]);

   return(dmu_drho);
}
/*************************************************************************/
