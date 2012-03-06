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

/* -----------------------------------------------------------------------------------------
dft_thermo_elec_MSAgeneral.c: Calculate the thermodynamic properties for a charged system 
where generalized MSA correlations are turned on.  From 
Oleksy and Hansen, Mol Phys. v.104, 2871-2883
--------------------------------------------------------------------------------------------*/
#include "dft_thermo_elec_MSAgeneral.h"

/********************************************************************
chempot_ELEC_MSA_GENERAL: Here we compute the chemical potential contribution due
      to cross correlations between the hard sphere and coulomb parts of 
      the potential. */
void chempot_ELEC_MSA_GENERAL(double *rho)
{
   int icomp,jcomp;

   Deltac_b = (double *) array_alloc (1, Ncomp, sizeof(double));
   for (icomp=0; icomp<Ncomp; icomp++) Deltac_b[icomp] = 0.0;

   for (icomp=0; icomp<Ncomp; icomp++) 
      for (jcomp=0; jcomp<Ncomp;jcomp++){
          Deltac_b[icomp]+= rho[jcomp]*int_stencil_bulk(THETA_CR_GENERAL_MSA,icomp,jcomp,NULL);
      }
   return;
}
/*******************************************************************************/
/*pressure_elec_MSA: Here we compute the chemical potential contribution due
      to cross correlations between the hard sphere and coulomb parts of
      the potential. */
double pressure_elec_MSA(double *rho)
{
   int icomp;
   double betap_elec;
   
   betap_elec=0.0;

   for (icomp=0; icomp<Ncomp; icomp++){ betap_elec -= 0.5*rho[icomp]*Deltac_b[icomp]; }
   return(betap_elec);
}
/*******************************************************************************/
/* deltaC_GENERAL_MSA:  calculate deltac_MSA generalized for ion size*/

double deltaC_GENERAL_MSA(double r,int i, int j)
{
  double deltac;
  double NplusGammaX_i,NplusGammaX_itemp;
  int itemp,jtemp;

  /* note that the parameters X_MSA, N_SMA, Gamma_MSA are computed in the routine
     precalc_GENmsa_params wich is called at the top of dft_stencil.c outside the quadrature 
     loops of the stencil computation */
 
  if (HS_diam[j]<HS_diam[i]){ itemp=j; jtemp=i; }
  else                      { itemp=i; jtemp=j; }

  NplusGammaX_i=N_MSA[i]+Gamma_MSA*X_MSA[i];
  NplusGammaX_itemp=N_MSA[itemp]+Gamma_MSA*X_MSA[itemp];

  if (r == 0.0) printf("trouble with deltaC term .... r=0");

  if (r <= (HS_diam[jtemp]-HS_diam[itemp])/2.0 && r>0) {
     deltac= (2.0/Temp_elec)*(
                 (Charge_f[itemp]*Charge_f[jtemp]/(2.0*r))
                 +Charge_f[itemp]*N_MSA[jtemp]
                 -X_MSA[itemp]*(NplusGammaX_itemp)
                 +(HS_diam[itemp]/3.0)*POW_DOUBLE_INT(NplusGammaX_itemp,2));
  }
  else if (r <= (HS_diam[i]+HS_diam[j])/2.0) {

     deltac = (1.0/Temp_elec)*
              ( (Charge_f[i]*Charge_f[j] +(HS_diam[i]-HS_diam[j])*MSAgen_term1[i][j])/r
               -MSAgen_term2[i][j]
               +r*MSAgen_term3[i][j]
               +r*r*r*MSAgen_term4[i][j]
              );
  }
  else deltac = 0.0;

  return deltac;
}
/*******************************************************************************/
/* deltaC_MSAgeneral_int:  given range of calculate the 
           integral of deltac_MSA at r */

double deltaC_GENERAL_MSA_int(double r,int i, int j)
{
  double deltac_int;
  double NplusGammaX_i,NplusGammaX_itemp;
  int itemp,jtemp;


  if (HS_diam[j]<HS_diam[i]){ itemp=j; jtemp=i; }
  else                      { itemp=i; jtemp=j; }

  NplusGammaX_itemp=N_MSA[i]+Gamma_MSA*X_MSA[i];   /* note that this needs to be checked */

  if (r <= (HS_diam[jtemp]-HS_diam[itemp])/2.0 && r>0) {
     deltac_int= (8.0*PI*r*r/Temp_elec)*(
              (Charge_f[itemp]*Charge_f[jtemp]/4.0)
             +(r/3.0)*( Charge_f[itemp]*N_MSA[jtemp]
                    -X_MSA[itemp]*(NplusGammaX_itemp)
                    +(HS_diam[itemp]/3.0)*POW_DOUBLE_INT(NplusGammaX_itemp,2)) );
  }
  else if (r <= (HS_diam[i]+HS_diam[j])/2.0) {

     deltac_int = (4*PI*r*r)*(1.0/Temp_elec)*
               ( 0.5*(Charge_f[i]*Charge_f[j])+0.5*(HS_diam[i]-HS_diam[j])*MSAgen_term1[i][j]
                -r*MSAgen_term2[i][j]/3.0
                +r*r*MSAgen_term3[i][j]/4.0
                +r*r*r*r*MSAgen_term4[i][j]/6.0
               );
  }
  return deltac_int;
}
/*******************************************************************************/
/* find_MSAgeneral_params(double *rho, double *x_msa, double *n_msa, double gamma)
/* Iterate to find parameters needed for generalized MSA approach.*/
void precalc_GENmsa_params(double *rho, double *x_msa, double *n_msa, double gamma)
{
  int i,j,iter=0;
  double c,densitySum,sum_num,sum_denom,xsum,nsum;
  double tol=1.e-12,error=1.0,xsum_old=0.0,nsum_old=0.0,gamma_old=0.0;
  double NplusGammaX_i,NplusGammaX_j;
   

  gamma=1.0;
  c = 0.0;
  for (i=0;i<Ncomp;i++){
     x_msa[i]=1.0;
     n_msa[i]=1.0;
     c += (PI/2.0)*(1./(1.-(PI/6.)*rho[i]*POW_DOUBLE_INT(HS_diam[i],3)));
  }

  while (error>tol && iter <10000){

      densitySum=0.0;
      for (i=0;i<Ncomp;i++){
         sum_num=0.0;
         sum_denom=0.0;
         for (j=0;j<Ncomp;j++){
             sum_num+=rho[j]*HS_diam[j]*Charge_f[j]/(1.0+gamma*HS_diam[j]);
             sum_denom+=rho[j]*POW_DOUBLE_INT(HS_diam[j],3)/(1.0+gamma*HS_diam[j]);
         }

         x_msa[i]=Charge_f[i]/(1+gamma*HS_diam[i])
                 -(c*HS_diam[i]*HS_diam[i]/(1+gamma*HS_diam[i]))*
                 sum_num/(1+c*sum_denom);

         n_msa[i]=(x_msa[i]-Charge_f[i])/HS_diam[i];

         densitySum+=rho[i]*x_msa[i]*x_msa[i];
      }
      gamma=sqrt((PI/Temp_elec)*densitySum);

      xsum=0.0; nsum=0.0;
      for (i=0;i<Ncomp;i++) {
            xsum+=x_msa[i];
            nsum+=n_msa[i];
      }

      error=sqrt( (gamma-gamma_old)*(gamma-gamma_old)
             + (xsum-xsum_old)*(xsum-xsum_old) 
             + (nsum-nsum_old)*(nsum-nsum_old) );
      
      gamma_old=gamma; xsum_old=xsum; nsum_old=nsum;
      iter++;
  }
  if (error >tol || iter ==10,000){
      printf("unable to converge general MSA properties\n");
      printf("iter=%d  tol=%g  \n",iter,tol);
      for (i=0;i<Ncomp;i++) printf("icomp=%d  gamma=%g  x_msa[%d]=%g  n_msa[%d]=%g\n",
                                                       i,gamma,i,x_msa[i],i,n_msa[i]);
      exit(-1);
  }
  else{
    if (Proc==0 &&Iwrite==VERBOSE){
        printf("converged values of MS params are:\n");
        for (i=0;i<Ncomp;i++)
        printf("icomp=%d  X_MSA[icomp]=%g  N_MSA[icomp]=%g Gamma=%g || iter=%d error=%g\n",i,x_msa[i],n_msa[i],gamma, iter,error);
    }
  }

  for (i=0;i<Ncomp;i++){
      for (j=0;j<Ncomp;j++){
        NplusGammaX_i=n_msa[i]+gamma*x_msa[i];
        NplusGammaX_j=n_msa[j]+gamma*x_msa[j];

        MSAgen_term1[i][j]=((x_msa[i]+x_msa[j])/4.0)*((NplusGammaX_i)-(NplusGammaX_j))
                           -((HS_diam[i]-HS_diam[j])/16.0)*
                           (POW_DOUBLE_INT(NplusGammaX_i+NplusGammaX_j,2)-4*n_msa[i]*n_msa[j]);

        MSAgen_term2[i][j]=(x_msa[i]-x_msa[j])*(n_msa[i]-n_msa[j])
                           +(HS_diam[i]+HS_diam[j])*n_msa[i]*n_msa[j]
                           +(x_msa[i]*x_msa[i] + x_msa[j]*x_msa[j])*gamma
                           -(1.0/3.0)*(HS_diam[i]*POW_DOUBLE_INT(NplusGammaX_i,2) +HS_diam[j]*POW_DOUBLE_INT(NplusGammaX_j,2));

        MSAgen_term3[i][j]= NplusGammaX_i* x_msa[i]/HS_diam[i]
                           + NplusGammaX_j* x_msa[j]/HS_diam[j]
                           + n_msa[i]*n_msa[j]
                           - 0.5*(POW_DOUBLE_INT(NplusGammaX_i,2)+ POW_DOUBLE_INT(NplusGammaX_j,2));

        MSAgen_term4[i][j]= (POW_DOUBLE_INT(NplusGammaX_i,2)/(6.0*HS_diam[i]*HS_diam[i]))
                           +(POW_DOUBLE_INT(NplusGammaX_j,2)/(6.0*HS_diam[j]*HS_diam[j]));
     
      }
  }

  return;
}
/*******************************************************************************/
