/* New version of dft_out_ads.c ---- testing a function pointer approach for the output */

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/**************************************************************************************/
void calc_adsorption(FILE *fp,double **x)
{
  int i,icomp,iunk;
  double ads[NCOMP_MAX],ads_ex[NCOMP_MAX],ads_b[NCOMP_MAX];;

   for (i=0;i<Imax;i++){

      for (icomp=0;icomp<Ncomp;icomp++) ads[icomp]=0.0;
      for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
         if (Type_poly==WTC) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
         else icomp = iunk-Phys2Unk_first[DENSITY];
         integrateInSpace(&integrand_adsorption,iunk,i,Nel_hit2,x);
         ads[icomp]+=Temporary_sum;
      }
      for (icomp=0;icomp<Ncomp;icomp++){
         if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen_comp(i,icomp,ads[icomp],"ADSORPTION");
/*         print_to_output_file(i,iunk,ads[icomp],"ADSORPTON");*/
      }    

      for (icomp=0;icomp<Ncomp;icomp++) ads_b[icomp]=0.0;
      for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
         if (Type_poly==WTC) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
         else icomp = iunk-Phys2Unk_first[DENSITY];
         integrateInSpace(&integrand_adsorption_bulk,iunk,i,Nel_hit,x);
         ads_b[icomp]+=Temporary_sum;
         ads_ex[icomp]=ads[icomp]-ads_b[icomp];
      }
      for (icomp=0;icomp<Ncomp;icomp++){
         if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen_comp(i,icomp,ads_ex[icomp],"EXCESS ADSORPTION");
/* print_to_output_file(i,iunk,ads_ex[icomp],"ADS_EX");*/
      }    
   }  
   if (Iwrite != NO_SCREEN) printf("=====================================================\n");
   return;
}
/**************************************************************************************/
void calc_fluid_charge(FILE *fp,double **x)
{
  int i;

  for (i=0;i<Imax;i++){
    integrateInSpace_SumInComp(&integrand_fluid_charge,i,Nel_hit2,x);
    if (Proc==0 && Iwrite != NO_SCREEN) print_to_screen(i,Temporary_sum,"CHARGE IN FLUID");
/*         print_to_output_file(i,val,"CHARGE_F");*/
  }    

  if (Iwrite != NO_SCREEN) printf("=====================================================\n");
  return;
}
/**************************************************************************************/
double integrand_adsorption(int iunk,int inode_box, double **x)
{
     double integrand;
     integrand = x[iunk][inode_box];
     return(integrand);
}
/**************************************************************************************/
double integrand_adsorption_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;

     if (Type_poly==WTC) rho_bulk = Rho_seg_b[iunk-Phys2Unk_first[DENSITY]];
     else                rho_bulk = Rho_b[iunk-Phys2Unk_first[DENSITY]];

     integrand = rho_bulk;
     return(integrand);
}
/**************************************************************************************/
double integrand_fluid_charge(int iunk, int inode_box, double **x)
{
     double integrand;
     int icomp;

     if (Type_poly==WTC) icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else                icomp = iunk-Phys2Unk_first[DENSITY];

     integrand = Charge_f[icomp]*x[iunk][inode_box];
     return(integrand);
}
/**************************************************************************************/

