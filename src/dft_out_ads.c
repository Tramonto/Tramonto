/* New version of dft_out_ads.c ---- testing a function pointer approach for the output */

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/**************************************************************************************/
void calc_adsorption(FILE *fp,double **x)
{
  int icomp,iunk;
  double ads[NCOMP_MAX],ads_ex[NCOMP_MAX],ads_b[NCOMP_MAX];
  static int first=TRUE;

  Integration_profile=NULL;

  if (Proc==0 &&Iwrite != NO_SCREEN) printf("-------------------- ADSORPTION -------------------------------\n");
  for (icomp=0;icomp<Ncomp;icomp++) ads[icomp]=0.0;
  for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
     if (Type_poly==WTC) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else icomp = iunk-Phys2Unk_first[DENSITY];
     integrateInSpace(&integrand_adsorption,iunk,Nel_hit2,x,Integration_profile);
     ads[icomp]+=Temporary_sum;
  }
 if (Proc==0 && Iwrite != NO_SCREEN){
     for (icomp=0;icomp<Ncomp;icomp++){
          print_to_screen_comp(icomp,ads[icomp],"ADSORPTION");
          print_to_file_comp(fp,icomp,ads[icomp],"ads",first);
      }
  }    

  for (icomp=0;icomp<Ncomp;icomp++) ads_b[icomp]=0.0;
  for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
     if (Type_poly==WTC) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else icomp = iunk-Phys2Unk_first[DENSITY];
     integrateInSpace(&integrand_adsorption_bulk,iunk,Nel_hit,x,Integration_profile);
     ads_b[icomp]+=Temporary_sum;
     ads_ex[icomp]=ads[icomp]-ads_b[icomp];
  }
  if (Proc==0 && Iwrite != NO_SCREEN){
     for (icomp=0;icomp<Ncomp;icomp++){
        print_to_screen_comp(icomp,ads_ex[icomp],"EXCESS ADSORPTION");
        print_to_file_comp(fp,icomp,ads_ex[icomp],"ads_ex",first);
     }    
  }
  if (first) first=FALSE;
  if (Proc==0 &&Iwrite != NO_SCREEN) printf("---------------------------------------------------------------\n");
  return;
}
/**************************************************************************************/
void calc_fluid_charge(FILE *fp,double **x)
{
 static int first=TRUE;
 if (Proc==0&&Iwrite != NO_SCREEN) printf("-------------------- CHARGE     -------------------------------\n");

 integrateInSpace_SumInComp(&integrand_fluid_charge,Nel_hit2,x);

 if (Proc==0 && Iwrite != NO_SCREEN){
      print_to_screen(Temporary_sum,"CHARGE IN FLUID");
      print_to_file(fp,Temporary_sum,"charge",first);
      printf("---------------------------------------------------------------\n");
 }
 if (first) first=FALSE;
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

