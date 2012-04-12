using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_functionals_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List )
{
  /**************************************/
  /* Define validators for this section.*/
  /**************************************/
  RCP<StringValidator> HSFuncValidator = rcp(
          new StringValidator(tuple<std::string>("Ideal Fluid / No Volume Exclusion","HS_Rosenfeld_Original", "HS_Rosenfeld_Zerocrossover","HS_White_Bear","HS_White_Bear2")));

  RCP<StringValidator> FuncPairPotValidator = rcp(
          new StringValidator(tuple<std::string>("No Mean Field Functional","UATT_core=[Umin; r=0,Rmin]","UATT_core=[0;r=0,Sigma]",
               "UATT_core=[U(sigma);r=0,Sigma]","UATT_core=[(0;r=0,Sigma),(Umin;r=Sigma,Rmin)]","UATTcore=[0;r=0,Rzero]")));

  RCP<StringValidator> CHARGEValidator = rcp(
          new StringValidator(tuple<std::string>("No Charge or No Poisson","Charge_Mean_Field","Charge_with_DeltaC_RPM", "Charge_with_DeltaC_General","Charge(MF)_with_Polarization")));

  RCP<StringValidator> POLYValidator = rcp(
          new StringValidator(tuple<std::string>("No Bonded Fluids","Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT","Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

  RCP<StringValidator> CalcTypeValidator = rcp(
         new StringValidator(tuple<std::string>("Equilibrium (homogeneous boundary conditions)","Equilibrium (inhomogeneous boudary conditions)","Steady State Diffusion (inhomogeneous boundaries)")));

  /*********************/
  /* set up parameters */
  /*********************/

  Functional_List->set("F0_Type_of_Calculation", "Equilibrium (homogeneous boundary conditions)", "Indicate the basic type of calculation to be performed.", CalcTypeValidator);
  Functional_List->set("F1_HS_Functional", "Ideal Fluid / No Volume Exclusion", "Type of Hard Sphere Functional (set to NONE for ideal gas or Poisson-Boltzman electrolyte)", HSFuncValidator);
  Functional_List->set("F2_PAIRPOTcore_Functional", "No Mean Field Functional", "Type of Core treatment for an extended pair potential \n using a strict mean field functional", FuncPairPotValidator);
  Functional_List->set("F3_CHARGE_Functional", "No Charge or No Poisson", "Functional for Charged systems (set to NONE if there are no charges in the system)", CHARGEValidator);
  Functional_List->set("F4_POLYMER_Functional", "No Bonded Fluids", "Functional for treating bonded systems", POLYValidator);

  return;
}
/************************************************************************************************************/
void dft_GUI_functionals_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List )
{
  string type_setting_tmp;

  /**************************************/
  /* Define validators for this section.*/
  /**************************************/
  RCP<StringValidator> HSFuncValidator = rcp(
         new StringValidator(tuple<std::string>("Ideal Fluid / No Volume Exclusion","HS_Rosenfeld_Original", "HS_Rosenfeld_Zerocrossover","HS_White_Bear","HS_White_Bear2")));

  RCP<StringValidator> FuncPairPotValidator = rcp(
         new StringValidator(tuple<std::string>("No Mean Field Functional","UATT_core=[Umin; r=0,Rmin]","UATT_core=[0;r=0,Sigma]",
              "UATT_core=[U(sigma);r=0,Sigma]","UATT_core=[(0;r=0,Sigma),(Umin;r=Sigma,Rmin)]","UATTcore=[0;r=0,Rzero]")));

  RCP<StringValidator> CHARGEValidator = rcp(
         new StringValidator(tuple<std::string>("No Charge or No Poisson","Charge_Mean_Field","Charge_with_DeltaC_RPM", "Charge_with_DeltaC_General","Charge(MF)_with_Polarization")));

  RCP<StringValidator> POLYValidator = rcp(
         new StringValidator(tuple<std::string>("No Bonded Fluids","Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT","Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

  RCP<StringValidator> CalcTypeValidator = rcp(
         new StringValidator(tuple<std::string>("Equilibrium (homogeneous boundary conditions)","Equilibrium (inhomogeneous boudary conditions)","Steady State Diffusion (inhomogeneous boundaries)")));

  /*********************/
  /* set up parameters */
  /*********************/

  if (Type_interface==UNIFORM_INTERFACE) type_setting_tmp="Equilibrium (homogeneous boundary conditions)";
  else if (Type_interface==PHASE_INTERFACE) type_setting_tmp="Equilibrium (inhomogeneous boudary conditions)";
  else if (Type_interface==DIFFUSIVE_INTERFACE) type_setting_tmp="Steady State Diffusion (inhomogeneous boundaries)";
  Functional_List->set("F0_Type_of_Calculation", type_setting_tmp, "Indicate the basic type of calculation to be performed.", CalcTypeValidator);

  if (Type_func==NONE) type_setting_tmp="Ideal Fluid / No Volume Exclusion";
  else if (Type_func==FMT1) type_setting_tmp="HS_Rosenfeld_Original";
  else if (Type_func==FMT2) type_setting_tmp="HS_Rosenfeld_Zerocrossover";
  else if (Type_func==FMT3) type_setting_tmp="HS_White_Bear";
  else if (Type_func==FMT4) type_setting_tmp="HS_White_Bear2";
  Functional_List->set("F1_HS_Functional", type_setting_tmp, "Type of Hard Sphere Functional (set to NONE for ideal gas or Poisson-Boltzman electrolyte)", HSFuncValidator);

  if (Type_attr==NONE) type_setting_tmp="No Mean Field Functional";
  else if (Type_attr==MFPAIR_RMIN_UMIN) type_setting_tmp="UATT_core=[Umin; r=0,Rmin]";
  else if (Type_attr==MFPAIR_SIGMA_ZERO) type_setting_tmp="UATT_core=[0;r=0,Sigma]";
  else if (Type_attr==MFPAIR_SIGMA_USIGMA) type_setting_tmp="UATT_core=[U(sigma);r=0,Sigma]";
  else if (Type_attr==MFPAIR_SIGTOUMIN_UMIN) type_setting_tmp="UATT_core=[(0;r=0,Sigma),(Umin;r=Sigma,Rmin)]";
  else if (Type_attr==MFPAIR_RCSZERO_ZERO) type_setting_tmp="UATTcore=[0;r=0,Rzero]";
  Functional_List->set("F2_PAIRPOTcore_Functional", type_setting_tmp, "Type of Core treatment for an extended pair potential \n using a strict mean field functional", FuncPairPotValidator);

  if (Type_coul==NONE) type_setting_tmp="No Charge or No Poisson";
  else if (Type_coul==BARE) type_setting_tmp="Charge_Mean_Field";
  else if (Type_coul==DELTAC_RPM) type_setting_tmp="Charge_with_DeltaC_RPM";
  else if (Type_coul==DELTAC_GENERAL) type_setting_tmp="Charge_with_DeltaC_General";
  else if (Type_coul==POLARIZE) type_setting_tmp="Charge(MF)_with_Polarization";
  Functional_List->set("F3_CHARGE_Functional", type_setting_tmp, "Functional for Charged systems (set to NONE if there are no charges in the system)", CHARGEValidator);

  if (Type_poly==NONE) type_setting_tmp="No Bonded Fluids";
  else if (Type_poly==CMS) type_setting_tmp="Polymer_CMS";
  else if (Type_poly==CMS_SCFT) type_setting_tmp="Polymer_CMS_SCFT";
  else if (Type_poly==WTC) type_setting_tmp="Polymer_TC_iSAFT";
  else if (Type_poly==WJDC) type_setting_tmp="Polymer_JDC_iSAFT(seg)";
  else if (Type_poly==WJDC2) type_setting_tmp="Polymer_JDC_iSAFT(segRho compField)";
  else if (Type_poly==WJDC3) type_setting_tmp="Polymer_JDC_iSAFT(comp)";
  Functional_List->set("F4_POLYMER_Functional", type_setting_tmp, "Functional for treating bonded systems", POLYValidator);

  return;
}
/************************************************************************************************************/
void dft_GUI_functionals_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List )
{
  /* No dependencies here */

   return;
}
/************************************************************************************************************/
