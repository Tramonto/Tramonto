//using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;
int func_testNcomp(int);

void dft_GUI_potentialsFF_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List)
{
  /**************************************/
  /* Define validators for this section.*/
  /**************************************/

  RCP<EnhancedNumberValidator<int> > NcompValidator = rcp(new EnhancedNumberValidator<int>(1,1000,1));

  RCP<StringValidator> HSDiamTypeValidator = rcp(
         new StringValidator(tuple<std::string>("HSdiam=Sigma","Manual Definition","Barker-Henderson")));

/*    RCP<EnhancedNumberValidator<double> > HSDiamValidator = rcp(new EnhancedNumberValidator<double>(0.0,100.,0.1));*/

  RCP<StringValidator> MixTypeValidator = rcp(
         new StringValidator(tuple<std::string>("Lorentz-Berthlot Mixing","Manual Definition")));

  RCP<StringValidator> PairPotValidator = rcp(
         new StringValidator(tuple<std::string>("none",
         "LJ 12-6 potential (cut/shift)",
         "Coulomb potential as mean field (cut/shift)",
         "Coulomb potential as mean field (cut only)",
         "Yukawa potential (cut/shift)",
         "Exponential potential (cut/shift)",
         "Square Well potential",
         "LJ 12-6 plus Yukawa potential (cut/shift)",
         "r^12 repulsion plus Yukawa potential (cut/shift)",
         "r^18 repulsion plus Yukawa potential (cut/shift)",
         "r^N repulsion plus Yukawa potential (cut/shift)")));

  /***************************************************************/
  /* set up fluid parameters */
  /***************************************************************/
  

  Fluid_List->set("F1_Ncomp", 1, "Number of fluid components in the problem \n(Note: For co-polymers, each bead type is a different component)", NcompValidator);

  Fluid_List->set("F2_HSDiamType", "HSdiam=Sigma", "Select method for defining hard sphere diameters", HSDiamTypeValidator);

  Array<double> HSDiam_Array( (Fluid_List->get<int>("F1_Ncomp")),1.0);
  Fluid_List->set("F3_HSDiam", HSDiam_Array, "Set hard sphere diameters manually: Ncomp entries");

  Fluid_List->set("F4_PairPotType", "none", "Type of pair interaction potential", PairPotValidator);

  Fluid_List->set("F5_Dimensionless_Energy_Entry", true, "Dimensionless energies must be provided in units of epsilon/kT.\n Otherwise, set to false for entry as epsilon/k (units of K)");
  Fluid_List->set("F6_Temperature",1.0, "Enter the Temperature (units of K)");


     /* items for PotentialsFF_sublist */

 PotentialsFF_List->set("PF0_Off_Diagonal_Definitions", "Manual Definition", "Options for definition of off-diagonal pair potential parameters", MixTypeValidator);

 Array<double> Mass_Array( Fluid_List->get<int>("F1_Ncomp"),1.0);
 PotentialsFF_List->set("PF7_Mass", Mass_Array, "Mass of each component");

 Array<double> Charge_Array( Fluid_List->get<int>("F1_Ncomp"),0.0);
 PotentialsFF_List->set("PF8_Charge", Charge_Array, "Charge of each component");

 Array<double> Polarization_Array( Fluid_List->get<int>("F1_Ncomp"),0.0);
 PotentialsFF_List->set("PF9_Polarization", Polarization_Array, "Polarizability of each component");

 Array<double> SigmaF_Array((Fluid_List->get<int>("F1_Ncomp")),1.0);
 PotentialsFF_List->set("PF1_SigmaF", SigmaF_Array, "Sigma - Characteristic diameter of fluid interactions (diagonal entries only)");

 Array<double> EpsF_Array( (Fluid_List->get<int>("F1_Ncomp")),1.0);
 PotentialsFF_List->set("PF2_EpsF", EpsF_Array, "Eps - Characteristic energy prefactor for fluid interactions (diagonal entries only)");

 Array<double> CutF_Array( (Fluid_List->get<int>("F1_Ncomp")),3.0);
 PotentialsFF_List->set("PF3_CutF", CutF_Array, "Cutoff distance for fluid-fluid pair interactions (diagonal entries only)");

 Array<double> EpsYukawaF_Array( (Fluid_List->get<int>("F1_Ncomp")),1.0);
 PotentialsFF_List->set("PF4_EpsYukawaF", EpsYukawaF_Array, "AYuk - Energy prefactor for Yukawa term in fluid-fluid pair interactions (diagonal entries only)");

 Array<double> YukawaKF_Array( (Fluid_List->get<int>("F1_Ncomp")),1.0);
 PotentialsFF_List->set("PF5_ExpDecayParamF", YukawaKF_Array, "alpha - exponential decay parameter in fluid-fluid pair interactions (diagonal entries only)");


     /* set up for 2D arrays needed when we do a manual definition of interactions */

 TwoDArray<double> SigmaFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),1.0);
 PotentialsFF_List->set("PF1_SigmaFF", SigmaFF_Array, "Sigma - Characteristic diameter fluid-fluid pair interactions");

 TwoDArray<double> EpsFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),1.0);
 PotentialsFF_List->set("PF2_EpsFF", EpsFF_Array, "Eps - Energy prefactor for fluid-fluid pair interactions");

 TwoDArray<double> CutFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),3.0);
 PotentialsFF_List->set("PF3_CutFF", CutFF_Array, "rcut - Cutoff distance for fluid-fluid pair interactions");

 TwoDArray<double> EpsYukawaFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),1.0);
 PotentialsFF_List->set("PF4_EpsYukawaFF", EpsYukawaFF_Array, "AYuk - Energy prefactor for Yukawa term in fluid-fluid pair interactions");

 TwoDArray<double> YukawaKFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),1.0);
 PotentialsFF_List->set("PF5_ExpDecayParamFF", YukawaKFF_Array, "alpha - exponential decay parameter in fluid-fluid pair interactions");

 TwoDArray<double> NpowFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),12.0);
 PotentialsFF_List->set("PF6_NpowFF", NpowFF_Array, "n - variable power for r^-n repulsive term in Fluid-fluid pair interactions");

 TwoDArray<double> BondFF_Array( (Fluid_List->get<int>("F1_Ncomp")),(Fluid_List->get<int>("F1_Ncomp")),0.0);
 PotentialsFF_List->set("PF10_BondFF", BondFF_Array, "Bond Lengths for fluid-fluid pair interactions");

 return;
}
/**************************************************************************************************************/
int func_testNcomp(int Ncomp)
{
   return Ncomp-1;
}
/**************************************************************************************************************/
void dft_GUI_potentialsFF_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List )
{
  string tmp_string_var;
  double tmp_1D[NCOMP_MAX];
  int i,j;

  /**************************************/
  /* Define validators for this section.*/
  /**************************************/

  RCP<EnhancedNumberValidator<int> > NcompValidator = rcp(new EnhancedNumberValidator<int>(1,1000,1));

  RCP<StringValidator> HSDiamTypeValidator = rcp(
         new StringValidator(tuple<std::string>("HSdiam=Sigma","Manual Definition","Barker-Henderson")));

/*    RCP<EnhancedNumberValidator<double> > HSDiamValidator = rcp(new EnhancedNumberValidator<double>(0.0,100.,0.1));*/

  RCP<StringValidator> MixTypeValidator = rcp(
         new StringValidator(tuple<std::string>("Lorentz-Berthlot Mixing","Manual Definition")));

  RCP<StringValidator> PairPotValidator = rcp(
         new StringValidator(tuple<std::string>("none",
         "LJ 12-6 potential (cut/shift)",
         "Coulomb potential as mean field (cut/shift)",
         "Coulomb potential as mean field (cut only)",
         "Yukawa potential (cut/shift)",
         "Exponential potential (cut/shift)",
         "Square Well potential",
         "LJ 12-6 plus Yukawa potential (cut/shift)",
         "r^12 repulsion plus Yukawa potential (cut/shift)",
         "r^18 repulsion plus Yukawa potential (cut/shift)",
         "r^N repulsion plus Yukawa potential (cut/shift)")));

  /***************************************/
  /* Set parameters from Old Format File.*/
  /***************************************/

  Fluid_List->set("F1_Ncomp", Ncomp, "Number of fluid components in the problem \n(Note: For co-polymers, each bead type is a different component)", NcompValidator);

  if(Type_hsdiam==SIGMA_DIAM)          tmp_string_var="HSdiam=Sigma";
  else if(Type_hsdiam==MANUAL_HS_DIAM) tmp_string_var="Manual Definition";
  else if(Type_hsdiam==BH_DIAM)        tmp_string_var="Barker-Henderson";
  Fluid_List->set("F2_HSDiamType", tmp_string_var, "Select method for defining hard sphere diameters", HSDiamTypeValidator);

  Array<double> HSDiam_Array(HS_diam,HS_diam+Ncomp);
  Fluid_List->set("F3_HSDiam", HSDiam_Array, "Set hard sphere diameters manually: Ncomp entries");

  if (Type_pairPot==PAIR_HARD)                 tmp_string_var="none";
  else if (Type_pairPot==PAIR_LJ12_6_CS)       tmp_string_var="LJ 12-6 potential (cut/shift)";
  else if (Type_pairPot==PAIR_COULOMB_CS)      tmp_string_var="Coulomb potential as mean field (cut/shift)";
  else if (Type_pairPot==PAIR_COULOMB)         tmp_string_var="Coulomb potential as mean field (cut only)";
  else if (Type_pairPot==PAIR_YUKAWA_CS)       tmp_string_var="Yukawa potential (cut/shift)";
  else if (Type_pairPot==PAIR_EXP_CS)          tmp_string_var="Exponential potential (cut/shift)";
  else if (Type_pairPot==PAIR_SW)              tmp_string_var="Square Well potential";
  else if (Type_pairPot==PAIR_LJandYUKAWA_CS)  tmp_string_var="LJ 12-6 plus Yukawa potential (cut/shift)";
  else if (Type_pairPot==PAIR_r12andYUKAWA_CS) tmp_string_var="r^12 repulsion plus Yukawa potential (cut/shift)";
  else if (Type_pairPot==PAIR_r18andYUKAWA_CS) tmp_string_var="r^18 repulsion plus Yukawa potential (cut/shift)";
  else if (Type_pairPot==PAIR_rNandYUKAWA_CS)  tmp_string_var="r^N repulsion plus Yukawa potential (cut/shift)";
  Fluid_List->set("F4_PairPotType", tmp_string_var, "Type of pair interaction potential", PairPotValidator);

  if (Temp>1.0){
      Fluid_List->set("F5_Dimensionless_Energy_Entry", false, "Dimensionless energies must be provided in units of epsilon/kT.\n Otherwise, set to false for entry as epsilon/k (units of K)");
      Fluid_List->set("F6_Temperature",Temp, "Enter the Temperature (units of K)");
  }
  else{
      Fluid_List->set("F5_Dimensionless_Energy_Entry", true, "Dimensionless energies must be provided in units of epsilon/kT.\n Otherwise, set to false for entry as epsilon/k (units of K)");
      Fluid_List->set("F6_Temperature",1.0, "Enter the Temperature (units of K)");
  }
  if (Mix_type==1) PotentialsFF_List->set("PF0_Off_Diagonal_Definitions", "Manual Definition", "Options for definition of off-diagonal pair potential parameters", MixTypeValidator);
  else             PotentialsFF_List->set("PF0_Off_Diagonal_Definitions", "Lorentz-Berthlot Mixing", "Options for definition of off-diagonal pair potential parameters", MixTypeValidator);

  for (i=0;i<Ncomp;i++) tmp_1D[i]=Sigma_ff[i][i];
  Array<double> SigmaF_Array(tmp_1D, tmp_1D+Ncomp);
  PotentialsFF_List->set("PF1_SigmaF", SigmaF_Array, "Sigma - Characteristic diameter of fluid interactions (diagonal entries only)");

  for (i=0;i<Ncomp;i++) tmp_1D[i]=Eps_ff[i][i];
  Array<double> EpsF_Array(tmp_1D,tmp_1D+Ncomp);
  PotentialsFF_List->set("PF2_EpsF", EpsF_Array, "Eps - Characteristic energy prefactor for fluid interactions (diagonal entries only)");

  for (i=0;i<Ncomp;i++) tmp_1D[i]=Cut_ff[i][i];
  Array<double> CutF_Array(tmp_1D,tmp_1D+Ncomp);
  PotentialsFF_List->set("PF3_CutF", CutF_Array, "Cutoff distance for fluid-fluid pair interactions (diagonal entries only)");

  for (i=0;i<Ncomp;i++) tmp_1D[i]=EpsYukawa_ff[i][i];
  Array<double> EpsYukawaF_Array(tmp_1D,tmp_1D+Ncomp);
  PotentialsFF_List->set("PF4_EpsYukawaF", EpsYukawaF_Array, "AYuk - Energy prefactor for Yukawa term in fluid-fluid pair interactions");

  for (i=0;i<Ncomp;i++) tmp_1D[i]=YukawaK_ff[i][i];
  Array<double> YukawaKF_Array(tmp_1D,tmp_1D+Ncomp);
  PotentialsFF_List->set("PF5_ExpDecayParamF", YukawaKF_Array, "alpha - exponential decay parameter in fluid-fluid pair interactions");



  TwoDArray<double> SigmaFF_Array(Ncomp,Ncomp);
  for (i=0; i<Ncomp;i++) for (j=0; j<Ncomp;j++) SigmaFF_Array(i,j)=Sigma_ff[i][j];
  PotentialsFF_List->set("PF1_SigmaFF", SigmaFF_Array, "Sigma - Characteristic diameter fluid-fluid pair interactions");

  TwoDArray<double> EpsFF_Array(Ncomp,Ncomp);
  for (i=0; i<Ncomp;i++) for (j=0; j<Ncomp;j++) EpsFF_Array(i,j)=Eps_ff[i][j]; 
  PotentialsFF_List->set("PF2_EpsFF", EpsFF_Array, "Eps - Energy prefactor for fluid-fluid pair interactions");

  TwoDArray<double> CutFF_Array(Ncomp,Ncomp);
  for (i=0; i<Ncomp;i++) for (j=0; j<Ncomp;j++) CutFF_Array(i,j)=Cut_ff[i][j]; 
  PotentialsFF_List->set("PF3_CutFF", CutFF_Array, "rcut - Cutoff distance for fluid-fluid pair interactions");

  TwoDArray<double> EpsYukawaFF_Array(Ncomp,Ncomp);
  for (i=0; i<Ncomp;i++) for (j=0; j<Ncomp;j++) EpsYukawaFF_Array(i,j)=EpsYukawa_ff[i][j]; 
  PotentialsFF_List->set("PF4_EpsYukawaFF", EpsYukawaFF_Array, "AYuk - Energy prefactor for Yukawa term in fluid-fluid pair interactions");

  TwoDArray<double> YukawaKFF_Array(Ncomp,Ncomp);
  for (i=0; i<Ncomp;i++) for (j=0; j<Ncomp;j++) YukawaKFF_Array(i,j)=YukawaK_ff[i][j]; 
  PotentialsFF_List->set("PF5_ExpDecayParamFF", YukawaKFF_Array, "alpha - exponential decay parameter in fluid-fluid pair interactions");

  TwoDArray<double> NpowFF_Array(Ncomp,Ncomp);
  for (i=0; i<Ncomp;i++) for (j=0; j<Ncomp;j++) NpowFF_Array(i,j)=Npow_ff[i][j]; 
  PotentialsFF_List->set("PF6_NpowFF", NpowFF_Array, "n - variable power for r^-n repulsive term in Fluid-fluid pair interactions");

  TwoDArray<double> BondFF_Array(Ncomp,Ncomp);
  for (i=0; i<Ncomp;i++) for (j=0; j<Ncomp;j++) BondFF_Array(i,j)=Bond_ff[i][j]; 
  PotentialsFF_List->set("PF10_BondFF", BondFF_Array, "Bond Lengths for fluid-fluid pair interactions");

  Array<double> Mass_Array( Mass,Mass+Ncomp);
  PotentialsFF_List->set("PF7_Mass", Mass_Array, "Mass of each component");

  Array<double> Charge_Array( Charge_f,Charge_f+Ncomp);
  PotentialsFF_List->set("PF8_Charge", Charge_Array, "Charge of each component");

  Array<double> Polarization_Array( Pol,Pol+Ncomp);
  PotentialsFF_List->set("PF9_Polarization", Polarization_Array, "Polarizability of each component");

  return;
}
/**************************************************************************************************************/
void dft_GUI_potentialsFF_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                    Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                    Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                    Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                    Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List )
{

  RCP<StringVisualDependency> HSDiamType_Dep = rcp(
        new StringVisualDependency(Functional_List->getEntryRCP("F1_HS_Functional"),
                                    Fluid_List->getEntryRCP("F2_HSDiamType"), 
                                    tuple<std::string>("Ideal Fluid / No Volume Exclusion"), false));

  RCP<StringVisualDependency> PotType_Dep =rcp(
        new StringVisualDependency(Functional_List->getEntryRCP("F2_PAIRPOTcore_Functional"), 
                                   Fluid_List->getEntryRCP("F4_PairPotType"), 
                                   tuple<std::string>("No Mean Field Functional"), false));

  RCP<StringVisualDependency> HSDiam_Dep = rcp(
        new StringVisualDependency(Fluid_List->getEntryRCP("F2_HSDiamType"),
                                   Fluid_List->getEntryRCP("F3_HSDiam"), 
                                   tuple<std::string>("Manual Definition"), true));

  RCP<BoolVisualDependency> Temp_Dep = rcp(
        new BoolVisualDependency(Fluid_List->getEntryRCP("F5_Dimensionless_Energy_Entry"),
                                 Fluid_List->getEntryRCP("F6_Temperature"), false));

  RCP<NumberArrayLengthDependency<int,double> > HSDiam_Dep2 = rcp(
       new NumberArrayLengthDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"), 
                                        Fluid_List->getEntryRCP("F3_HSDiam")));


                                                   /* SET UP A VARIETY OF CONDITIONS THAT NEED TO BE COMBINED AS AND CONDITIONS */
  RCP<StringCondition> PotTypeStringCon = rcp(
        new StringCondition(PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"), "Manual Definition"));

  RCP<StringCondition> PotTypeStringConF = rcp(
        new StringCondition(PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"), "Lorentz-Berthlot Mixing"));

  RCP<StringCondition> SigmaStringCon1 = rcp(
     new StringCondition(Functional_List->getEntryRCP("F1_HS_Functional"),tuple<std::string>(
             "HS_Rosenfeld_Original","HS_Rosenfeld_Zerocrossover","HS_White_Bear","HS_White_Bear2")));

  RCP<StringCondition> ATTFuncStringCon1 = rcp(
     new StringCondition(Functional_List->getEntryRCP("F2_PAIRPOTcore_Functional"),tuple<std::string>(
             "UATT_core=[Umin; r=0,Rmin]","UATT_core=[0;r=0,Sigma]","UATT_core=[U(sigma);r=0,Sigma]",
             "UATT_core=[(0;r=0,Sigma),(Umin;r=Sigma,Rmin)]","UATTcore=[0;r=0,Rzero]")));

      /* set up conditions for making potential type selection */
  Condition::ConstConditionList PotEntry_conList=tuple<RCP<const Condition> >(SigmaStringCon1,ATTFuncStringCon1);
  RCP<OrCondition> PotEntry_allowed = rcp(new OrCondition(PotEntry_conList));
  RCP<ConditionVisualDependency> PotFFentry_Dep = rcp(
       new ConditionVisualDependency(PotEntry_allowed, PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"), true));

      /* set up SigmaFF conditions */
  Condition::ConstConditionList SigmaFF_conList=tuple<RCP<const Condition> >(SigmaStringCon1,PotTypeStringCon);
  RCP<AndCondition> SigmaFF_allowed = rcp(new AndCondition(SigmaFF_conList));
  RCP<ConditionVisualDependency> SigmaFF_Dep = rcp(
       new ConditionVisualDependency(SigmaFF_allowed, PotentialsFF_List->getEntryRCP("PF1_SigmaFF"), true));

      /* now set up SigmaF conditions */
  Condition::ConstConditionList SigmaF_conList=tuple<RCP<const Condition> >(SigmaStringCon1,PotTypeStringConF);
  RCP<AndCondition> SigmaF_allowed = rcp(new AndCondition(SigmaF_conList));
  RCP<ConditionVisualDependency> SigmaF_Dep = rcp(
       new ConditionVisualDependency(SigmaF_allowed, PotentialsFF_List->getEntryRCP("PF1_SigmaF"), true));

  RCP<StringCondition> EpsStringCon1 = rcp(
     new StringCondition(Fluid_List->getEntryRCP("F4_PairPotType"), 
                                 tuple<std::string>("LJ 12-6 potential (cut/shift)",
                                    "Exponential potential (cut/shift)","Square Well potential",
                                    "LJ 12-6 plus Yukawa potential (cut/shift)",
                                    "r^12 repulsion plus Yukawa potential (cut/shift)",
                                    "r^18 repulsion plus Yukawa potential (cut/shift)",
                                    "r^N repulsion plus Yukawa potential (cut/shift)")));

      /* set up EpsFF conditions */
  Condition::ConstConditionList EpsFF_conList=tuple<RCP<const Condition> >(ATTFuncStringCon1,EpsStringCon1,PotTypeStringCon);
  RCP<AndCondition> EpsFF_allowed = rcp(new AndCondition(EpsFF_conList));
  RCP<ConditionVisualDependency> EpsFF_Dep = rcp(
       new ConditionVisualDependency(EpsFF_allowed, PotentialsFF_List->getEntryRCP("PF2_EpsFF"), true));

      /* now set up EpsF conditions */
  Condition::ConstConditionList EpsF_conList=tuple<RCP<const Condition> >(ATTFuncStringCon1,EpsStringCon1,PotTypeStringConF);
  RCP<AndCondition> EpsF_allowed = rcp(new AndCondition(EpsF_conList));
  RCP<ConditionVisualDependency> EpsF_Dep = rcp(
       new ConditionVisualDependency(EpsF_allowed, PotentialsFF_List->getEntryRCP("PF2_EpsF"), true));

  RCP<StringCondition> CutStringCon1 = rcp(
      new StringCondition(Fluid_List->getEntryRCP("F4_PairPotType"),
                                 tuple<std::string>("LJ 12-6 potential (cut/shift)","Exponential potential (cut/shift)",
                                 "Coulomb potential as mean field (cut/shift)","Coulomb potential as mean field (cut only)",
                                 "Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
                                 "r^12 repulsion plus Yukawa potential (cut/shift)",
                                 "r^18 repulsion plus Yukawa potential (cut/shift)",
                                 "r^N repulsion plus Yukawa potential (cut/shift)")));

      /* set up CutFF conditions */
  Condition::ConstConditionList CutFF_conList=tuple<RCP<const Condition> >(ATTFuncStringCon1,CutStringCon1,PotTypeStringCon);
  RCP<AndCondition> CutFF_allowed = rcp(new AndCondition(CutFF_conList));
  RCP<ConditionVisualDependency> CutFF_Dep = rcp(
       new ConditionVisualDependency(CutFF_allowed, PotentialsFF_List->getEntryRCP("PF3_CutFF"), true));

      /* now set up CutF conditions */
  Condition::ConstConditionList CutF_conList=tuple<RCP<const Condition> >(ATTFuncStringCon1,CutStringCon1,PotTypeStringConF);
  RCP<AndCondition> CutF_allowed = rcp(new AndCondition(CutF_conList));
  RCP<ConditionVisualDependency> CutF_Dep = rcp(
       new ConditionVisualDependency(CutF_allowed, PotentialsFF_List->getEntryRCP("PF3_CutF"), true));


       /* HSdiam entry */
  RCP<StringVisualDependency> HSFuncVis_Dep = rcp(
      new StringVisualDependency(Functional_List->getEntryRCP("F1_HS_Functional"), Fluid_List->getEntryRCP("F2_HSDiamType"),
                                 "Ideal Fluid / No Volume Exclusion",false));

  RCP<StringCondition> EpsYukawaStringCon1 = rcp(
      new StringCondition(Fluid_List->getEntryRCP("F4_PairPotType"),
                                 tuple<std::string>("Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
                                 "r^12 repulsion plus Yukawa potential (cut/shift)", "r^18 repulsion plus Yukawa potential (cut/shift)",
                                 "r^N repulsion plus Yukawa potential (cut/shift)")));

      /* set up EpsYukawaFF conditions */
  Condition::ConstConditionList EpsYukawaFF_conList=tuple<RCP<const Condition> >(ATTFuncStringCon1,EpsYukawaStringCon1,PotTypeStringCon);
  RCP<AndCondition> EpsYukawaFF_allowed = rcp(new AndCondition(EpsYukawaFF_conList));
  RCP<ConditionVisualDependency> EpsYukawaFF_Dep = rcp(
       new ConditionVisualDependency(EpsYukawaFF_allowed, PotentialsFF_List->getEntryRCP("PF4_EpsYukawaFF"), true));

      /* now set up EpsYukawaF conditions */
  Condition::ConstConditionList EpsYukawaF_conList=tuple<RCP<const Condition> >(ATTFuncStringCon1,EpsYukawaStringCon1,PotTypeStringConF);
  RCP<AndCondition> EpsYukawaF_allowed = rcp(new AndCondition(EpsYukawaF_conList));
  RCP<ConditionVisualDependency> EpsYukawaF_Dep = rcp(
       new ConditionVisualDependency(EpsYukawaF_allowed, PotentialsFF_List->getEntryRCP("PF4_EpsYukawaF"), true));

  RCP<StringCondition> YukawaKStringCon1 = rcp(
      new StringCondition(Fluid_List->getEntryRCP("F4_PairPotType"),
                                 tuple<std::string>("Exponential potential (cut/shift)",
                                 "Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
                                 "r^12 repulsion plus Yukawa potential (cut/shift)",
                                 "r^18 repulsion plus Yukawa potential (cut/shift)",
                                 "r^N repulsion plus Yukawa potential (cut/shift)")));

      /* set up YukawaKFF conditions */
  Condition::ConstConditionList YukawaKFF_conList=tuple<RCP<const Condition> >(ATTFuncStringCon1,YukawaKStringCon1,PotTypeStringCon);
  RCP<AndCondition> YukawaKFF_allowed = rcp(new AndCondition(YukawaKFF_conList));
  RCP<ConditionVisualDependency> YukawaKFF_Dep = rcp(
       new ConditionVisualDependency(YukawaKFF_allowed, PotentialsFF_List->getEntryRCP("PF5_ExpDecayParamFF"), true));

      /* now set up YukawaKF conditions */
  Condition::ConstConditionList YukawaKF_conList=tuple<RCP<const Condition> >(ATTFuncStringCon1,YukawaKStringCon1,PotTypeStringConF);
  RCP<AndCondition> YukawaKF_allowed = rcp(new AndCondition(YukawaKF_conList));
  RCP<ConditionVisualDependency> YukawaKF_Dep = rcp(
       new ConditionVisualDependency(YukawaKF_allowed, PotentialsFF_List->getEntryRCP("PF5_ExpDecayParamF"), true));


   RCP<StringVisualDependency> NpowFFArray_Dep = rcp(
      new StringVisualDependency(Fluid_List->getEntryRCP("F4_PairPotType"), 
                                 PotentialsFF_List->getEntryRCP("PF6_NpowFF"), 
                                 tuple<std::string>("r^N repulsion plus Yukawa potential (cut/shift)"),true));

   RCP<StringVisualDependency> ChargeArray_Dep = rcp(
      new StringVisualDependency(Functional_List->getEntryRCP("F3_CHARGE_Functional"), 
                                 PotentialsFF_List->getEntryRCP("PF8_Charge"), 
                                 tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM",
                                 "Charge_with_DeltaC_General","Charge(MF)_with_Polarization"),true));

   RCP<StringVisualDependency> ChargeArray_Dep2 = rcp(
      new StringVisualDependency(Fluid_List->getEntryRCP("F4_PairPotType"),
                                 PotentialsFF_List->getEntryRCP("PF8_Charge"), 
                                 tuple<std::string>("Coulomb potential as mean field (cut/shift)",
                                  "Coulomb potential as mean field (cut only)"),true));

   RCP<StringVisualDependency> PolarizationArray_Dep = rcp(
      new StringVisualDependency(Functional_List->getEntryRCP("F3_CHARGE_Functional"),
                                 PotentialsFF_List->getEntryRCP("PF9_Polarization"), 
                                 tuple<std::string>("Charge(MF)_with_Polarization"),true));

   RCP<StringVisualDependency> BondArray_Dep = rcp(
      new StringVisualDependency(Functional_List->getEntryRCP("F4_POLYMER_Functional"),
                                 PotentialsFF_List->getEntryRCP("PF10_BondFF"), 
                                 tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT",
                                      "Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)",
                                      "Polymer_JDC_iSAFT(comp)"),true));
   
   Dependency::ParameterEntryList PotFF2DRowNumber_Dependents;
   PotFF2DRowNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF1_SigmaFF"));
   PotFF2DRowNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF2_EpsFF"));
   PotFF2DRowNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF3_CutFF"));
   PotFF2DRowNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF4_EpsYukawaFF"));
   PotFF2DRowNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF5_ExpDecayParamFF"));
   PotFF2DRowNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF6_NpowFF"));
   PotFF2DRowNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF10_BondFF"));
   RCP<TwoDRowDependency<int,double> > PotFF2DRowNumber_Dep = rcp(
             new TwoDRowDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"),PotFF2DRowNumber_Dependents));

   Dependency::ParameterEntryList PotFF2DColNumber_Dependents;
   PotFF2DColNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF1_SigmaFF"));
   PotFF2DColNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF2_EpsFF"));
   PotFF2DColNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF3_CutFF"));
   PotFF2DColNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF4_EpsYukawaFF"));
   PotFF2DColNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF5_ExpDecayParamFF"));
   PotFF2DColNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF6_NpowFF"));
   PotFF2DColNumber_Dependents.insert(PotentialsFF_List->getEntryRCP("PF10_BondFF"));
   RCP<TwoDColDependency<int,double> > PotFF2DColNumber_Dep = rcp(
             new TwoDColDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"),PotFF2DColNumber_Dependents));

   Dependency::ParameterEntryList PotFFArrayLength_Dependents;
   PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF7_Mass"));
   PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF8_Charge"));
   PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF9_Polarization"));
   PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF1_SigmaF"));
   PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF2_EpsF"));
   PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF3_CutF"));
   PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF4_EpsYukawaF"));
   PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF5_ExpDecayParamF"));
   RCP<NumberArrayLengthDependency<int,double> > PotFFArrayLength_Dep = rcp(
             new NumberArrayLengthDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"), PotFFArrayLength_Dependents));


 /*****************************************/
 /* add the dependencies for this section.*/
 /*****************************************/
    depSheet_Tramonto->addDependency(PotType_Dep);
    depSheet_Tramonto->addDependency(HSDiamType_Dep);
    depSheet_Tramonto->addDependency(HSDiam_Dep);
    depSheet_Tramonto->addDependency(HSDiam_Dep2);
    depSheet_Tramonto->addDependency(Temp_Dep);
    depSheet_Tramonto->addDependency(EpsFF_Dep);
    depSheet_Tramonto->addDependency(EpsF_Dep);
    depSheet_Tramonto->addDependency(CutF_Dep);
    depSheet_Tramonto->addDependency(CutFF_Dep);
    depSheet_Tramonto->addDependency(EpsYukawaFF_Dep);
    depSheet_Tramonto->addDependency(EpsYukawaF_Dep);
    depSheet_Tramonto->addDependency(YukawaKF_Dep);
    depSheet_Tramonto->addDependency(YukawaKFF_Dep);
    depSheet_Tramonto->addDependency(NpowFFArray_Dep);
    depSheet_Tramonto->addDependency(ChargeArray_Dep);
    depSheet_Tramonto->addDependency(ChargeArray_Dep2);
    depSheet_Tramonto->addDependency(PolarizationArray_Dep);
    depSheet_Tramonto->addDependency(BondArray_Dep);
    depSheet_Tramonto->addDependency(HSFuncVis_Dep);
    depSheet_Tramonto->addDependency(SigmaF_Dep);
    depSheet_Tramonto->addDependency(SigmaFF_Dep);
    depSheet_Tramonto->addDependency(PotFFentry_Dep);
    depSheet_Tramonto->addDependency(PotFFArrayLength_Dep);
    depSheet_Tramonto->addDependency(PotFF2DRowNumber_Dep);
    depSheet_Tramonto->addDependency(PotFF2DColNumber_Dep);

    return;
}
