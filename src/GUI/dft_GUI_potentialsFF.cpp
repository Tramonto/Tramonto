using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;
bool func_testNcomp(Teuchos::RCP<Teuchos::ParameterList> Fluid_List);

void dft_GUI_potentialsFF(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List )
{
    /****************************************************************************************************************/
  /****************************** FUNCTIONAL CONTROL PARAMETER SECTION ********************************************/
  /****************************************************************************************************************/
    /**************************************/
    /* Define validators for this section.*/
    /**************************************/

    RCP<EnhancedNumberValidator<int> > NcompValidator = rcp(new EnhancedNumberValidator<int>(1,1000,1));

    RCP<StringValidator> HSDiamTypeValidator = rcp(
           new StringValidator(tuple<std::string>("SIGMA","Manual Definition","BARKER_HENDERSON")));

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
    /* set up independent parameters that do not have dependencies */
    /***************************************************************/
    Fluid_List->set("F1_Ncomp", 1, "Number of fluid components in the problem \n(Note: For co-polymers, each bead type is a different component)", NcompValidator);
    

    /*******************************/
    /* define dependent parameters */
    /*******************************/

    
    Fluid_List->set("F2_HSDiamType", "SIGMA", "Select method for defining hard sphere diameters", HSDiamTypeValidator);

    Array<double> HSDiam_Array( (Fluid_List->get<int>("F1_Ncomp")),1.0);
    Fluid_List->set("F3_HSDiam", HSDiam_Array, "Set hard sphere diameters manually: Ncomp entries");

    Fluid_List->set("F4_PairPotType", "none", "Type of pair interaction potential", PairPotValidator);

    PotentialsFF_List->set("PF0_Off_Diagonal_Definitions", "Manual Definition", "Options for definition of off-diagonal pair potential parameters", MixTypeValidator);

           /* items for PotentialsFF_sublist */
    Array<double> SigmaFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),1.0);
    PotentialsFF_List->set("PF1_SigmaFF", SigmaFF_Array, "Sigma - Characteristic diameter fluid-fluid pair interactions");

    Array<double> EpsFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),1.0);
    PotentialsFF_List->set("PF2_EpsFF", EpsFF_Array, "Eps - Energy prefactor for fluid-fluid pair interactions");

    Array<double> CutFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),3.0);
    PotentialsFF_List->set("PF3_CutFF", CutFF_Array, "rcut - Cutoff distance for fluid-fluid pair interactions");

    Array<double> EpsYukawaFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),1.0);
    PotentialsFF_List->set("PF4_EpsYukawaFF", EpsYukawaFF_Array, "AYuk - Energy prefactor for Yukawa term in fluid-fluid pair interactions");

    Array<double> YukawaKFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),1.0);
    PotentialsFF_List->set("PF5_ExpDecayParamFF", YukawaKFF_Array, "alpha - exponential parameter in Yukawa term of fluid-fluid pair interactions");

    Array<double> NpowFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),12.0);
    PotentialsFF_List->set("PF6_NpowFF", NpowFF_Array, "n - variable power for r^-n repulsive term in Fluid-fluid pair interactions");

    Array<double> Mass_Array( Fluid_List->get<int>("F1_Ncomp"),1.0);
    PotentialsFF_List->set("PF7_Mass", Mass_Array, "Mass of each component");

    Array<double> Charge_Array( Fluid_List->get<int>("F1_Ncomp"),0.0);
    PotentialsFF_List->set("PF8_Charge", Charge_Array, "Charge of each component");

    Array<double> Polarization_Array( Fluid_List->get<int>("F1_Ncomp"),0.0);
    PotentialsFF_List->set("PF9_Polarization", Polarization_Array, "Polarizability of each component");

/*    Array<double> BondFF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Fluid_List->get<int>("F1_Ncomp")),0.0);
    PotentialsFF_List->set("P9_BondFF", BondFF_Array, "Bond Lengths for fluid-fluid pair interactions");*/

/*    Array<double> VextMembrane_Array( (Surface_List->get<int>("S3_Nsurface_types"))*(Fluid_List->get<int>("F1_Ncomp")),SurfaceInteraction_List->get<double>("SI1_VEXT_MAX"));*/
/*    PotentialsFF_List->set("PW1_Vext_membrane", VextMembrane_Array, "constant Vext values inside semipermiable surfaces.  Set to VEXT_MAX to make surface impenetrable.");*/


    /***************************************************************************************************************/
    /* define dependencies (show dependent parameters only if the independent parameters has a particular setting. */
    /***************************************************************************************************************/

     RCP<StringVisualDependency> PotType_Dep =rcp(
           new StringVisualDependency("F2_PAIRPOTcore_Functional", Functional_List, "F4_PairPotType", Fluid_List, tuple<std::string>("No Mean Field Functional"), false));

      /* trying to set up dependency for a sublist based on the value of F4_PairPotType .... doesn't work */
/*     RCP<StringVisualDependency> PairPotList_Dep =rcp(
           new StringVisualDependency("F4_PairPotType", Fluid_List, "PotentialsFF_List", Fluid_List, "NO_PAIRPOT", false));*/

     RCP<StringVisualDependency> HSDiamType_Dep = rcp(
           new StringVisualDependency( "F1_HS_Functional", Functional_List, "F2_HSDiamType", Fluid_List, tuple<std::string>("Ideal Fluid / No Volume Exclusion"), false));
     RCP<StringVisualDependency> HSDiam_Dep = rcp(
           new StringVisualDependency( "F2_HSDiamType", Fluid_List, "F3_HSDiam", Fluid_List, tuple<std::string>("Manual Definition"), true));
    RCP<NumberArrayLengthDependency> HSDiam_Dep2 = rcp(
          new NumberArrayLengthDependency( "F1_Ncomp", Fluid_List, "F3_HSDiam",Fluid_List));


      /* trying to set up dependency based on the value of F1_Ncomp-1 rather than F1_Ncomp .... doesn't work */
/*     RCP<NumberVisualDependency<int> > MixType_Dep = rcp(
            new NumberVisualDependency<int>( "F1_Ncomp", Fluid_List, "PF0_Off_Diagonal_Definitions", PotentialsFF_List,func_testNcomp(Fluid_List)));*/

     RCP<NumberVisualDependency<int> > MixType_Dep = rcp(
           new NumberVisualDependency<int>( "F1_Ncomp", Fluid_List,"PF0_Off_Diagonal_Definitions", PotentialsFF_List));

     RCP<StringVisualDependency> MixType2_Dep =rcp(
           new StringVisualDependency("F2_PAIRPOTcore_Functional", Functional_List, "PF0_Off_Diagonal_Definitions", PotentialsFF_List, 
           tuple<std::string>("No Mean Field Functional"), false));

      RCP<StringVisualDependency> SigmaFFArray_Dep2 = rcp(
         new StringVisualDependency("F1_HS_Functional", Functional_List, "PF1_SigmaFF", PotentialsFF_List, "Ideal Fluid / No Volume Exclusion",false));

      RCP<StringVisualDependency> EpsFFArray_Dep2 = rcp(
         new StringVisualDependency("F4_PairPotType", Fluid_List, "PF2_EpsFF", PotentialsFF_List, 
             tuple<std::string>("LJ 12-6 potential (cut/shift)",
                "Exponential potential (cut/shift)","Square Well potential",
                "LJ 12-6 plus Yukawa potential (cut/shift)",
                "r^12 repulsion plus Yukawa potential (cut/shift)",
                "r^18 repulsion plus Yukawa potential (cut/shift)",
		"r^N repulsion plus Yukawa potential (cut/shift)"),true));


      RCP<StringVisualDependency> CutFFArray_Dep2 = rcp(
         new StringVisualDependency("F4_PairPotType", Fluid_List, "PF3_CutFF", PotentialsFF_List, 
             tuple<std::string>("LJ 12-6 potential (cut/shift)","Exponential potential (cut/shift)",
		"Coulomb potential as mean field (cut/shift)","Coulomb potential as mean field (cut only)",
		"Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
		"r^12 repulsion plus Yukawa potential (cut/shift)",
		"r^18 repulsion plus Yukawa potential (cut/shift)",
		"r^N repulsion plus Yukawa potential (cut/shift)"),true));


      RCP<StringVisualDependency> EpsYukawaFFArray_Dep2 = rcp(
         new StringVisualDependency("F4_PairPotType", Fluid_List, "PF4_EpsYukawaFF", PotentialsFF_List, 
             tuple<std::string>("Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
		"r^12 repulsion plus Yukawa potential (cut/shift)",
		"r^18 repulsion plus Yukawa potential (cut/shift)",
		"r^N repulsion plus Yukawa potential (cut/shift)"),true));

      RCP<StringVisualDependency> YukawaKFFArray_Dep2 = rcp(
         new StringVisualDependency("F4_PairPotType", Fluid_List, "PF5_ExpDecayParamFF", PotentialsFF_List, 
             tuple<std::string>("Exponential potential (cut/shift)",
 		"Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
		"r^12 repulsion plus Yukawa potential (cut/shift)",
		"r^18 repulsion plus Yukawa potential (cut/shift)",
		"r^N repulsion plus Yukawa potential (cut/shift)"),true));

      RCP<StringVisualDependency> NpowFFArray_Dep2 = rcp(
         new StringVisualDependency("F4_PairPotType", Fluid_List, "PF6_NpowFF", PotentialsFF_List, 
             tuple<std::string>("r^N repulsion plus Yukawa potential (cut/shift)"),true));

      RCP<StringVisualDependency> ChargeArray_Dep2 = rcp(
         new StringVisualDependency("F3_CHARGE_Functional", Functional_List, "PF8_Charge", PotentialsFF_List, 
             tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM","Charge_with_DeltaC_General","Charge(MF)_with_Polarization"),true));

      RCP<StringVisualDependency> ChargeArray_Dep3 = rcp(
         new StringVisualDependency("F4_PairPotType", Fluid_List, "PF8_Charge", PotentialsFF_List, 
             tuple<std::string>("Coulomb potential as mean field (cut/shift)","Coulomb potential as mean field (cut only)"),true));

      RCP<StringVisualDependency> PolarizationArray_Dep2 = rcp(
         new StringVisualDependency("F3_CHARGE_Functional", Functional_List, "PF9_Polarization", PotentialsFF_List, 
              tuple<std::string>("Charge(MF)_with_Polarization"),true));


      Dependency::ParameterParentMap PotFFArrayLength_Dependents;
      PotFFArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF1_SigmaFF", PotentialsFF_List));
      PotFFArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF2_EpsFF", PotentialsFF_List));
      PotFFArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF3_CutFF", PotentialsFF_List));
      PotFFArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF4_EpsYukawaFF", PotentialsFF_List));
      PotFFArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF5_ExpDecayParamFF", PotentialsFF_List));
      PotFFArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF6_NpowFF", PotentialsFF_List));

      RCP<NumberArrayLengthDependency> PotFFArrayLength_Dep = rcp(
                new NumberArrayLengthDependency( "F1_Ncomp", Fluid_List, PotFFArrayLength_Dependents));

      Dependency::ParameterParentMap PotFF2ArrayLength_Dependents;
      PotFF2ArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF7_Mass", PotentialsFF_List));
      PotFF2ArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF8_Charge", PotentialsFF_List));
      PotFF2ArrayLength_Dependents.insert(std::pair<std::string, RCP<ParameterList> >("PF9_Polarization", PotentialsFF_List));

      RCP<NumberArrayLengthDependency> PotFFArrayLength2_Dep = rcp(
                new NumberArrayLengthDependency( "F1_Ncomp", Fluid_List, PotFF2ArrayLength_Dependents));

/*      RCP<BoolVisualDependency> VextSemiperm_Dep = rcp(
           new BoolVisualDependency( "SI5_LSemiperm", SurfaceInteraction_List, "PW1_Vext_membrane", PotentialsFF_List,true));*/


/* need to figure out how to do 2D arrays ---- if 1D array of NxN length, how do we make the NumberArrayLengthDependency work?*/

    /*****************************************/
    /* add the dependencies for this section.*/
    /*****************************************/
       depSheet_Tramonto->addDependency(PotType_Dep);
       depSheet_Tramonto->addDependency(HSDiamType_Dep);
       depSheet_Tramonto->addDependency(HSDiam_Dep);
       depSheet_Tramonto->addDependency(HSDiam_Dep2);
       depSheet_Tramonto->addDependency(MixType_Dep);
       depSheet_Tramonto->addDependency(MixType2_Dep);
       depSheet_Tramonto->addDependency(SigmaFFArray_Dep2);
       depSheet_Tramonto->addDependency(EpsFFArray_Dep2);
       depSheet_Tramonto->addDependency(CutFFArray_Dep2);
       depSheet_Tramonto->addDependency(EpsYukawaFFArray_Dep2);
       depSheet_Tramonto->addDependency(YukawaKFFArray_Dep2);
       depSheet_Tramonto->addDependency(NpowFFArray_Dep2);
       depSheet_Tramonto->addDependency(ChargeArray_Dep2);
       depSheet_Tramonto->addDependency(ChargeArray_Dep3);
       depSheet_Tramonto->addDependency(PolarizationArray_Dep2);

       depSheet_Tramonto->addDependency(PotFFArrayLength_Dep);
       depSheet_Tramonto->addDependency(PotFFArrayLength2_Dep);
     
      
/*       depSheet_Tramonto->addDependency(VextSemiperm_Dep);*/
/*       depSheet_Tramonto->addDependency(PairPotList_Dep);*/

  /****************************************************************************************************************/
  /****************************** END FUNCTIONAL CONTROL PARAMETER SECTION ****************************************/
  /****************************************************************************************************************/
 
  /****************************************************************************************************************/
  /****************************** ANY NEW SUBSECTION PARAMETER SECTION ********************************************/
  /****************************************************************************************************************/
    /* SUBLIST(S) */
    /* DEFINE VALIDATORS*/ 
    /* INDEPENDENT PARAMETERS */
    /* DEPENDENT PARAMETERS */
    /* DEPENDENCIES */
    /* DEPENDENCY SHEET ENTRIES*/
  /****************************************************************************************************************/
  /****************************** END ANY NEW SUBSECTION PARAMETER SECTION ****************************************/
  /****************************************************************************************************************/

  return;
}
/**************************************************************************************************************/
bool func_testNcomp(Teuchos::RCP<Teuchos::ParameterList> Fluid_List)
{
  int test;
  bool logical;
   test=Fluid_List->get<int>("F1_Ncomp")-1;
   if (test>0) logical=true;
   else logical=false;
   return logical;
}

