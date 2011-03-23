using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;
int func_testNcomp(int);

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

    Fluid_List->set("F5_Dimensionless_Energy_Entry", true, "Dimensionless energies must be provided in units of epsilon/kT.\n Otherwise, set to false for entry as epsilon/k (units of K)");

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
           new StringVisualDependency(
              Functional_List->getEntryRCP("F2_PAIRPOTcore_Functional"), 
              Fluid_List->getEntryRCP("F4_PairPotType"), 
              tuple<std::string>("No Mean Field Functional"), 
              false)
     );

      /* trying to set up dependency for a sublist based on the value of F4_PairPotType .... doesn't work */
/*     RCP<StringVisualDependency> PairPotList_Dep =rcp(
           new StringVisualDependency("F4_PairPotType", Fluid_List, "PotentialsFF_List", Fluid_List, "NO_PAIRPOT", false));*/

     RCP<StringVisualDependency> HSDiamType_Dep = rcp(
           new StringVisualDependency( 
               Functional_List->getEntryRCP("F1_HS_Functional"),
               Fluid_List->getEntryRCP("F2_HSDiamType"), 
               tuple<std::string>("Ideal Fluid / No Volume Exclusion"), 
               false)
     );

     RCP<StringVisualDependency> HSDiam_Dep = rcp(
           new StringVisualDependency( 
               Fluid_List->getEntryRCP("F2_HSDiamType"), 
               Fluid_List->getEntryRCP("F3_HSDiam"), 
               tuple<std::string>("Manual Definition"), 
               true)
     );

    RCP<NumberArrayLengthDependency<int, double> > HSDiam_Dep2 = rcp(
          new NumberArrayLengthDependency<int, double> (
               Fluid_List->getEntryRCP("F1_Ncomp"),
               Fluid_List->getEntryRCP("F3_HSDiam"))
    );


      /* trying to set up dependency based on the value of F1_Ncomp-1 rather than F1_Ncomp .... doesn't work */
/*     RCP<NumberVisualDependency<int> > MixType_Dep = rcp(
            new NumberVisualDependency<int>(
               Fluid_List->getEntryRCP("F1_Ncomp"), 
               PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"),
               func_testNcomp)
     );*/

     RCP<NumberVisualDependency<int> > MixType_Dep = rcp(
        new NumberVisualDependency<int>( 
            Fluid_List->getEntryRCP("F1_Ncomp"), 
            PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"))
     );

     RCP<StringVisualDependency> MixType2_Dep =rcp(
        new StringVisualDependency(
            Functional_List->getEntryRCP("F2_PAIRPOTcore_Functional"), 
            PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"),  
            tuple<std::string>("No Mean Field Functional"), 
            false)
      );

      RCP<StringVisualDependency> SigmaFFArray_Dep2 = rcp(
         new StringVisualDependency(
             Functional_List->getEntryRCP("F1_HS_Functional"), 
             PotentialsFF_List->getEntryRCP("PF1_SigmaFF"), 
             "Ideal Fluid / No Volume Exclusion",
             false)
      );

      RCP<StringVisualDependency> EpsFFArray_Dep2 = rcp(
         new StringVisualDependency(
             Fluid_List->getEntryRCP("F4_PairPotType"), 
             PotentialsFF_List->getEntryRCP("PF2_EpsFF"),  
             tuple<std::string>("LJ 12-6 potential (cut/shift)",
                "Exponential potential (cut/shift)","Square Well potential",
                "LJ 12-6 plus Yukawa potential (cut/shift)",
                "r^12 repulsion plus Yukawa potential (cut/shift)",
                "r^18 repulsion plus Yukawa potential (cut/shift)",
		"r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
     );


      RCP<StringVisualDependency> CutFFArray_Dep2 = rcp(
         new StringVisualDependency(
             Fluid_List->getEntryRCP("F4_PairPotType"), 
             PotentialsFF_List->getEntryRCP("PF3_CutFF"), 
             tuple<std::string>("LJ 12-6 potential (cut/shift)","Exponential potential (cut/shift)",
		"Coulomb potential as mean field (cut/shift)","Coulomb potential as mean field (cut only)",
		"Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
		"r^12 repulsion plus Yukawa potential (cut/shift)",
		"r^18 repulsion plus Yukawa potential (cut/shift)",
		"r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
      );


      RCP<StringVisualDependency> EpsYukawaFFArray_Dep2 = rcp(
         new StringVisualDependency(
             Fluid_List->getEntryRCP("F4_PairPotType"),
             PotentialsFF_List->getEntryRCP("PF4_EpsYukawaFF"), 
             tuple<std::string>("Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
		"r^12 repulsion plus Yukawa potential (cut/shift)",
		"r^18 repulsion plus Yukawa potential (cut/shift)",
		"r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> YukawaKFFArray_Dep2 = rcp(
         new StringVisualDependency(
             Fluid_List->getEntryRCP("F4_PairPotType"), 
             PotentialsFF_List->getEntryRCP("PF5_ExpDecayParamFF"),  
             tuple<std::string>("Exponential potential (cut/shift)",
 		"Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
		"r^12 repulsion plus Yukawa potential (cut/shift)",
		"r^18 repulsion plus Yukawa potential (cut/shift)",
		"r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> NpowFFArray_Dep2 = rcp(
         new StringVisualDependency(
             Fluid_List->getEntryRCP("F4_PairPotType"), 
             PotentialsFF_List->getEntryRCP("PF6_NpowFF"), 
             tuple<std::string>("r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> ChargeArray_Dep2 = rcp(
         new StringVisualDependency(
             Functional_List->getEntryRCP("F3_CHARGE_Functional"),
             PotentialsFF_List->getEntryRCP("PF8_Charge"), 
             tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM",
                  "Charge_with_DeltaC_General","Charge(MF)_with_Polarization"),
             true)
      );

      RCP<StringVisualDependency> ChargeArray_Dep3 = rcp(
         new StringVisualDependency(
             Fluid_List->getEntryRCP("F4_PairPotType"), 
             PotentialsFF_List->getEntryRCP("PF8_Charge"),  
             tuple<std::string>("Coulomb potential as mean field (cut/shift)",
                                "Coulomb potential as mean field (cut only)"),
             true)
      );

      RCP<StringVisualDependency> PolarizationArray_Dep2 = rcp(
         new StringVisualDependency(
             Functional_List->getEntryRCP("F3_CHARGE_Functional"), 
             PotentialsFF_List->getEntryRCP("PF9_Polarization"),  
             tuple<std::string>("Charge(MF)_with_Polarization"),
             true)
      );


      Dependency::ParameterEntryList PotFFArrayLength_Dependents;
      PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF1_SigmaFF"));
      PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF2_EpsFF"));
      PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF3_CutFF"));
      PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF4_EpsYukawaFF"));
      PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF5_ExpDecayParamFF"));
      PotFFArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF6_NpowFF"));

      RCP<NumberArrayLengthDependency<int, double> > PotFFArrayLength_Dep = rcp(
         new NumberArrayLengthDependency<int, double> ( 
             Fluid_List->getEntryRCP("F1_Ncomp"), 
             PotFFArrayLength_Dependents)
      );

      Dependency::ParameterEntryList PotFF2ArrayLength_Dependents;
      PotFF2ArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF7_Mass"));
      PotFF2ArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF8_Charge"));
      PotFF2ArrayLength_Dependents.insert(PotentialsFF_List->getEntryRCP("PF9_Polarization"));

      RCP<NumberArrayLengthDependency<int, double> > PotFFArrayLength2_Dep = rcp(
         new NumberArrayLengthDependency<int, double> ( 
             Fluid_List->getEntryRCP("F1_Ncomp"),
             PotFF2ArrayLength_Dependents)
      );

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
int func_testNcomp(int Ncomp)
{
   return Ncomp-1;
}

