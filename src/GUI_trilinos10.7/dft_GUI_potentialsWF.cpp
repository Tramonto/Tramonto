using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_potentialsWF(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Surface_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List )
{
    /****************************************************************************************************************/
  /****************************** FUNCTIONAL CONTROL PARAMETER SECTION ********************************************/
  /****************************************************************************************************************/
    /**************************************/
    /* Define validators for this section.*/
    /**************************************/

    RCP<StringValidator> HSDiamValidator = rcp(
           new StringValidator(tuple<std::string>("SIGMA","MANUAL","BARKER_HENDERSON")));

    RCP<StringValidator> MixTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Lorentz-Berthlot Mixing","Manual Definition")));

    RCP<StringValidator> PairPotValidator = rcp(
           new StringValidator(tuple<std::string>("NO_PAIRPOT","LJ12_6_CS","COULOMB_CS","COULOMB_noShift",
                          "YUKAWA_CS","EXP_CS","SQUARE_WELL","LJandYUKAWA_CS","r12andYUKAWA_CS","r18andYUKAWA_CS","rNandYUKAWA_CS")));


    /***************************************************************/
    /* set up independent parameters that do not have dependencies */
    /***************************************************************/
    

    /*******************************/
    /* define dependent parameters */
    /*******************************/

           /* items for PotentialsWF_sublist */
    Array<double> SigmaWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),1.0);
    PotentialsWF_List->set("SP1_SigmaWF", SigmaWF_Array, "Sigma - Characteristic diameter fluid-fluid pair interactions");

    Array<double> EpsWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),1.0);
    PotentialsWF_List->set("SP2_EpsWF", EpsWF_Array, "Eps - Energy prefactor for fluid-fluid pair interactions");

    Array<double> CutWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),3.0);
    PotentialsWF_List->set("SP3_CutWF", CutWF_Array, "rcut - Cutoff distance for fluid-fluid pair interactions");

    Array<double> EpsYukawaWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),1.0);
    PotentialsWF_List->set("SP4_EpsYukawaWF", EpsYukawaWF_Array, "AYuk - Energy prefactor for Yukawa term in fluid-fluid pair interactions");

    Array<double> YukawaKWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),1.0);
    PotentialsWF_List->set("SP5_ExpDecayParamWF", YukawaKWF_Array, "alpha - exponential parameter in Yukawa term of fluid-fluid pair interactions");

    Array<double> VextMembrane_Array( (Surface_List->get<int>("S3: Number of surface types"))*(Surface_List->get<int>("S3: Number of surface types")),SurfaceInteraction_List->get<double>("SI1_VEXT_MAX"));
    PotentialsWF_List->set("PW1_Vext_membrane", VextMembrane_Array, "constant Vext values inside semipermiable surfaces.  Set to VEXT_MAX to make surface impenetrable.");


    /***************************************************************************************************************/
    /* define dependencies (show dependent parameters only if the independent parameters has a particular setting. */
    /***************************************************************************************************************/


               /* Visual Dependencies */

      RCP<StringVisualDependency> EpsWFArray_Dep2 = rcp(
         new StringVisualDependency(
             SurfaceInteraction_List->getEntryRCP("SI4_U3D_TO_BE_INTEGRATED"), 
             PotentialsWF_List->getEntryRCP("SP2_EpsWF"), 
             tuple<std::string>("LJ 12-6 potential (cut/shift)",
                "Exponential potential (cut/shift)","Square Well potential",
                "LJ 12-6 plus Yukawa potential (cut/shift)",
                "r^12 repulsion plus Yukawa potential (cut/shift)",
                "r^18 repulsion plus Yukawa potential (cut/shift)",
                "r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> EpsWFArray_Dep3 = rcp(
         new StringVisualDependency(
             SurfaceInteraction_List->getEntryRCP("SI2_Vext_1D"), 
             PotentialsWF_List->getEntryRCP("SP2_EpsWF"),
             tuple<std::string>("LJ 9-3 potential (cut/shift)","LJ 9-3 potential (version 2) (cut/shift)",
				"LJ 9-3 potential (no c/s)","LJ 9-3 potential with shifted x",
				"r^9 repulsive potential (no c/s)","Expotential potential (no c/s)",
				"R^7 repulsion plus Yukawa field (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> CutWFArray_Dep2 = rcp(
         new StringVisualDependency(
             SurfaceInteraction_List->getEntryRCP("SI4_U3D_TO_BE_INTEGRATED"), 
             PotentialsWF_List->getEntryRCP("SP3_CutWF"), 
             tuple<std::string>("LJ 12-6 potential (cut/shift)","Exponential potential (cut/shift)",
                "Coulomb potential as mean field (cut/shift)","Coulomb potential as mean field (cut only)",
                "Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
                "r^12 repulsion plus Yukawa potential (cut/shift)",
                "r^18 repulsion plus Yukawa potential (cut/shift)",
                "r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> CutWFArray_Dep3 = rcp(
         new StringVisualDependency(
             SurfaceInteraction_List->getEntryRCP("SI2_Vext_1D"), 
             PotentialsWF_List->getEntryRCP("SP3_CutWF"),
             tuple<std::string>("LJ 9-3 potential (cut/shift)","LJ 9-3 potential (version 2) (cut/shift)",
				"LJ 9-3 potential with shifted x",
				"R^7 repulsion plus Yukawa field (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> EpsYukawaWFArray_Dep2 = rcp(
         new StringVisualDependency(
             SurfaceInteraction_List->getEntryRCP("SI4_U3D_TO_BE_INTEGRATED"), 
             PotentialsWF_List->getEntryRCP("SP4_EpsYukawaWF"), 
             tuple<std::string>("Yukawa potential (cut/shift)",
				"LJ 12-6 plus Yukawa potential (cut/shift)",
		                "r^12 repulsion plus Yukawa potential (cut/shift)",
       			        "r^18 repulsion plus Yukawa potential (cut/shift)",
				"r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> EpsYukawaWFArray_Dep3 = rcp(
         new StringVisualDependency(
             SurfaceInteraction_List->getEntryRCP("SI2_Vext_1D"), 
             PotentialsWF_List->getEntryRCP("SP4_EpsYukawaWF"), 
             tuple<std::string>("R^7 repulsion plus Yukawa field (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> YukawaKWFArray_Dep2 = rcp(
         new StringVisualDependency(
             SurfaceInteraction_List->getEntryRCP("SI4_U3D_TO_BE_INTEGRATED"), 
             PotentialsWF_List->getEntryRCP("SP5_ExpDecayParamWF"),  
             tuple<std::string>("Exponential potential (cut/shift)",
				"Yukawa potential (cut/shift)",
				"LJ 12-6 plus Yukawa potential (cut/shift)",
		                "r^12 repulsion plus Yukawa potential (cut/shift)",
       			        "r^18 repulsion plus Yukawa potential (cut/shift)",
				"r^N repulsion plus Yukawa potential (cut/shift)"),
             true)
      );

      RCP<StringVisualDependency> YukawaKWFArray_Dep3 = rcp(
         new StringVisualDependency(
             SurfaceInteraction_List->getEntryRCP("SI2_Vext_1D"),
             PotentialsWF_List->getEntryRCP("SP5_ExpDecayParamWF"),
             tuple<std::string>("Expotential potential (no c/s)",
				"R^7 repulsion plus Yukawa field (cut/shift)"),
             true)
      );


             /* Array Length Dependencies */
      RCP<NumberArrayLengthDependency<int, double> > SigmaWFArray_Dep = rcp(
         new NumberArrayLengthDependency<int, double> ( 
             Fluid_List->getEntryRCP("F1_Ncomp"), 
             PotentialsWF_List->getEntryRCP("SP1_SigmaWF"))
      );

      RCP<NumberArrayLengthDependency<int, double> > SigmaWFArray_Dep2 = rcp(
         new NumberArrayLengthDependency<int, double> ( 
             Surface_List->getEntryRCP("S3: Number of surface types"), 
             PotentialsWF_List->getEntryRCP("SP1_SigmaWF"))
      );

      RCP<NumberArrayLengthDependency<int, double> > EpsWFArray_Dep = rcp(
         new NumberArrayLengthDependency<int, double> ( 
             Fluid_List->getEntryRCP("F1_Ncomp"), 
             PotentialsWF_List->getEntryRCP("SP2_EpsWF"))
      );

      RCP<NumberArrayLengthDependency<int, double> > CutWFArray_Dep = rcp(
         new NumberArrayLengthDependency<int, double> (
             Fluid_List->getEntryRCP("F1_Ncomp"),
             PotentialsWF_List->getEntryRCP("SP3_CutWF"))
      );

      RCP<NumberArrayLengthDependency<int, double> > EpsYukawaWFArray_Dep = rcp(
         new NumberArrayLengthDependency<int, double> ( 
             Fluid_List->getEntryRCP("F1_Ncomp"), 
             PotentialsWF_List->getEntryRCP("SP4_EpsYukawaWF"))
      );

      RCP<NumberArrayLengthDependency<int, double> > YukawaKWFArray_Dep = rcp(
         new NumberArrayLengthDependency<int, double> ( 
             Fluid_List->getEntryRCP("F1_Ncomp"), 
             PotentialsWF_List->getEntryRCP("SP5_ExpDecayParamWF"))
      );

      RCP<BoolVisualDependency> VextSemiperm_Dep = rcp(
         new BoolVisualDependency( 
             SurfaceInteraction_List->getEntryRCP("SI5_LSemiperm"), 
             PotentialsWF_List->getEntryRCP("PW1_Vext_membrane"), 
             true)
      );


/* need to figure out how to do 2D arrays ---- if 1D array of NxN length, how do we make the NumberArrayLengthDependency work?*/

    /*****************************************/
    /* add the dependencies for this section.*/
    /*****************************************/
       depSheet_Tramonto->addDependency(SigmaWFArray_Dep2);
       depSheet_Tramonto->addDependency(SigmaWFArray_Dep);

       depSheet_Tramonto->addDependency(EpsWFArray_Dep3);
       depSheet_Tramonto->addDependency(EpsWFArray_Dep2);
       depSheet_Tramonto->addDependency(EpsWFArray_Dep);

       depSheet_Tramonto->addDependency(CutWFArray_Dep3);
       depSheet_Tramonto->addDependency(CutWFArray_Dep2);
       depSheet_Tramonto->addDependency(CutWFArray_Dep);

       depSheet_Tramonto->addDependency(EpsYukawaWFArray_Dep3);
       depSheet_Tramonto->addDependency(EpsYukawaWFArray_Dep2);
       depSheet_Tramonto->addDependency(EpsYukawaWFArray_Dep);

       depSheet_Tramonto->addDependency(YukawaKWFArray_Dep3);
       depSheet_Tramonto->addDependency(YukawaKWFArray_Dep2);
       depSheet_Tramonto->addDependency(YukawaKWFArray_Dep);

       depSheet_Tramonto->addDependency(VextSemiperm_Dep);

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
