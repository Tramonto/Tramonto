using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_vextType(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List)
{

  /****************************************************************************************************************/
  /********************************* SURFACE DEFINITIONS SECTION **************************************************/
  /****************************************************************************************************************/

     /* VALIDATORS*/

        /* 1D options */
        RCP<StringValidator> Vext1Dtype_Valid = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "None",
                       "hard surface",
                       "Vext 1D",
                       "Numerically integrated Vext from Uij")));

        /* 2D options */
        RCP<StringValidator> Vext2Dtype_Valid = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "None",
                       "hard surface",
                       "Vext 1D_XMIN",
                       "Vext 1D_ORIENTATION",
                       "Numerically integrated Vext from Uij")));
        /* 3D options */
        RCP<StringValidator> Vext3Dtype_Valid = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "None",
                       "hard surface",
                       "Vext 1D_XMIN",
                       "Vext 1D_ORIENTATION",
                       "Numerically integrated Vext from Uij",
                       "Vext for atomic surfaces")));

        RCP<StringValidator> U3D_Valid = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "none",
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

        RCP<StringValidator> Vext1D_Valid = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "none",
                       "LJ 9-3 potential (cut/shift)",
                       "LJ 9-3 potential (version 2) (cut/shift)",
                       "LJ 9-3 potential (no c/s)",
                       "LJ 9-3 potential with shifted x",
                       "r^9 repulsive potential (no c/s)",
                       "Expotential potential (no c/s)",
                       "Linear external field (no c/s)",
                       "R^7 repulsion plus Yukawa field (cut/shift)")));

     /* INDEPENDENT PARAMEBERS */

        switch(Mesh_List->get("M1_Ndim", (int)NULL)){
             case 1:
                SurfaceInteraction_List->set("SI0_Vext_type","None", "Type of external field to apply with the surfaces", Vext1Dtype_Valid); break;
             case 2:
                SurfaceInteraction_List->set("SI0_Vext_type","None", "Type of external field to apply with the surfaces", Vext2Dtype_Valid); break;
             case 3:
                SurfaceInteraction_List->set("SI0_Vext_type","None", "Type of external field to apply with the surfaces", Vext3Dtype_Valid); break;
             default:
                SurfaceInteraction_List->set("SI0_Vext_type","None", "Type of external field to apply with the surfaces", Vext1Dtype_Valid); break;
        }

        SurfaceInteraction_List->set("SI1_VEXT_MAX",20.0, "Set point for Vext inside surfaces. (VEXT_MAX typically between 10. and 20.)");

        SurfaceInteraction_List->set("SI2_Vext_1D","LJ 9-3 potential (cut/shift)", "Specific analytic external field to apply (1D)", Vext1D_Valid);

        SurfaceInteraction_List->set("SI3_CAREFUL_BOUNDARIES",true,"Use careful element by element approach to integrating through surfaces.  \n Should be true to properly integrate across hard walls, \n but can turn it off if desired.  Should not matter for soft walls.");

        SurfaceInteraction_List->set("SI4_U3D_TO_BE_INTEGRATED","none",
                                     "Specific 3D pair potential to use in numerically integrated potential", U3D_Valid);

        SurfaceInteraction_List->set("SI5_LSemiperm",false,"set to true if any of the surface types are semipermeable membranes.");

        SurfaceInteraction_List->set("SI6_UWW_TF",false,"Calculate wall-wall interactions? Useful for computing g(r).");

        SurfaceInteraction_List->set("SI7_U3D_FOR_WW_INTERACTIONS","LJ 12-6 potential (cut/shift)",
                                     "Specific 3D pair potential to use for wall-wall interactions.\n Useful in computing g(r).", U3D_Valid);

     /* DEPENDENCIES */
        RCP<NumberVisualDependency<int> > VextTypeDep = rcp(
           new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), 
                                            SurfaceInteraction_List->getEntryRCP("SI0_Vext_type")));

        RCP<NumberVisualDependency<int> > VextMaxDep = rcp(
           new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), 
                                            SurfaceInteraction_List->getEntryRCP("SI1_VEXT_MAX")));

        RCP<NumberVisualDependency<int> > LhardDep = rcp(
           new NumberVisualDependency<int>(Surface_List->getEntryRCP("S1: Number of Surfaces"),
                                           SurfaceInteraction_List->getEntryRCP("SI3_CAREFUL_BOUNDARIES")));

        RCP<NumberVisualDependency<int> > Lsemiperm_Dep = rcp(
           new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), 
                                            SurfaceInteraction_List->getEntryRCP("SI5_LSemiperm")));

        RCP<StringVisualDependency> Vext1D_Dep = rcp(
           new StringVisualDependency(SurfaceInteraction_List->getEntryRCP("SI0_Vext_type"),
                                      SurfaceInteraction_List->getEntryRCP("SI2_Vext_1D"), 
                                      tuple<std::string>( "Vext 1D", "Vext 1D_XMIN", "Vext 1D_ORIENTATION"),true));

        RCP<StringVisualDependency> VextU3D_Dep = rcp(
           new StringVisualDependency( SurfaceInteraction_List->getEntryRCP("SI0_Vext_type"),
                                       SurfaceInteraction_List->getEntryRCP("SI4_U3D_TO_BE_INTEGRATED"),
                                       tuple<std::string>( "Numerically integrated Vext from Uij", "Vext for atomic surfaces"),true));

        RCP<StringVisualDependency> UwwTF_Dep = rcp(
           new StringVisualDependency(SurfaceInteraction_List->getEntryRCP("SI0_Vext_type"), 
                                      SurfaceInteraction_List->getEntryRCP("SI6_UWW_TF"),
                                      tuple<std::string>("Vext for atomic surfaces"),true));

        RCP<BoolVisualDependency> VextUWW_Dep = rcp(
           new BoolVisualDependency( SurfaceInteraction_List->getEntryRCP("SI6_UWW_TF"),
                                     SurfaceInteraction_List->getEntryRCP("SI7_U3D_FOR_WW_INTERACTIONS"),true));


     /* DEPENDENCY SHEET ENTRIES*/
        depSheet_Tramonto->addDependency(VextTypeDep);
        depSheet_Tramonto->addDependency(VextMaxDep);
        depSheet_Tramonto->addDependency(LhardDep);
        depSheet_Tramonto->addDependency(Vext1D_Dep);
        depSheet_Tramonto->addDependency(VextU3D_Dep);
        depSheet_Tramonto->addDependency(Lsemiperm_Dep);
        depSheet_Tramonto->addDependency(VextUWW_Dep);
        depSheet_Tramonto->addDependency(UwwTF_Dep);

  /****************************************************************************************************************/
  /****************************** END ANY NEW SUBSECTION PARAMETER SECTION ****************************************/
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

