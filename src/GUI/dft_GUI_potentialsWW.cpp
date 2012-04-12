using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_potentialsWW(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List)
{

  /****************************************************************************************************************/
  /********************************* SURFACE DEFINITIONS SECTION **************************************************/
  /****************************************************************************************************************/

     /* VALIDATORS*/

   RCP<StringValidator> SurfElecBCType_Validator = rcp(new StringValidator(tuple<std::string>(
	   "none - neutral surfaces",
           "constant potential surfaces",
           "constant surface charge density",
           "discrete atomic charges")));

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


     /* INDEPENDENT PARAMEBERS */
       PotentialsWW_List->set("WW1: Compute wall-wall interactions?",false,"Indicate whether to compute wall-wall interactions (only for 3D systems with atomic surfaces)" );
       PotentialsWW_List->set("WW2: Type of wall-wall interactions","none","Identify type of wall-wall interactions.",PairPotValidator);

       Array<int> UWWArray(Surface_List->get<int>("S3: Number of surface types"),1.0);
       PotentialsWW_List->set("WW3: Sigma_ww[iwall_type][jwall_type]",UWWArray,"Array for characteristic length in wall-wall interactions.");
       PotentialsWW_List->set("WW4: Eps_ww[iwall_type][jwall_type]",UWWArray,"Array for energy parameters in wall-wall interactions.");
       PotentialsWW_List->set("WW5: Cut_ww[iwall_type][jwall_type]",UWWArray,"Array for cutoff distances in wall-wall interactions.");
       PotentialsWW_List->set("WW6: EpsYukawa_ww[iwall_type][jwall_type]",UWWArray,"Array for prefactor of Yukawa term in wall-wall interactions.");
       PotentialsWW_List->set("WW7: YukawaK_ww[iwall_type][jwall_type]",UWWArray,"Array for Yukawa Decay Parameter in wall-wall interactions.");

       SurfaceParamCharge_List->set("SC1: Type_elec_BC","none - neutral surfaces","Select type of boundary condition you would like to apply at the charged surfaces.",SurfElecBCType_Validator);
       SurfaceParamCharge_List->set("SC2: N_Localized_Charges",0,"Enter number of localized fixed source charges in the system.");


     /* DEPENDENCIES */
        /* need to write a dependency for "WW1: Compute wall-wall interactions?" to show up if Ndim=3D and at least one surface type is atomic */
        RCP<BoolVisualDependency> UWWType_Dep = rcp(new BoolVisualDependency(
               PotentialsWW_List->getEntryRCP("WW1: Compute wall-wall interactions?"), 
	       PotentialsWW_List->getEntryRCP("WW2: Type of wall-wall interactions"), true));

     /* DEPENDENCY SHEET ENTRIES*/
  return;
}

