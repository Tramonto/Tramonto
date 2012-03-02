using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_OutputParams(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Output_List)
{
    /****************************************************************************************************************/
  /****************************** FUNCTIONAL CONTROL PARAMETER SECTION ********************************************/
  /****************************************************************************************************************/
    /**************************************/
    /* Define validators for this section.*/
    /**************************************/
    RCP<StringValidator> OutputTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Minimal","Densities", 
						  "Extended","Verbose","No Screen")));

    RCP<StringValidator> StateOutputValidator = rcp( new StringValidator(tuple<std::string>("Density[icont]","Betamu[icont]","kappa","Density[icont], Betamu[icont], and kappa","Density[Ncomp], Betamu[Ncomp], kappa")));

    RCP<StringValidator> StateOutputNoChargeValidator = rcp( new StringValidator(tuple<std::string>("Density[icont]","Betamu[icont]","Density[icont] and Betamu[icont]","Density[Ncomp] and Betamu[Ncomp]")));

    RCP<StringValidator> AdsOutputValidator = rcp( new StringValidator(tuple<std::string>("total adsorption","excess adsorption", "excess and total adsorption","adsorption/volume (density in bulk fluid)")));
    RCP<StringValidator> EnergyOutputValidator = rcp( new StringValidator(tuple<std::string>("total free energy","excess surface free energy", "excess and total free energy","free energy/volume (pressure in bulk fluid)")));

    RCP<StringValidator> MeshOutputValidator = rcp( new StringValidator(tuple<std::string>("Positions","Separations")));



    /***************************************************************/
    /* set up independent parameters that do not have dependencies */
    /***************************************************************/

    Output_List->set("O1: Output Type", "Extended", "Select how much output you would like to generate.\n Extended output will give external fields and segment densities (for polymers).\n Verbose output will give all ancillary fields in the calculation.", OutputTypeValidator);

    Output_List->set("O2: Energies per unit area?", false, "Set to true for output (energy, adsorption, force) per unit area.\n If false, output will be energy (3D), energy/Length (2D), and energy/Area (1D).");
    Output_List->set("O3: Count reflections?", true, "Set to true if you want output (energy/adsorption/force) to use reflected images in calculations.");
    Output_List->set("O4: Print radial correlation function: g(r)?", false, "Set to true to produce a radial correlation function.");
    Output_List->set("O5: Print surface-surface interactions?", false, "Set to true to produce output for surface-surface interactions.");

    if (Functional_List->get<string>("F3_CHARGE_Functional") !="No Charge or No Poisson" || 
        Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut/shift)" ||
        Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut only)" ){ 
        Output_List->set("O6: Type of State point Output", "Density[icont]", "Select output type for state point in file dft_output.dat.\n Options are densities, chemical potentials, or kappa(for ionic systems)",StateOutputValidator);
    }
    else{
        Output_List->set("O6: Type of State point Output", "Density[icont]", "Select output type for state point in file dft_output.dat.\n Options are densities or chemical potentials",StateOutputNoChargeValidator);
    }
    Output_List->set("O7: Type of Mesh Output", "Positions", "Select output type for mesh continuation.\n Options are to output positions of surfaces or separations of surfaces.",MeshOutputValidator);

    /*******************************/
    /* define dependent parameters */
    /*******************************/

    /**********************************************************************************************/
    /* show the dependent parameters only if the independent parameters has a particular setting. */
    /**********************************************************************************************/
      RCP<StringVisualDependency> GofRprint_Dep = rcp(
           new StringVisualDependency( SurfaceInteraction_List->getEntryRCP("SI0_Vext_type"),Output_List->getEntryRCP("O4: Print radial correlation function: g(r)?"),  
               tuple<std::string>("Vext for atomic surfaces")));

      RCP<StringVisualDependency> UWWprint_Dep = rcp(
           new StringVisualDependency(SurfaceInteraction_List->getEntryRCP("SI0_Vext_type"),Output_List->getEntryRCP("O5: Print surface-surface interactions?"),
               tuple<std::string>("Vext for atomic surfaces")));

    /*****************************************/
    /* add the dependencies for this section.*/
    /*****************************************/
      depSheet_Tramonto->addDependency(GofRprint_Dep);
      depSheet_Tramonto->addDependency(UWWprint_Dep);

  /****************************************************************************************************************/
  /****************************** END FUNCTIONAL CONTROL PARAMETER SECTION ****************************************/
  /****************************************************************************************************************/
 
  return;
}

