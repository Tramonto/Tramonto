/*#include "Teuchos_ParameterList.hpp"*/
using namespace std;
#include <iostream>
/*#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"*/
#include "Teuchos_Version.hpp"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;


extern "C" void dft_OptikaGUI()
{
 /* this file is built based on the examples provided with the optika package with the intent
    of beginning GUI development for Tramonto.  This is a first working version that handles only
    the first few parameters in the Tramonto input file.  In order to even implement the first few
    parameters, dependencies and validators are needed */

           /* Create the empty parameter list called Tramonto_List. *
            * All Tramonto parameters will be added to this list.   */

  RCP<ParameterList> Tramonto_List = RCP<ParameterList>((new ParameterList("Root Tramonto List")));

            /* Create a  dependency sheet for the GUI*/

  RCP<DependencySheet> depSheet_Tramonto = rcp(new DependencySheet(Tramonto_List));

            /* Create sublists that must be passed around */




            /* Call different sections of the GUI.   *
             * Note that each file contains the bits *
             * that were formerly found in different *
             * sections of dft_input.dat             */

  RCP<ParameterList> Mesh_List = sublist(Tramonto_List,"Sect. 1: Computational Domain");
  dft_GUI_mesh(Tramonto_List,depSheet_Tramonto,Mesh_List);

  RCP<ParameterList> Functional_List = sublist(Tramonto_List,"Sect. 2: Functional Types");
  dft_GUI_functionals(Tramonto_List,depSheet_Tramonto,Functional_List);

  RCP<ParameterList> Fluid_List = sublist(Tramonto_List, "Sect. 3: Fluid");
  RCP<ParameterList> PotentialsFF_List = sublist(Fluid_List, "Fluid-Fluid Potential Parameters");
  dft_GUI_potentialsFF(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,PotentialsFF_List);
  RCP<ParameterList> Polymer_List = sublist(Fluid_List, "Bonded Fluid Parameters");
  RCP<ParameterList> ChargedFluid_List = sublist(Fluid_List, "Charged Fluid Parameters");
  RCP<ParameterList> Diffusion_List = sublist(Fluid_List, "Charged Fluid Parameters");

  RCP<ParameterList> Surface_List = sublist(Tramonto_List, "Sect. 4: Surfaces");
  RCP<ParameterList> SurfaceGeometry_List = sublist(Surface_List, "Surface Geometry Parameters");
  dft_GUI_surfaces(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfaceGeometry_List);

  RCP<ParameterList> SurfaceInteraction_List = sublist(Surface_List, "Surface-Fluid Interaction Types");
  dft_GUI_vextType(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfaceInteraction_List);

  RCP<ParameterList> PotentialsWF_List = sublist(Surface_List, "Surface-Fluid Potential Parameters");
  dft_GUI_potentialsWF(Tramonto_List,depSheet_Tramonto,Functional_List, Surface_List,
                       SurfaceInteraction_List,Fluid_List,PotentialsWF_List);

  RCP<ParameterList> ChargedSurface_List = sublist(Surface_List, "Charged Surface Parameters");
  RCP<ParameterList> PotentialsWW_List = sublist(Surface_List, "Surface-Fluid Potential Parameters");

  RCP<ParameterList> Thermodynamics_List = sublist(Tramonto_List, "Sect. 5: State Point(s)");
  RCP<ParameterList> Continuation_List = sublist(Tramonto_List, "Sect. 6: Continuation");
  RCP<ParameterList> Startup_List = sublist(Tramonto_List, "Sect. 7: Startup Control");

  RCP<ParameterList> Output_List = sublist(Tramonto_List, "Sect. 8: Output Control");
  dft_GUI_OutputParams(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,SurfaceInteraction_List,Output_List);

  RCP<ParameterList> Solver_List = sublist(Tramonto_List, "Sect. 9: Numerical Methods");
  RCP<ParameterList> Coarsening_List = sublist(Solver_List, "Jacobian & Mesh Coarsening");
  RCP<ParameterList> LinearSolver_List = sublist(Solver_List, "Linear Solver Options");
  RCP<ParameterList> NonlinearSolver_List = sublist(Solver_List, "Nonlinear Solver Options");
  dft_GUI_NumericalMethods(Tramonto_List,depSheet_Tramonto,Functional_List,Solver_List,
                            Coarsening_List,NonlinearSolver_List,LinearSolver_List);

             /* The getInput function starts up an Optika GUI and      *
              * lets the user start to input parameter values. When    *
              * the user has completed their data entry, the function  *
              * will finish right after all of the input values are    *
              * stored in My_List.                                     */

  Optika::getInput(Tramonto_List,depSheet_Tramonto);

             /* Print out what the user entered in nice XML format.  */

  RCP<FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  writeParameterListToXmlOStream(*Tramonto_List, *out);

  dft_GUI_toTramonto(Tramonto_List,Mesh_List,Functional_List,Surface_List,SurfaceGeometry_List);

             /* Here save parameter a to return to C code --- a fully 
                functioning GUI will need to return all parameters entered
                by the user to Tramonto properly. */

  return;
}

