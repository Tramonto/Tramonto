using namespace std;
#include <iostream>
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

/*  RCP<DependencySheet> depSheet_Tramonto = rcp(new DependencySheet(Tramonto_List));*/
  RCP<DependencySheet> depSheet_Tramonto = rcp(new DependencySheet());

            /* Create sublists that must be passed around */




            /* Call different sections of the GUI.   *
             * Note that each file contains the bits *
             * that were formerly found in different *
             * sections of dft_input.dat             */

  RCP<ParameterList> Mesh_List = sublist(Tramonto_List,"Sect. 1: Computational Domain");
  dft_GUI_mesh(Tramonto_List,depSheet_Tramonto,Mesh_List);

  RCP<ParameterList> Functional_List = sublist(Tramonto_List,"Sect. 2: Functionals");
  dft_GUI_functionals(Tramonto_List,depSheet_Tramonto,Functional_List);

  RCP<ParameterList> Fluid_List = sublist(Tramonto_List, "Sect. 3: Fluid");
  RCP<ParameterList> PotentialsFF_List = sublist(Fluid_List, "F7_Fluid-Fluid Potential Parameters");
  dft_GUI_potentialsFF(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,PotentialsFF_List);

  RCP<ParameterList> Polymer_List = sublist(Fluid_List, "F8_Polymer(Bonded) Fluid Parameters");
  RCP<ParameterList> PolymerCMS_List = sublist(Polymer_List, "P9: Chandler-McCoy-Singer Theory Parameters");
  RCP<ParameterList> PolymerGraft_List = sublist(Polymer_List, "P9: Grafted Polymer Parameters Parameters");
  RCP<ParameterList> PolymerArch_List = sublist(Polymer_List, "P9: Architecture Entry for Bonded Systems");
  dft_GUI_Polymer(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List,PolymerCMS_List,PolymerArch_List,PolymerGraft_List);

  RCP<ParameterList> StatePoint_List = sublist(Fluid_List, "F9_Bulk Fluid Properties");
  RCP<ParameterList> Diffusion_List = sublist(Fluid_List, "Parameters for Diffusing Systems");
  RCP<ParameterList> ChargedFluid_List = sublist(StatePoint_List, "BF7: Parameters for Charged Systems");
  dft_GUI_StatePoint(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List,
                                 StatePoint_List,Diffusion_List,ChargedFluid_List);


  RCP<ParameterList> Surface_List = sublist(Tramonto_List, "Sect. 4: Surfaces");
  RCP<ParameterList> SurfaceGeometry_List = sublist(Surface_List, "4.1: Geometry Parameters");
  RCP<ParameterList> PotentialsWW_List = sublist(Surface_List, "4.2: Interaction Parameters");
  RCP<ParameterList> SurfaceParamCharge_List = sublist(Surface_List, "4.3 Charged Surface Parameters");
  dft_GUI_surfaces(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfaceGeometry_List,PotentialsWW_List,SurfaceParamCharge_List);


  RCP<ParameterList> SurfaceInteraction_List = sublist(Surface_List, "4.4 Fluid-Surface Interactions");
  dft_GUI_vextType(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfaceInteraction_List);

  RCP<ParameterList> PotentialsWF_List = sublist(Surface_List, "4.5 Fluid-Surface Potential Parameters");
  dft_GUI_potentialsWF(Tramonto_List,depSheet_Tramonto,Functional_List, Surface_List,
                       SurfaceInteraction_List,Fluid_List,PotentialsWF_List);


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

/*  RCP<FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();*/
/*  writeParameterListToXmlOStream(*Tramonto_List, *out);*/

/*  cout<<"calling dft_GUI_toTramonto\n"<<endl;*/
  dft_GUI_toTramonto(Tramonto_List,Mesh_List,Functional_List,Fluid_List,
                     PotentialsFF_List,Polymer_List,PolymerGraft_List,PolymerArch_List,PolymerCMS_List,
                     StatePoint_List,Diffusion_List,ChargedFluid_List,
                     Surface_List,SurfaceGeometry_List);

             /* Here save parameter a to return to C code --- a fully 
                functioning GUI will need to return all parameters entered
                by the user to Tramonto properly. */

  return;
}

