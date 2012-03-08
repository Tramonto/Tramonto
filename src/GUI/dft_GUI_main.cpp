using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;


extern "C" void dft_OptikaGUI()
{

  bool test_input_file=false;
 /* this file is built based on the examples provided with the optika package with the intent
    of beginning GUI development for Tramonto.  This is a first working version that handles only
    the first few parameters in the Tramonto input file.  In order to even implement the first few
    parameters, dependencies and validators are needed */

           /* Create the empty parameter list called Tramonto_List. *
            * All Tramonto parameters will be added to this list.   */

  RCP<ParameterList> Tramonto_List = RCP<ParameterList>((new ParameterList("Root Tramonto List")));

            /* Create a  dependency sheet for the GUI*/

  RCP<DependencySheet> depSheet_Tramonto = rcp(new DependencySheet());

            /* Create sublists that must be passed around */
  RCP<ParameterList> Mesh_List = sublist(Tramonto_List,"Sect. 1: Computational Domain");
  RCP<ParameterList> Functional_List = sublist(Tramonto_List,"Sect. 2: Functionals");
  RCP<ParameterList> Fluid_List = sublist(Tramonto_List, "Sect. 3: Fluid");
  RCP<ParameterList> PotentialsFF_List = sublist(Fluid_List, "Parameters: Fluid-Fluid Interaction Potentials");
  RCP<ParameterList> Polymer_List = sublist(Fluid_List, "Parameters: Polymers / Bonded Fluids");
  RCP<ParameterList> PolymerCMS_List = sublist(Polymer_List, "P9: Chandler-McCoy-Singer Theory Parameters");
  RCP<ParameterList> PolymerGraft_List = sublist(Polymer_List, "P9: Grafted Polymer Parameters Parameters");
  RCP<ParameterList> PolymerArch_List = sublist(Polymer_List, "P9: Architecture Entry for Bonded Systems");
  RCP<ParameterList> StatePoint_List = sublist(Fluid_List, "Parameters: Properties of the Bulk Fluid");
  RCP<ParameterList> Diffusion_List = sublist(Fluid_List, "Special (optional) Parameters for Diffusing Systems");
  RCP<ParameterList> ChargedFluid_List = sublist(StatePoint_List, "Parameters: Fluids with Charge / Electrostatics");
  RCP<ParameterList> Surface_List = sublist(Tramonto_List, "Sect. 4: Surfaces");
  RCP<ParameterList> SurfaceGeometry_List = sublist(Surface_List, "4.1: Geometry Parameters");
  RCP<ParameterList> PotentialsWW_List = sublist(Surface_List, "4.2: Interaction Parameters");
  RCP<ParameterList> SurfaceParamCharge_List = sublist(Surface_List, "4.3 Charged Surface Parameters");
  RCP<ParameterList> SurfaceInteraction_List = sublist(Surface_List, "4.4 Fluid-Surface Interactions");
  RCP<ParameterList> PotentialsWF_List = sublist(Surface_List, "4.5 Fluid-Surface Potential Parameters");
  RCP<ParameterList> Continuation_List = sublist(Tramonto_List, "Sect. 6: Continuation");
  RCP<ParameterList> Startup_List = sublist(Tramonto_List, "Sect. 7: Startup Control");
  RCP<ParameterList> Output_List = sublist(Tramonto_List, "Sect. 8: Output Control");
  RCP<ParameterList> Solver_List = sublist(Tramonto_List, "Sect. 9: Numerical Methods");
  RCP<ParameterList> Coarsening_List = sublist(Solver_List, "Numerical Coarsening Options");
  RCP<ParameterList> LoadBalance_List = sublist(Solver_List, "Load Balancing Options");
  RCP<ParameterList> PhysicsMethod_List = sublist(Solver_List, "Physics based Solver Options");
  RCP<ParameterList> LinearSolver_List = sublist(Solver_List, "Linear Solver Options");
  RCP<ParameterList> NonlinearSolver_List = sublist(Solver_List, "Nonlinear Solver Options");

  if (test_input_file){

/*      Teuchos::RCP<Teuchos::ParameterList> myParameters;*/

  /**
   * Then we just call getInput! There's a little more to it, so let's
   * head on over to the inputs.xml file to see what's going on there.
   */
        Optika::getInput("inputs.xml", Tramonto_List);
  }
  else{

            /* Call different sections of the GUI.   *
             * Note that each file contains the bits *
             * that were formerly found in different *
             * sections of dft_input.dat             */

  dft_GUI_mesh(Tramonto_List,depSheet_Tramonto,Mesh_List);

  dft_GUI_functionals(Tramonto_List,depSheet_Tramonto,Functional_List);

  dft_GUI_potentialsFF(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,PotentialsFF_List);

  dft_GUI_Polymer(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List,PolymerCMS_List,PolymerArch_List,PolymerGraft_List);

  dft_GUI_StatePoint(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List,
                                 StatePoint_List,Diffusion_List,ChargedFluid_List);

  dft_GUI_surfaces(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfaceGeometry_List,PotentialsWW_List,SurfaceParamCharge_List);


  dft_GUI_vextType(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfaceInteraction_List);

  dft_GUI_potentialsWF(Tramonto_List,depSheet_Tramonto,Functional_List, Surface_List,
                       SurfaceInteraction_List,Fluid_List,PotentialsWF_List);

  dft_GUI_Continuation(Tramonto_List,depSheet_Tramonto,Functional_List, PotentialsFF_List, StatePoint_List,
                       Continuation_List);

  dft_GUI_NumericalMethods(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List,Solver_List,
                            Coarsening_List,LoadBalance_List,PhysicsMethod_List,NonlinearSolver_List,LinearSolver_List);

  dft_GUI_OutputParams(Tramonto_List,depSheet_Tramonto,Mesh_List,Functional_List,Fluid_List,Continuation_List,SurfaceInteraction_List,Output_List);

             /* The getInput function starts up an Optika GUI and      *
              * lets the user start to input parameter values. When    *
              * the user has completed their data entry, the function  *
              * will finish right after all of the input values are    *
              * stored in My_List.                                     */

  }

  Optika::getInput(Tramonto_List,depSheet_Tramonto);

             /* Print out what the user entered in nice XML format.  */

/*  RCP<FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  writeParameterListToXmlOStream(*Tramonto_List, *out);*/

/*  cout<<"calling dft_GUI_toTramonto\n"<<endl;*/
  dft_GUI_toTramonto(Tramonto_List,Mesh_List,Functional_List,Fluid_List,
                     PotentialsFF_List,Polymer_List,PolymerGraft_List,PolymerArch_List,PolymerCMS_List,
                     StatePoint_List,Diffusion_List,ChargedFluid_List,Continuation_List,
                     Solver_List,Coarsening_List,LoadBalance_List,PhysicsMethod_List,LinearSolver_List,NonlinearSolver_List,
                     Surface_List,SurfaceGeometry_List);

             /* Here save parameter a to return to C code --- a fully 
                functioning GUI will need to return all parameters entered
                by the user to Tramonto properly. */

  return;
}

