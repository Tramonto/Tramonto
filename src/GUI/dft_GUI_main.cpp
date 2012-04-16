using namespace std;
#include <iostream>
#include "dft_globals_const.h" 
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

/*int test_NsurfTypes(int,int);*/

extern "C" void dft_OptikaGUI()
{
  int i;
  bool test_input_file=true;
  bool read_input_old_format=true;

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
  RCP<ParameterList> PotentialsFF_List = sublist(Fluid_List, "FL1: Fluid-Fluid Interactions");
  RCP<ParameterList> Polymer_List = sublist(Fluid_List, "FL2: Polymers or Bonded Fluids");
  RCP<ParameterList> StatePoint_List = sublist(Fluid_List, "FL3: Properties of the Bulk Fluid");
  RCP<ParameterList> Diffusion_List = sublist(Fluid_List, "FL4: Diffusion Params (optional)");


  RCP<ParameterList> PolymerCMS_List = sublist(Polymer_List, "PL1: Chandler-McCoy-Singer Theory Params");
  RCP<ParameterList> PolymerGraft_List = sublist(Polymer_List, "PL2: Grafted Polymer Parameters Params");
  RCP<ParameterList> PolymerArch_List = sublist(Polymer_List, "PL3: Architecture Entry Bonded Systems");

  RCP<ParameterList> ChargedFluid_List = sublist(StatePoint_List, "BFL1: Fluids with Charge and Electrostatics");

  RCP<ParameterList> Surface_List = sublist(Tramonto_List, "Sect. 4: Surfaces");

  RCP<ParameterList> SurfGeom0_List = sublist(Surface_List, "SL1: GeomParams SurfType 0");
  RCP<ParameterList> SurfGeom1_List = sublist(Surface_List, "SL1: GeomParams SurfType 1");
  RCP<ParameterList> SurfGeom2_List = sublist(Surface_List, "SL1: GeomParams SurfType 2");
  RCP<ParameterList> SurfGeom3_List = sublist(Surface_List, "SL1: GeomParams SurfType 3");
  RCP<ParameterList> SurfGeom4_List = sublist(Surface_List, "SL1: GeomParams SurfType 4");
  RCP<ParameterList> SurfGeom5_List = sublist(Surface_List, "SL1: GeomParams SurfType 5");
  RCP<ParameterList> SurfGeom6_List = sublist(Surface_List, "SL1: GeomParams SurfType 6");
  RCP<ParameterList> SurfGeom7_List = sublist(Surface_List, "SL1: GeomParams SurfType 7");
  RCP<ParameterList> SurfGeom8_List = sublist(Surface_List, "SL1: GeomParams SurfType 8");
  RCP<ParameterList> SurfGeom9_List = sublist(Surface_List, "SL1: GeomParams SurfType 9");
  RCP<ParameterList> SurfGeom10_List = sublist(Surface_List, "SL1: GeomParams SurfType 9");
  RCP<ParameterList> SurfGeom11_List = sublist(Surface_List, "SL1: GeomParams SurfType 9");

  RCP<ParameterList> PotentialsWW_List = sublist(Surface_List, "SL2: Interaction Params");
  RCP<ParameterList> SurfaceParamCharge_List = sublist(Surface_List, "SL3: Charged Surface Params");
  RCP<ParameterList> SurfaceInteraction_List = sublist(Surface_List, "SL4: Fluid-Surface Interactions");
  RCP<ParameterList> PotentialsWF_List = sublist(Surface_List, "SL5: Fluid-Surface Potential Params");

  RCP<ParameterList> Continuation_List = sublist(Tramonto_List, "Sect. 6: Continuation");
  RCP<ParameterList> DensProfile_List = sublist(Tramonto_List, "Sect. 7: Initial Guess Options");
  RCP<ParameterList> Output_List = sublist(Tramonto_List, "Sect. 8: Output Control");
  RCP<ParameterList> Solver_List = sublist(Tramonto_List, "Sect. 9: Numerical Methods");

  RCP<ParameterList> NonlinearSolver_List = sublist(Solver_List, "SL1: Nonlinear Solver Options");
  RCP<ParameterList> LinearSolver_List = sublist(Solver_List, "SL2: Linear Solver Options");
  RCP<ParameterList> PhysicsMethod_List = sublist(Solver_List, "SL3: Physics based Solver Options");
  RCP<ParameterList> Coarsening_List = sublist(Solver_List, "SL4: Numerical Coarsening Options");
  RCP<ParameterList> LoadBalance_List = sublist(Solver_List, "SL5: Load Balancing Options");

  if (test_input_file){

  /**
   * Then we just call getInput! There's a little more to it, so let's
   * head on over to the inputs.xml file to see what's going on there.
   */
cout << "reading inputs.xml" << "test_input_file="<<test_input_file<<endl;
        Optika::getInput("inputs.xml", Tramonto_List);
  }
  else{

    /* set defaults for each section of the GUI */
     dft_GUI_mesh_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List);
     dft_GUI_functionals_set_defaults(Tramonto_List,depSheet_Tramonto,Functional_List);
     dft_GUI_potentialsFF_set_defaults(Tramonto_List,depSheet_Tramonto,Fluid_List,PotentialsFF_List);
     dft_GUI_Polymer_set_defaults(Tramonto_List,depSheet_Tramonto,Fluid_List,Polymer_List,PolymerCMS_List,PolymerArch_List,PolymerGraft_List);
     dft_GUI_StatePoint_set_defaults(Tramonto_List,depSheet_Tramonto,Fluid_List,Polymer_List, StatePoint_List,Diffusion_List,ChargedFluid_List);
     dft_GUI_surfaces_set_defaults(Tramonto_List,depSheet_Tramonto,Surface_List);
     for (i=0;i<10;i++){
             if (i==0) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom0_List);
        else if (i==1) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom1_List);
        else if (i==2) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom2_List);
        else if (i==3) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom3_List);
        else if (i==4) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom4_List);
        else if (i==5) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom5_List);
        else if (i==6) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom6_List);
        else if (i==7) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom7_List);
        else if (i==8) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom8_List);
        else if (i==9) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom9_List);
        else if (i==9) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom10_List);
        else if (i==9) dft_GUI_surface_geometry_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom11_List);
     }
     dft_GUI_Continuation_set_defaults(Tramonto_List,depSheet_Tramonto,Continuation_List);
     dft_GUI_NumericalMethods_set_defaults(Tramonto_List,depSheet_Tramonto,Fluid_List,Solver_List,
                              Coarsening_List,LoadBalance_List,PhysicsMethod_List,NonlinearSolver_List,LinearSolver_List);

    /* replace defauls if requested for each section of the GUI */
     if (read_input_old_format){ 
        dft_GUI_mesh_set_OldFormat(Tramonto_List,depSheet_Tramonto,Mesh_List);
        dft_GUI_functionals_set_OldFormat(Tramonto_List,depSheet_Tramonto,Functional_List);
        dft_GUI_potentialsFF_set_OldFormat(Tramonto_List,depSheet_Tramonto,Fluid_List,PotentialsFF_List);
        dft_GUI_Polymer_set_OldFormat(Tramonto_List,depSheet_Tramonto,Fluid_List,Polymer_List,PolymerCMS_List,PolymerArch_List,PolymerGraft_List);
        dft_GUI_StatePoint_set_OldFormat(Tramonto_List,depSheet_Tramonto,Fluid_List,Polymer_List, StatePoint_List,Diffusion_List,ChargedFluid_List);
        dft_GUI_surfaces_set_OldFormat(Tramonto_List,depSheet_Tramonto,Surface_List);
        for (i=0;i<Nwall_type;i++){
                if (i==0) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom0_List);
           else if (i==1) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom1_List);
           else if (i==2) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom2_List);
           else if (i==3) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom3_List);
           else if (i==4) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom4_List);
           else if (i==5) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom5_List);
           else if (i==6) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom6_List);
           else if (i==7) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom7_List);
           else if (i==8) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom8_List);
           else if (i==9) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom9_List);
           else if (i==9) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom10_List);
           else if (i==9) dft_GUI_surface_geometry_set_OldFormat(i,Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom11_List);
        }
        dft_GUI_Continuation_set_OldFormat(Tramonto_List,depSheet_Tramonto,Continuation_List);
        dft_GUI_NumericalMethods_set_OldFormat(Tramonto_List,depSheet_Tramonto,Solver_List,
                                 Coarsening_List,LoadBalance_List,PhysicsMethod_List,NonlinearSolver_List,LinearSolver_List);
     }

    /* set dependencies for each section of the GUI */
     dft_GUI_mesh_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List);
     dft_GUI_functionals_dependencies(Tramonto_List,depSheet_Tramonto,Functional_List);
     dft_GUI_potentialsFF_dependencies(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,PotentialsFF_List);
     dft_GUI_Polymer_dependencies(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List,PolymerCMS_List,PolymerArch_List,PolymerGraft_List);
     dft_GUI_StatePoint_dependencies(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List, StatePoint_List,Diffusion_List,ChargedFluid_List);
     dft_GUI_surfaces_dependencies(Tramonto_List,depSheet_Tramonto,Surface_List);
     for (i=0;i<10;i++){
             if (i==0) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom0_List);
        else if (i==1) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom1_List);
        else if (i==2) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom2_List);
        else if (i==3) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom3_List);
        else if (i==4) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom4_List);
        else if (i==5) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom5_List);
        else if (i==6) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom6_List);
        else if (i==7) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom7_List);
        else if (i==8) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom8_List);
        else if (i==9) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom9_List);
        else if (i==9) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom10_List);
        else if (i==9) dft_GUI_surface_geometry_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfGeom11_List);
     }
     dft_GUI_Continuation_dependencies(Tramonto_List,depSheet_Tramonto,PotentialsFF_List,Continuation_List);
     dft_GUI_NumericalMethods_dependencies(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List,Solver_List,
                              Coarsening_List,LoadBalance_List,PhysicsMethod_List,NonlinearSolver_List,LinearSolver_List);

  dft_GUI_potentialsWW(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,PotentialsWW_List,SurfaceParamCharge_List);
  dft_GUI_vextType(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,SurfaceInteraction_List);
  dft_GUI_potentialsWF(Tramonto_List,depSheet_Tramonto,Functional_List, Surface_List,
                       SurfaceInteraction_List,Fluid_List,PotentialsWF_List); 


  dft_GUI_DensityStartupParams(Tramonto_List,depSheet_Tramonto,Fluid_List,Continuation_List,DensProfile_List);
  dft_GUI_OutputParams(Tramonto_List,depSheet_Tramonto,Mesh_List,Functional_List,Fluid_List,Continuation_List,SurfaceInteraction_List,Output_List);

             /* The getInput function starts up an Optika GUI and      *
              * lets the user start to input parameter values. When    *
              * the user has completed their data entry, the function  *
              * will finish right after all of the input values are    *
              * stored in My_List.                                     */

  }

/* set up list dependencies */
  
   RCP<StringVisualDependency> PolyList_Dep = rcp(new StringVisualDependency(
          Functional_List->getEntryRCP("F4_POLYMER_Functional"), Fluid_List->getEntryRCP("FL2: Polymers or Bonded Fluids"), 
          tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT",
                             "Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

   RCP<StringVisualDependency> DiffusionList_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F0_Type_of_Calculation"),Fluid_List->getEntryRCP("FL4: Diffusion Params (optional)"),
           tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)","Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringVisualDependency> ChargedFluidList_Dep = rcp( new StringVisualDependency(
           Functional_List->getEntryRCP("F3_CHARGE_Functional"),
           StatePoint_List->getEntryRCP("BFL1: Fluids with Charge and Electrostatics"),
           tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM", "Charge_with_DeltaC_General","Charge(MF)_with_Polarization")));

/*   RCP<SimpleFunctionObject<int> >test_NsurfTypes(int,int);
   RCP<NumberVisualDependency<int> > SurfaceTypeList_DepN8=rcp(new NumberVisualDependency<int>(
         Surface_List->getEntry("S3: Number of surface types"),
         Surface_List->getEntry("SL1:GeomParams SurfType 7"),true,
         test_NsurfTypes(Surface_List->get<int>("S3: Number of surface types"),8) ));*/
  
   depSheet_Tramonto->addDependency(PolyList_Dep);
   depSheet_Tramonto->addDependency(DiffusionList_Dep);
   depSheet_Tramonto->addDependency(ChargedFluidList_Dep);

             /* Greate the GUI windows for data entry with all dependencies.  */
    Optika::getInput(Tramonto_List,depSheet_Tramonto);

             /* Print out what the user entered in nice XML format.  */

  RCP<FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  writeParameterListToXmlOStream(*Tramonto_List, *out);

  dft_GUI_toTramonto(Tramonto_List,Mesh_List,Functional_List,Fluid_List,
                     PotentialsFF_List,Polymer_List,PolymerGraft_List,PolymerArch_List,PolymerCMS_List,
                     StatePoint_List,Diffusion_List,ChargedFluid_List,Continuation_List,
                     Solver_List,Coarsening_List,LoadBalance_List,PhysicsMethod_List,LinearSolver_List,NonlinearSolver_List,Output_List,
                     DensProfile_List,Surface_List,SurfGeom0_List, SurfGeom1_List, SurfGeom2_List, SurfGeom3_List, SurfGeom4_List, SurfGeom5_List,
                     SurfGeom6_List, SurfGeom7_List, SurfGeom8_List, SurfGeom9_List, SurfGeom10_List,SurfGeom11_List);

             /* Here save parameter a to return to C code --- a fully 
                functioning GUI will need to return all parameters entered
                by the user to Tramonto properly. */

  return;
}
/*********************************************************************************************************************/
/*int test_NsurfTypes(int Ninput,int Ntest){

   if (Ninput >= Ntest) return Ntest;
   else return -1;
}*/
/*********************************************************************************************************************/

