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
  int i,n;
  bool read_xml_file=false, read_input_old_format=false,start_GUI=true;
  string runpath, InputOLD_File,InputXML_File, InputGUI_File;


  /* set up a parameter list for a very simple GUI to select the type of run we are interested in */
  RCP<ParameterList> RunType_List = RCP<ParameterList>((new ParameterList("Root RunType List")));
  RCP<DependencySheet> depSheet_RunType = rcp(new DependencySheet());

  RCP<StringValidator> RunTypeValidator = rcp(
        new StringValidator(tuple<std::string>(
           "GUI - Defaults", 
           "GUI - Old Format File Parameters",
           "GUI - XML File Parameters",
           "No GUI / Old Format File",
           "No GUI / XML File")));
 
  RCP<FileNameValidator> InputXMLFileVal = rcp (new FileNameValidator);
  RCP<FileNameValidator> InputOLDFileVal = rcp (new FileNameValidator);
  RCP<FileNameValidator> InputGUIFileVal = rcp (new FileNameValidator);

  RunType_List->set("R1: Run Type","GUI - Defaults","Indicate how you would like to run Tramonto for this job.",RunTypeValidator);

  RunType_List->set("R2: Input File (OLD format)","","select the OLD format input file. This choice also sets run directory.",InputOLDFileVal);
  RunType_List->set("R2: Input File (XML)","","select the XML input file. This choice also sets run directory.",InputXMLFileVal);
  RunType_List->set("R2: Select Any File in Desired Directory","","select any file in desired directory.  The runpath will be extracted from the filename.",InputGUIFileVal);

   RCP<StringVisualDependency> RunType_Dep1 = rcp(new StringVisualDependency(
          RunType_List->getEntryRCP("R1: Run Type"), RunType_List->getEntryRCP("R2: Select Any File in Desired Directory"), 
          tuple<std::string>("GUI - Defaults")));

   RCP<StringVisualDependency> RunType_Dep2 = rcp(new StringVisualDependency(
          RunType_List->getEntryRCP("R1: Run Type"), RunType_List->getEntryRCP("R2: Input File (XML)"), 
          tuple<std::string>("GUI - XML File Parameters","No GUI / XML File")));

   RCP<StringVisualDependency> RunType_Dep3 = rcp(new StringVisualDependency(
          RunType_List->getEntryRCP("R1: Run Type"), RunType_List->getEntryRCP("R2: Input File (OLD format)"), 
          tuple<std::string>("GUI - Old Format File Parameters","No GUI / Old Format File")));

   depSheet_RunType->addDependency(RunType_Dep1);
   depSheet_RunType->addDependency(RunType_Dep2);
   depSheet_RunType->addDependency(RunType_Dep3);
  
 /* Greate the GUI window for run startup for Tramonto  */
  Optika::getInput(RunType_List,depSheet_RunType);


  if (RunType_List->get<string>("R1: Run Type")=="GUI - Defaults"){
     InputGUI_File=RunType_List->get<string>("R2: Select Any File in Desired Directory");
     string::size_type n=InputGUI_File.find_last_of("/");
     runpath=InputGUI_File.substr(0,n);
  }
  else if (RunType_List->get<string>("R1: Run Type")=="GUI - Old Format File Parameters"){
     InputOLD_File=RunType_List->get<string>("R2: Input File (OLD format)");
     string::size_type n=InputOLD_File.find_last_of("/");
     runpath=InputOLD_File.substr(0,n);
     read_input_old_format=true;
  }
  else if (RunType_List->get<string>("R1: Run Type")=="GUI - XML File Parameters"){
     InputXML_File=RunType_List->get<string>("R2: Input File (XML)");
     string::size_type n=InputXML_File.find_last_of("/");
     runpath=InputXML_File.substr(0,n);
     read_xml_file=true;
  }
  else if (RunType_List->get<string>("R1: Run Type")=="No GUI / Old Format File"){
     InputOLD_File=RunType_List->get<string>("R2: Input File (OLD format)");
     string::size_type n=InputOLD_File.find_last_of("/");
     runpath=InputOLD_File.substr(0,n);
     read_input_old_format=true;
     start_GUI=false;
  }
  else if (RunType_List->get<string>("R1: Run Type")=="No GUI / XML File"){
     InputXML_File=RunType_List->get<string>("R2: Input File (XML)");
     string::size_type n=InputXML_File.find_last_of("/");
     runpath=InputXML_File.substr(0,n);
     read_xml_file=true;
     start_GUI=false;
  }
  
           /* Create the empty parameter list called Tramonto_List. *
            * All Tramonto parameters will be added to this list.   */

  RCP<ParameterList> Tramonto_List = RCP<ParameterList>((new ParameterList("Root Tramonto List")));


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

  RCP<ParameterList> Geometry_List = sublist(Surface_List, "S4: Surface Geometry");
  RCP<ParameterList> SurfGeom0_List = sublist(Geometry_List, "SurfType 1 : Geometry");
  RCP<ParameterList> SurfGeom1_List = sublist(Geometry_List, "SurfType 2 : Geometry");
  RCP<ParameterList> SurfGeom2_List = sublist(Geometry_List, "SurfType 3 : Geometry");
  RCP<ParameterList> SurfGeom3_List = sublist(Geometry_List, "SurfType 4 : Geometry");
  RCP<ParameterList> SurfGeom4_List = sublist(Geometry_List, "SurfType 5 : Geometry");
  RCP<ParameterList> SurfGeom5_List = sublist(Geometry_List, "SurfType 6 : Geometry");
  RCP<ParameterList> SurfGeom6_List = sublist(Geometry_List, "SurfType 7 : Geometry");
  RCP<ParameterList> SurfGeom7_List = sublist(Geometry_List, "SurfType 8 : Geometry");
  RCP<ParameterList> SurfGeom8_List = sublist(Geometry_List, "SurfType 9 : Geometry");
  RCP<ParameterList> SurfGeom9_List = sublist(Geometry_List, "SurfType 10: Geometry");
  RCP<ParameterList> SurfGeom10_List = sublist(Geometry_List, "SurfType 11: Geometry");
  RCP<ParameterList> SurfGeom11_List = sublist(Geometry_List, "SurfType 12: Geometry");

  RCP<ParameterList> PotentialsWW_List = sublist(Surface_List, "S5: Surface Interactions");
  RCP<ParameterList> SurfaceParamCharge_List = sublist(Surface_List, "S6: Charged Surfaces and Sources");
  RCP<ParameterList> SurfaceInteraction_List = sublist(Surface_List, "S7: Fluid-Surface Interaction Types");
  RCP<ParameterList> PotentialsWF_List = sublist(Surface_List, "S7: Fluid-Surface Potential Parameterss");

  RCP<ParameterList> Continuation_List = sublist(Tramonto_List, "Sect. 6: Continuation");
  RCP<ParameterList> DensProfile_List = sublist(Tramonto_List, "Sect. 7: Initial Guess Options");
  RCP<ParameterList> Output_List = sublist(Tramonto_List, "Sect. 8: Output Control");
  RCP<ParameterList> Solver_List = sublist(Tramonto_List, "Sect. 9: Numerical Methods");

  RCP<ParameterList> NonlinearSolver_List = sublist(Solver_List, "SL1: Nonlinear Solver Options");
  RCP<ParameterList> LinearSolver_List = sublist(Solver_List, "SL2: Linear Solver Options");
  RCP<ParameterList> PhysicsMethod_List = sublist(Solver_List, "SL3: Physics based Solver Options");
  RCP<ParameterList> Coarsening_List = sublist(Solver_List, "SL4: Numerical Coarsening Options");
  RCP<ParameterList> LoadBalance_List = sublist(Solver_List, "SL5: Load Balancing Options");

  if (read_xml_file){

  /**
   * Then we just call getInput! There's a little more to it, so let's
   * head on over to the inputs.xml file to see what's going on there.
   */
        Optika::getInput(RunType_List->get<string>("R2: Input File (XML)"), Tramonto_List);
cout << "done with GUI - heading over to Tramonto data format translation....."<<endl;
        dft_GUI_toTramonto(Tramonto_List,Mesh_List,Functional_List,Fluid_List,
                     PotentialsFF_List,Polymer_List,PolymerGraft_List,PolymerArch_List,PolymerCMS_List,
                     StatePoint_List,Diffusion_List,ChargedFluid_List,Continuation_List,
                     Solver_List,Coarsening_List,LoadBalance_List,PhysicsMethod_List,LinearSolver_List,NonlinearSolver_List,Output_List,
                     DensProfile_List,Surface_List,SurfGeom0_List, SurfGeom1_List, SurfGeom2_List, SurfGeom3_List, SurfGeom4_List, SurfGeom5_List,
                     SurfGeom6_List, SurfGeom7_List, SurfGeom8_List, SurfGeom9_List, SurfGeom10_List,SurfGeom11_List);
  }
  else{
            /* Create a  dependency sheet for the GUI*/
     RCP<DependencySheet> depSheet_Tramonto = rcp(new DependencySheet());

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
     dft_GUI_potentialsWW_set_defaults(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,PotentialsFF_List,PotentialsWW_List,SurfaceParamCharge_List);
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
        dft_GUI_potentialsWW_set_OldFormat(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,PotentialsFF_List,PotentialsWW_List,SurfaceParamCharge_List);
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
     dft_GUI_potentialsWW_dependencies(Tramonto_List,depSheet_Tramonto,Mesh_List,Surface_List,PotentialsFF_List,PotentialsWW_List,SurfaceParamCharge_List);
     dft_GUI_Continuation_dependencies(Tramonto_List,depSheet_Tramonto,PotentialsFF_List,Continuation_List);
     dft_GUI_NumericalMethods_dependencies(Tramonto_List,depSheet_Tramonto,Functional_List,Fluid_List,Polymer_List,Solver_List,
                              Coarsening_List,LoadBalance_List,PhysicsMethod_List,NonlinearSolver_List,LinearSolver_List);

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
        if (start_GUI){
            Optika::getInput(Tramonto_List,depSheet_Tramonto);
            dft_GUI_toTramonto(Tramonto_List,Mesh_List,Functional_List,Fluid_List,
                     PotentialsFF_List,Polymer_List,PolymerGraft_List,PolymerArch_List,PolymerCMS_List,
                     StatePoint_List,Diffusion_List,ChargedFluid_List,Continuation_List,
                     Solver_List,Coarsening_List,LoadBalance_List,PhysicsMethod_List,LinearSolver_List,NonlinearSolver_List,Output_List,
                     DensProfile_List,Surface_List,SurfGeom0_List, SurfGeom1_List, SurfGeom2_List, SurfGeom3_List, SurfGeom4_List, SurfGeom5_List,
                     SurfGeom6_List, SurfGeom7_List, SurfGeom8_List, SurfGeom9_List, SurfGeom10_List,SurfGeom11_List);
        }
  }
  
  /* Print out what the user entered in nice XML format.  */
/*  RCP<FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  writeParameterListToXmlOStream(*Tramonto_List, *out);*/

  return;
}
/*********************************************************************************************************************/
/*int test_NsurfTypes(int Ninput,int Ntest){

   if (Ninput >= Ntest) return Ntest;
   else return -1;
}*/
/*********************************************************************************************************************/

