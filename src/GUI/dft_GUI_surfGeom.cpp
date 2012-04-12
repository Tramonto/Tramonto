using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_surface_geometry(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfGeom_List)
{

  /****************************************************************************************************************/
  /********************************* SURFACE DEFINITIONS SECTION **************************************************/
  /****************************************************************************************************************/

     /* VALIDATORS*/

         RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));
         RCP<EnhancedNumberValidator<int> > NumFuncValidator = rcp(new EnhancedNumberValidator<int>(0,4,1));

        /* 1D options */
        RCP<StringValidator> surfTypeVali1 = rcp(new StringValidator(tuple<std::string>(
                       "PLANE: 1D-2D-3D :Infinite in two dimensions","1-ELEMENT SURFACE: 1D-2D-3D :Membrane,Rod,Point")));
/*        RCP<ArrayStringValidator> arraySurfTypeVali1 = rcp(new ArrayStringValidator(surfTypeVali1));*/

        /* 2D options */
       RCP<StringValidator> surfTypeVali2 = rcp(new StringValidator(tuple<std::string>(
                       "PLANE: 1D-2D-3D :Infinite in two dimensions",
                       "PLANE: 2D-3D :Finite length all dimensions",
                       "CYLINDER: 2D :Infinite length",
                       "1-ELEMENT SURFACE: 1D-2D-3D :Membrane,Rod,Point",
                       "PORE CYLINDRICAL: 2D :Infinite Length",
                       "PORE SLIT: 2D :Finite Length")));
/*        RCP<ArrayStringValidator> arraySurfTypeVali2 = rcp(new ArrayStringValidator(surfTypeVali2));*/

        /* 3D options */
       RCP<StringValidator> surfTypeVali3 = rcp(new StringValidator(tuple<std::string>(
                       "PLANE: 1D-2D-3D :Infinite in two dimensions",
                       "PLANE: 2D-3D :Finite length all dimensions",
                       "CYLINDER: 2D :Infinite length",
                       "CYLINDER: 3D :Finite Length",
                       "SPHERE/COLLOID: 3D :Radius>Sigma_w typically",
                       "ATOMIC: 3D :Radius=Sigma_w",
                       "1-ELEMENT SURFACE: 1D-2D-3D :Membrane,Rod,Point",
                       "PORE CYLINDRICAL: 2D :Infinite Length",
                       "PORE SPHERICAL: 3D :",
                       "PORE SLIT: 2D :Finite Length",
                       "PORE CYLINDRICAL: 3D :Finite Length")));
/*        RCP<ArrayStringValidator> arraySurfTypeVali3 = rcp(new ArrayStringValidator(surfTypeVali3));*/

     /* INDEPENDENT PARAMEBERS */

        switch(Mesh_List->get("M1_Ndim", (int)NULL)){    /* note that all are set to Vali3 right now because the RangeValidatorDependency isn't working right */
             case 1: SurfGeom_List->set("SG1: Surface Type", "PLANE: 1D-2D-3D :Infinite in two dimensions", "Set base surface type", surfTypeVali3); break;
             case 2: SurfGeom_List->set("SG1: Surface Type", "PLANE: 1D-2D-3D :Infinite in two dimensions", "Set base surface type", surfTypeVali3); break;
             case 3: SurfGeom_List->set("SG1: Surface Type", "PLANE: 1D-2D-3D :Infinite in two dimensions", "Set base surface type", surfTypeVali3); break;
             default: SurfGeom_List->set("SG1: Surface Type", "PLANE: 1D-2D-3D :Infinite in two dimensions", "Set base surface type", surfTypeVali3); break;
        }

        SurfGeom_List->set("SG2: Orientation Plane",0,"Enter the direction of the normal to the planar base surface",DimValidator);
        SurfGeom_List->set("SG2: Orientation Cylinder",0,"Enter the direction of the long axis of the cylindrical base surface",DimValidator);

        SurfGeom_List->set("SG3: Radius",1.0,"Enter the base radius of the cylindrical or spherical surface");
        SurfGeom_List->set("SG3: Half-width",1.0,"Enter the distance from the center to the edge of the surface.");
        Array<double>HalfWidth_Array(Mesh_List->get<int>("M1_Ndim"),1.0);
        SurfGeom_List->set("SG3: Half-width array",HalfWidth_Array,"Enter the distance from the center to the edge of the surface in all Ndim.");
        SurfGeom_List->set("SG3: Half-Length",1.0,"Enter the distance from the center to the edge of the surface down long axis of cylinder.");
       
           /* surface roughness */ 
        SurfGeom_List->set("SG4 : Apply Roughness to Radius?",false,"Set to true for Radius=R+tau where tau is a random roughness parameter.");
        SurfGeom_List->set("SG4 : Apply Roughness to HalfLength?",false,"Set to true for HalfLength=HalfLength+tau where tau is a random roughness parameter to the length parameter.");
        SurfGeom_List->set("SG4 : Apply Roughness to HalfWidth?",false,"Set to true for HalfWidth=HalfWidth+tau where tau is a random roughness parameter to the length parameter.");

        SurfGeom_List->set("SG4.1: Tau_max",0.0,"Enter maximum value for random roughness.  Surface param (L or R) will become (L or R) +/- tau where tau<=Tau_max");
        SurfGeom_List->set("SG4.2: TauLength",0.0,"Enter characteristic lateral length scale for roughness.  This allows comparison of roughness on different length scales.");


           /* periodic function modifications to the surface  */ 
        SurfGeom_List->set("SG5 :Add Periodic Func. to HalfWidth?",false,"Set to true to superimpose one or more periodic function(s) on the surface.  HalfWidth=HalfWidth+A*cos(2 pi (x-x_o)/lambda.");
        SurfGeom_List->set("SG5 :Add Periodic Func. to HalfLength?",false,"Set to true to superimpose one or more periodic function(s) on the surface.  HalfLength=HalfLength+A*cos(2 pi (x-x_o)/lambda.");

        SurfGeom_List->set("SG5.1: Number of Periodic Fncs.",0,"Set number of periodic functions that will be superimposed on this surface.",NumFuncValidator);
           Array<int>OrientationPeriodic_Array(SurfGeom_List->get<int>("SG5.1: Number of Periodic Fncs."),0);
        SurfGeom_List->set("SG5.2: Orientation Periodic Fnc.",OrientationPeriodic_Array,"Identify the Orientation of each of the periodic functions of interest. ");
           Array<double>AmplitudePeriodic_Array(SurfGeom_List->get<int>("SG5.1: Number of Periodic Fncs."),0.0);
        SurfGeom_List->set("SG5.3: Amplitude Periodic Fnc.",AmplitudePeriodic_Array,"Set the amplitude of each of the periodic functions of interest.");
           Array<double>PeriodPeriodic_Array(SurfGeom_List->get<int>("SG5.1: Number of Periodic Fncs."),0.0);
        SurfGeom_List->set("SG5.4: Period of Periodic Fnc.",PeriodPeriodic_Array,"Set the period of each of the periodic functions of interest.");
           Array<double>OriginPeriodic_Array(SurfGeom_List->get<int>("SG5.1: Number of Periodic Fncs."),0.0);
        SurfGeom_List->set("SG5.5: Origin of Periodic Fnc.",OriginPeriodic_Array,"Set the origin in the defined Orientation direction for each of the periodic functions of interest.");


           /* linear function modifications to the surface  */ 
        SurfGeom_List->set("SG6 :Add Linear Func. to HalfWidth?",false,"Set to true to superimpose one or more linear function(s) on the surface.  HalfWidth=HalfWidth+A*(x-x_o).");
        SurfGeom_List->set("SG6 :Add Linear Func. to HalfLength?",false,"Set to true to superimpose one or more linear function(s) on the surface.  HalfLength=HalfLength+A*(x-x_o).");

        SurfGeom_List->set("SG6.1: Number of Linear Fncs.",0,"Set number of linear functions that will be superimposed on this surface.",NumFuncValidator);
           Array<int>OrientationLinear_Array(SurfGeom_List->get<int>("SG6.1: Number of Linear Fncs."),0);
        SurfGeom_List->set("SG6.2: Orientation Linear Fnc.",OrientationLinear_Array,"Identify the Orientation of each of the linear functions of interest. ");
           Array<double>SlopeLinear_Array(SurfGeom_List->get<int>("SG6.1: Number of Linear Fncs."),0.0);
        SurfGeom_List->set("SG6.3: Slope of Linear Fnc.",SlopeLinear_Array,"Set the slope for each of the linear functions of interest.");
           Array<double>OriginLinear_Array(SurfGeom_List->get<int>("SG6.1: Number of Linear Fncs."),0.0);
        SurfGeom_List->set("SG6.4: Origin of Linear Fnc.",OriginLinear_Array,"Set the origin of each of the linear functions of interest.");
           Array<double>EndptLinear_Array(SurfGeom_List->get<int>("SG6.1: Number of Linear Fncs."),0.0);
        SurfGeom_List->set("SG6.5: End Point for Linear Fnc.",EndptLinear_Array,"Set an end position for the linear function.  Can use this to create divots in smooth surfaces because the function will only be applied between the origin and the end point.");

        SurfGeom_List->set("SG7 : Angular cutouts from the surface?",false,"Set to true if you want to create a surface with missing wedge. This can be helpful for creating heterogeneous surfaces of various kinds");
        SurfGeom_List->set("SG7.1: Angle start",0.0,"Enter the angle where the cutout should begin in degrees.  The origin of the angular cutout is always taken to be the defined position of center of the surface. The angles are measured from the x axis in 2D systems and an appropriate permutation is applied for cylinders in 3D depending on the orientation of the cylinder.  For example a cylinder with orientation=2 (long axis parallel to z axis) will have angles measured from the x axis of the x-y plane.");
        SurfGeom_List->set("SG7.2: Angle end",0.0,"Enter the angle where the cutout should end in degrees.  Same conventions as the Angle start parameter.");


     /* DEPENDENCIES */
        RangeValidatorDependency<int>::RangeToValidatorMap dimranges;
        dimranges[std::pair<int,int>(1,1)] = surfTypeVali1;
        dimranges[std::pair<int,int>(2,2)] = surfTypeVali2;
        dimranges[std::pair<int,int>(3,3)] = surfTypeVali3;

        RCP<RangeValidatorDependency<int> >
        surfTypeValiDep = rcp(
                new RangeValidatorDependency<int>(Mesh_List->getEntryRCP("M1_Ndim"), SurfGeom_List->getEntryRCP("SG1: Surface Type"), dimranges, surfTypeVali3)
        );

       RCP<StringVisualDependency> OrientPlane_Dep = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG2: Orientation Plane"),
                tuple<std::string>("PLANE: 1D-2D-3D :Infinite in two dimensions")));

       RCP<StringVisualDependency> OrientCyl_Dep = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG2: Orientation Cylinder"),
                tuple<std::string>("CYLINDER: 3D :Finite Length", "PORE CYLINDRICAL: 2D :Infinite Length", "PORE SPHERICAL: 3D :",
                       "PORE SLIT: 2D :Finite Length", "PORE CYLINDRICAL: 3D :Finite Length")));

       RCP<StringVisualDependency> Radius_Dep = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG3: Radius"),
                tuple<std::string>("CYLINDER: 2D :Infinite length","CYLINDER: 3D :Finite Length",
                       "SPHERE/COLLOID: 3D :Radius>Sigma_w typically","PORE CYLINDRICAL: 2D :Infinite Length",
                       "PORE SPHERICAL: 3D :","PORE SLIT: 2D :Finite Length","PORE CYLINDRICAL: 3D :Finite Length")));

       RCP<StringVisualDependency> HalfWidth_Dep = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG3: Half-width"),
                tuple<std::string>("PLANE: 1D-2D-3D :Infinite in two dimensions")));

       RCP<StringVisualDependency> HalfLength_Dep = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG3: Half-Length"),
                tuple<std::string>("CYLINDER: 3D :Finite Length","PORE SLIT: 2D :Finite Length","PORE CYLINDRICAL: 3D :Finite Length")));

       RCP<StringVisualDependency> HalfWidthArray_Dep = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG3: Half-width array"),
                tuple<std::string>("PLANE: 2D-3D :Finite length all dimensions")));

       RCP<NumberArrayLengthDependency<int,double> > HalfWidthArrayLength_Dep = rcp(
               new NumberArrayLengthDependency<int,double>(Mesh_List->getEntryRCP("M1_Ndim"), SurfGeom_List->getEntryRCP("SG3: Half-width array")));

       RCP<StringVisualDependency> Rough_Dep = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG4 : Apply Roughness to Radius?"),
                tuple<std::string>("CYLINDER: 3D :Finite Length","CYLINDER: 2D :Infinite length",
                       "PORE CYLINDRICAL: 2D :Infinite Length", "PORE SLIT: 2D :Finite Length", "PORE CYLINDRICAL: 3D :Finite Length")));

       RCP<StringVisualDependency> Rough_Dep2 = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG4 : Apply Roughness to HalfLength?"),
                tuple<std::string>("CYLINDER: 3D :Finite Length", "PORE SLIT: 2D :Finite Length", "PORE CYLINDRICAL: 3D :Finite Length"))); 

       RCP<StringVisualDependency> Rough_Dep3 = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG4 : Apply Roughness to HalfWidth?"),
                tuple<std::string>("PLANE: 1D-2D-3D :Infinite in two dimensions","PLANE: 2D-3D :Finite length all dimensions")));

        RCP<BoolCondition> RoughTFCond1=rcp(new BoolCondition(SurfGeom_List->getEntryRCP("SG4 : Apply Roughness to Radius?")));
        RCP<BoolCondition> RoughTFCond2=rcp(new BoolCondition(SurfGeom_List->getEntryRCP("SG4 : Apply Roughness to HalfLength?")));
        RCP<BoolCondition> RoughTFCond3=rcp(new BoolCondition(SurfGeom_List->getEntryRCP("SG4 : Apply Roughness to HalfWidth?")));
        Condition::ConstConditionList RoughTF_ConList=tuple<RCP<const Condition> >(RoughTFCond1,RoughTFCond2,RoughTFCond3);
        RCP<OrCondition>RoughTF_orCon=rcp(new OrCondition(RoughTF_ConList));

        Dependency::ParameterEntryList RoughParams_Deps;
        RoughParams_Deps.insert(SurfGeom_List->getEntryRCP("SG4.1: Tau_max"));
        RoughParams_Deps.insert(SurfGeom_List->getEntryRCP("SG4.2: TauLength"));
        RCP<ConditionVisualDependency> RoughParams_Dep = rcp(new ConditionVisualDependency(RoughTF_orCon,RoughParams_Deps, true));


       RCP<StringVisualDependency> Periodic_Dep1 = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG5 :Add Periodic Func. to HalfLength?"),
                tuple<std::string>("CYLINDER: 3D :Finite Length"))); 

       RCP<StringVisualDependency> Periodic_Dep2 = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG5 :Add Periodic Func. to HalfWidth?"),
                tuple<std::string>("PLANE: 1D-2D-3D :Infinite in two dimensions","PLANE: 2D-3D :Finite length all dimensions")));

       RCP<BoolCondition> PeriodicTFCond1=rcp(new BoolCondition(SurfGeom_List->getEntryRCP("SG5 :Add Periodic Func. to HalfWidth?")));
       RCP<BoolCondition> PeriodicTFCond2=rcp(new BoolCondition(SurfGeom_List->getEntryRCP("SG5 :Add Periodic Func. to HalfLength?")));
       Condition::ConstConditionList PeriodicTF_ConList=tuple<RCP<const Condition> >(PeriodicTFCond1,PeriodicTFCond2);
       RCP<OrCondition>PeriodicTF_orCon=rcp(new OrCondition(PeriodicTF_ConList));

       Dependency::ParameterEntryList PeriodicParams_Deps;
       PeriodicParams_Deps.insert(SurfGeom_List->getEntryRCP("SG5.1: Number of Periodic Fncs."));
       PeriodicParams_Deps.insert(SurfGeom_List->getEntryRCP("SG5.2: Orientation Periodic Fnc."));
       PeriodicParams_Deps.insert(SurfGeom_List->getEntryRCP("SG5.3: Amplitude Periodic Fnc."));
       PeriodicParams_Deps.insert(SurfGeom_List->getEntryRCP("SG5.4: Period of Periodic Fnc."));
       PeriodicParams_Deps.insert(SurfGeom_List->getEntryRCP("SG5.5: Origin of Periodic Fnc."));
       RCP<ConditionVisualDependency> PeriodicParams_Dep = rcp(new ConditionVisualDependency(PeriodicTF_orCon,PeriodicParams_Deps, true));


       Dependency::ParameterEntryList PeriodicArrayLength_Deps;
       PeriodicArrayLength_Deps.insert(SurfGeom_List->getEntryRCP("SG5.3: Amplitude Periodic Fnc."));
       PeriodicArrayLength_Deps.insert(SurfGeom_List->getEntryRCP("SG5.4: Period of Periodic Fnc."));
       PeriodicArrayLength_Deps.insert(SurfGeom_List->getEntryRCP("SG5.5: Origin of Periodic Fnc."));
       RCP<NumberArrayLengthDependency<int,double> > PeriodicParamsArrayLength_Dep = rcp(
             new NumberArrayLengthDependency<int,double>(SurfGeom_List->getEntryRCP("SG5.1: Number of Periodic Fncs."), PeriodicArrayLength_Deps));
       RCP<NumberArrayLengthDependency<int,int> > PeriodicParamsArrayLength_Dep2 = rcp(
             new NumberArrayLengthDependency<int,int>(SurfGeom_List->getEntryRCP("SG5.1: Number of Periodic Fncs."), SurfGeom_List->getEntryRCP("SG5.2: Orientation Periodic Fnc.")));


       RCP<StringVisualDependency> Linear_Dep1 = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG6 :Add Linear Func. to HalfLength?"),
                tuple<std::string>("CYLINDER: 3D :Finite Length","PORE SLIT: 2D :Finite Length", "PORE CYLINDRICAL: 3D :Finite Length"))); 

        RCP<StringVisualDependency> Linear_Dep2 = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG6 :Add Linear Func. to HalfWidth?"),
                tuple<std::string>("PLANE: 1D-2D-3D :Infinite in two dimensions","PLANE: 2D-3D :Finite length all dimensions")));

        RCP<BoolCondition> LinearTFCond1=rcp(new BoolCondition(SurfGeom_List->getEntryRCP("SG6 :Add Linear Func. to HalfWidth?")));
        RCP<BoolCondition> LinearTFCond2=rcp(new BoolCondition(SurfGeom_List->getEntryRCP("SG6 :Add Linear Func. to HalfLength?")));
        Condition::ConstConditionList LinearTF_ConList=tuple<RCP<const Condition> >(LinearTFCond1,LinearTFCond2);
        RCP<OrCondition>LinearTF_orCon=rcp(new OrCondition(LinearTF_ConList));

        Dependency::ParameterEntryList LinearParams_Deps;
        LinearParams_Deps.insert(SurfGeom_List->getEntryRCP("SG6.1: Number of Linear Fncs."));
        LinearParams_Deps.insert(SurfGeom_List->getEntryRCP("SG6.2: Orientation Linear Fnc."));
        LinearParams_Deps.insert(SurfGeom_List->getEntryRCP("SG6.3: Slope of Linear Fnc."));
        LinearParams_Deps.insert(SurfGeom_List->getEntryRCP("SG6.4: Origin of Linear Fnc."));
        LinearParams_Deps.insert(SurfGeom_List->getEntryRCP("SG6.5: End Point for Linear Fnc."));
        RCP<ConditionVisualDependency> LinearParams_Dep = rcp(new ConditionVisualDependency(LinearTF_orCon,LinearParams_Deps, true));

       Dependency::ParameterEntryList LinearArrayLength_Deps;
       LinearArrayLength_Deps.insert(SurfGeom_List->getEntryRCP("SG6.3: Slope of Linear Fnc."));
       LinearArrayLength_Deps.insert(SurfGeom_List->getEntryRCP("SG6.4: Origin of Linear Fnc."));
       LinearArrayLength_Deps.insert(SurfGeom_List->getEntryRCP("SG6.5: End Point for Linear Fnc."));
       RCP<NumberArrayLengthDependency<int,double> > LinearParamsArrayLength_Dep = rcp(
             new NumberArrayLengthDependency<int,double>(SurfGeom_List->getEntryRCP("SG6.1: Number of Linear Fncs."), LinearArrayLength_Deps));
       RCP<NumberArrayLengthDependency<int,int> > LinearParamsArrayLength_Dep2 = rcp(
             new NumberArrayLengthDependency<int,int>(SurfGeom_List->getEntryRCP("SG6.1: Number of Linear Fncs."), SurfGeom_List->getEntryRCP("SG6.2: Orientation Linear Fnc.")));


        RCP<StringVisualDependency> Angle_Dep = rcp( new StringVisualDependency( 
                SurfGeom_List->getEntryRCP("SG1: Surface Type"), SurfGeom_List->getEntryRCP("SG7 : Angular cutouts from the surface?"),
                tuple<std::string>("CYLINDER: 3D :Finite Length","PORE SLIT: 2D :Finite Length", "PORE CYLINDRICAL: 3D :Finite Length",
                               "PORE CYLINDRICAL: 2D :Infinite Length",
                               "PLANE: 1D-2D-3D :Infinite in two dimensions","PLANE: 2D-3D :Finite length all dimensions",
                               "SPHERE/COLLOID: 3D :Radius>Sigma_w typically","CYLINDER: 2D :Infinite length","PORE SPHERICAL: 3D :"))); 

        Dependency::ParameterEntryList AngleParams_Deps;
        AngleParams_Deps.insert(SurfGeom_List->getEntryRCP("SG7.1: Angle start"));
        AngleParams_Deps.insert(SurfGeom_List->getEntryRCP("SG7.2: Angle end"));
        RCP<BoolVisualDependency> AngleParams_Dep = rcp(new BoolVisualDependency(SurfGeom_List->getEntryRCP("SG7 : Angular cutouts from the surface?"),AngleParams_Deps, true));



     /* DEPENDENCY SHEET ENTRIES*/
        depSheet_Tramonto->addDependency(surfTypeValiDep);
        depSheet_Tramonto->addDependency(OrientPlane_Dep);
        depSheet_Tramonto->addDependency(OrientCyl_Dep);
        depSheet_Tramonto->addDependency(Radius_Dep);
        depSheet_Tramonto->addDependency(HalfWidth_Dep);
        depSheet_Tramonto->addDependency(HalfLength_Dep);
        depSheet_Tramonto->addDependency(HalfWidthArray_Dep);
        depSheet_Tramonto->addDependency(HalfWidthArrayLength_Dep);

        depSheet_Tramonto->addDependency(Rough_Dep);
        depSheet_Tramonto->addDependency(Rough_Dep2);
        depSheet_Tramonto->addDependency(Rough_Dep3);
        depSheet_Tramonto->addDependency(RoughParams_Dep);
        depSheet_Tramonto->addDependency(Periodic_Dep1);
        depSheet_Tramonto->addDependency(Periodic_Dep2);
        depSheet_Tramonto->addDependency(PeriodicParams_Dep);
        depSheet_Tramonto->addDependency(PeriodicParamsArrayLength_Dep);
        depSheet_Tramonto->addDependency(PeriodicParamsArrayLength_Dep2);
        depSheet_Tramonto->addDependency(Linear_Dep1);
        depSheet_Tramonto->addDependency(Linear_Dep2);
        depSheet_Tramonto->addDependency(LinearParams_Dep);
        depSheet_Tramonto->addDependency(LinearParamsArrayLength_Dep);
        depSheet_Tramonto->addDependency(LinearParamsArrayLength_Dep2);
        depSheet_Tramonto->addDependency(Angle_Dep);
        depSheet_Tramonto->addDependency(AngleParams_Dep);


  return;
}

