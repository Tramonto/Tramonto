using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_surfaces_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List)
{
  /* VALIDATORS*/
  RCP<EnhancedNumberValidator<int> > NwallValidator = rcp(new EnhancedNumberValidator<int>(0,1000,1));
  RCP<EnhancedNumberValidator<int> > NwallTypeValidator = rcp(new EnhancedNumberValidator<int>(0,12,1));
  RCP<FileNameValidator> SurfaceFile_Validator = rcp(new FileNameValidator);
  RCP<StringValidator> SurfacePos_Validator = rcp<StringValidator>(new StringValidator(tuple<std::string>(
        "Read From File","Set up in GUI","Random placement of walls")));


  /* Set Surface PARAMETERS */
  Surface_List->set("S1: Number of Surfaces", 0, "Number of surfaces (or total surface subunits) in the problem",NwallValidator);
  Surface_List->set("S2: Number of macro surfaces",  Surface_List->get<int>("S1: Number of Surfaces"), "Indicates if groups of surface form macrosurfaces.\n For example atoms that are part of a single molecule would have Nsurf_macro=1 and Nsurf>1)",NwallValidator);
  Surface_List->set("S3: Number of surface types", 1, "Number of different types surfaces (or total surface subunits) in the problem",NwallTypeValidator);

  SurfacePosition_List->set("SP0: Type of surface position/charge entry", "Read From File", "Identify how the positions and charges of individual files will be set up",SurfacePos_Validator);
  SurfacePosition_List->set("SP1: file containing surface positions and charges","dft_surfaces.dat","Enter the file name where the polymer architecture can be found.", SurfaceFile_Validator);

  Array<int> SurfaceType_Array(Surface_List->get<int>("S1: Number of Surfaces"),0);
  SurfacePosition_List->set("SP2: Surface Type IDs", SurfaceType_Array, "Type ID for each surface in the system.  Surface type IDs identify geometry and properties of surfaces.  Many surfaces may have the same type ID. Indexing starts with 0.");

  Array<int> SurfaceMacro_Array(Surface_List->get<int>("S1: Number of Surfaces"),0);
  SurfacePosition_List->set("SP3: Surface Macro IDs", SurfaceMacro_Array, "Macrosurface ID for each surface in the system.  Surfaces with the same macrosurface ID will be treated as one larger macrosurfaces for computation of forces.");

  TwoDArray<double> SurfacePos_Array( (Surface_List->get<int>("S1: Number of Surfaces"))*(Mesh_List->get<int>("M1_Ndim")),0.0);
  SurfacePosition_List->set("SP4: Surface Positions", SurfacePos_Array, "Position of each surface in the system in the coordinate system with the origin at the center of the computational domain (-Size_x[idim]/2,+Size_x[idim]/2).");

   return;
}
/*******************************************************************************************/
void dft_GUI_surfaces_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List)
{
  int i,j;
  int tmp_1D[NWALL_MAX];
  /* VALIDATORS*/
  RCP<EnhancedNumberValidator<int> > NwallValidator = rcp(new EnhancedNumberValidator<int>(0,1000,1));
  RCP<EnhancedNumberValidator<int> > NwallTypeValidator = rcp(new EnhancedNumberValidator<int>(0,12,1));
  RCP<FileNameValidator> SurfaceFile_Validator = rcp(new FileNameValidator);
  RCP<StringValidator> SurfacePos_Validator = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                                      "Read From File","Set up in GUI","Random placement of walls")));

  Surface_List->set("S1: Number of Surfaces", Nwall, "Number of surfaces (or total surface subunits) in the problem",NwallValidator);
  Surface_List->set("S2: Number of macro surfaces",  Nlink, "Indicates if groups of surface form macrosurfaces.\n For example atoms that are part of a single molecule would have Nsurf_macro=1 and Nsurf>1)",NwallValidator);
  Surface_List->set("S3: Number of surface types", Nwall_type, "Number of different types surfaces (or total surface subunits) in the problem",NwallTypeValidator);

  SurfacePosition_List->set("SP0: Type of surface position/charge entry", "Set up in GUI", "Identify how the positions and charges of individual files will be set up",SurfacePos_Validator);
  SurfacePosition_List->set("SP1: file containing surface positions and charges","dft_surfaces.dat","Enter the file name where the polymer architecture can be found.", SurfaceFile_Validator);

  for (i=0;i<Nwall;i++) tmp_1D[i]=WallType[i]; 
  Array<int> SurfaceType_Array(tmp_1D, tmp_1D+Nwall);
  SurfacePosition_List->set("SP2: Surface Type IDs", SurfaceType_Array, "Type ID for each surface in the system.  Surface type IDs identify geometry and properties of surfaces.  Many surfaces may have the same type ID. Indexing starts with 0.");

  for (i=0;i<Nwall;i++) tmp_1D[i]=Link[i]; 
  Array<int> SurfaceMacro_Array(tmp_1D, tmp_1D+Nwall);
  SurfacePosition_List->set("SP3: Surface Macro IDs", SurfaceMacro_Array, "Macrosurface ID for each surface in the system.  Surfaces with the same macrosurface ID will be treated as one larger macrosurfaces for computation of forces.");

  TwoDArray<double> SurfacePos_Array(Nwall,Ndim);
  for (i=0;i<Nwall;i++) for (j=0;j<Ndim;j++) SurfacePos_Array[i][j]=WallPos[j][i]; 
  SurfacePosition_List->set("SP4: Surface Positions", SurfacePos_Array, "Position of each surface in the system in the coordinate system with the origin at the center of the computational domain (-Size_x[idim]/2,+Size_x[idim]/2).");
  return;
}
/*******************************************************************************************/
void dft_GUI_surfaces_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List)
{
  /* DEPENDENCIES */
  RCP<NumberVisualDependency<int> > NsurfType_Dep = rcp(
        new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), Surface_List->getEntryRCP("S3: Number of surface types")));

  RCP<NumberVisualDependency<int> > NsurfMacro_Dep = rcp(
        new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), Surface_List->getEntryRCP("S2: Number of macro surfaces")));

  Dependency::ParameterEntryList ShowSurfacePosition_Deps;
  ShowSurfacePosition_Deps.insert(SurfacePosition_List->getEntryRCP("SP3: Surface Macro IDs"));
  ShowSurfacePosition_Deps.insert(SurfacePosition_List->getEntryRCP("SP2: Surface Type IDs"));

  RCP<StringVisualDependency> SurfType_Dep=rcp(new StringVisualDependency(SurfacePosition_List->getEntryRCP("SP0: Type of surface position/charge entry"),
                                              ShowSurfacePosition_Deps,tuple<std::string>("Set up in GUI","Random placement of walls")));

  RCP<StringVisualDependency> SurfPos_Dep=rcp(new StringVisualDependency(SurfacePosition_List->getEntryRCP("SP0: Type of surface position/charge entry"),
                                              SurfacePosition_List->getEntryRCP("SP4: Surface Positions"),"Set up in GUI"));

  RCP<StringVisualDependency> SurfFile_Dep=rcp(new StringVisualDependency(SurfacePosition_List->getEntryRCP("SP0: Type of surface position/charge entry"),
         SurfacePosition_List->getEntryRCP("SP1: file containing surface positions and charges"),"Read From File"));

  RCP<TwoDRowDependency<int,double> > SurfPosRowNumber_Dep = rcp(
                new TwoDRowDependency<int,double>(Surface_List->getEntryRCP("S1: Number of Surfaces"),SurfacePosition_List->getEntryRCP("SP4: Surface Positions")));

  RCP<TwoDColDependency<int,double> > SurfPosColNumber_Dep = rcp(
                new TwoDColDependency<int,double>(Mesh_List->getEntryRCP("M1_Ndim"),SurfacePosition_List->getEntryRCP("SP4: Surface Positions")));

  RCP<NumberArrayLengthDependency<int,int> > SurfTypeIDLength_Dep = rcp(
             new NumberArrayLengthDependency<int,int>(Surface_List->getEntryRCP("S1: Number of Surfaces"), SurfacePosition_List->getEntryRCP("SP2: Surface Type IDs")));

  RCP<NumberArrayLengthDependency<int,int> > SurfMacroIDLength_Dep = rcp(
             new NumberArrayLengthDependency<int,int>(Surface_List->getEntryRCP("S1: Number of Surfaces"), SurfacePosition_List->getEntryRCP("SP3: Surface Macro IDs")));


  /* DEPENDENCY SHEET ENTRIES*/
  depSheet_Tramonto->addDependency(NsurfType_Dep);
  depSheet_Tramonto->addDependency(NsurfMacro_Dep);
  depSheet_Tramonto->addDependency(SurfPos_Dep);
  depSheet_Tramonto->addDependency(SurfFile_Dep);
  depSheet_Tramonto->addDependency(SurfPosRowNumber_Dep);
  depSheet_Tramonto->addDependency(SurfPosColNumber_Dep);
  depSheet_Tramonto->addDependency(SurfTypeIDLength_Dep);
  depSheet_Tramonto->addDependency(SurfMacroIDLength_Dep);
  
  return;
}
/*******************************************************************************************/

