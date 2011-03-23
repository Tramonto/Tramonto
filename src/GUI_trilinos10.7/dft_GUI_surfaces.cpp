using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_surfaces(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceGeometry_List,
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


        RCP<EnhancedNumberValidator<int> > NwallValidator = rcp(new EnhancedNumberValidator<int>(0,1000,1));

        /* 1D options */
        RCP<StringValidator> surfTypeVali1 = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "PLANE: Infinite in two dimensions")));
        RCP<ArrayStringValidator> arraySurfTypeVali1 = rcp(new ArrayStringValidator(surfTypeVali1));

        /* 2D options */
       RCP<StringValidator> surfTypeVali2 = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "PLANE: Infinite in two dimensions",
                       "PLANE: Infinite in one dimension",
                       "CYLINDERICAL SURFACE: Infinite_L",
                       "CYLINDRICAL PORE: Infinite_L",
                       "SLIT PORE: Finite_L",
                       "SLIT TAPERED PORE: Finite_L")));
        RCP<ArrayStringValidator> arraySurfTypeVali2 = rcp(new ArrayStringValidator(surfTypeVali2));

        /* 3D options */
       RCP<StringValidator> surfTypeVali3 = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "PLANE: Infinite in two dimensions",
                       "PLANE: Infinite in one dimensions",
                       "PLANE: Finite three dimensions",
                       "SPHERICAL SURFACE: COLLOID",
                       "SPHERICAL SURFACE: ATOM",
                       "SPHERICAL SURFACE: POINT ATOM",
                       "SPHERICAL CAVITY",
                       "CYLINDERICAL SURFACE: Infinite_L",
                       "CYLINDERICAL SURFACE: Finite_L",
                       "CYLINDERICAL SURFACE PERIODIC: Infinite_L",
                       "CYLINDRICAL PORE: Infinite_L",
                       "CYLINDRICAL PORE: Finite_L",
                       "CYLINDRICAL TAPERED PORE: Finite_L")));
        RCP<ArrayStringValidator> arraySurfTypeVali3 = rcp(new ArrayStringValidator(surfTypeVali3));

     /* INDEPENDENT PARAMEBERS */

       Surface_List->set("S1: Number of Surfaces", 0, "Number of surfaces (or total surface subunits) in the problem",NwallValidator);
     /* DEPENDENT PARAMETERS */
        Surface_List->set("S2: Number of macro surfaces",  Surface_List->get("S1: Number of Surfaces", (int)NULL), "Indicates if groups of surface form macrosurfaces.\n For example atoms that are part of a single molecule would have Nsurf_macro=1 and Nsurf>1)",NwallValidator);
        Surface_List->set("S3: Number of surface types", 1, "Number of different types surfaces (or total surface subunits) in the problem",NwallValidator);

        Array<std::string> surfTypeArray(Surface_List->get("S3: Number of surface types", (int)NULL),"PLANE: Infinite in two dimensions");
        switch(Mesh_List->get("M1_Ndim", (int)NULL)){
             case 1:
                SurfaceGeometry_List->set("SG1: Surface Type Array", surfTypeArray, "Select types of surfaces in the problem.", arraySurfTypeVali1); break;
             case 2:
                SurfaceGeometry_List->set("SG1: Surface Type Array", surfTypeArray, "Select types of surfaces in the problem.", arraySurfTypeVali2); break;
             case 3:
                SurfaceGeometry_List->set("SG1: Surface Type Array", surfTypeArray, "Select types of surfaces in the problem.", arraySurfTypeVali3); break;
             default:
                SurfaceGeometry_List->set("SG1: Surface Type Array", surfTypeArray, "Select types of surfaces in the problem.", arraySurfTypeVali1); break;
        }

/* need to figure out what to do about orientation and wall param arrays.  Currently we have an array of length surf_type for each param.
   we then enter parameters that matter to the particular surface type. In the optika world it would be great to make this a lot better
   by only showing options that matter.  for example some surface type need an orientation while others do not.  We would like to 
    probe the surface type array change the possible entries in the wall parameter arrays depending on surface type. Note that some array 
    entries are not even needed*/

/* as an alternative, could we set up an array, but then use the set command on one entry of the array at a time so that we could have
   different options for each of the array entries? */


/* maybe the right thing to do is just set up parameters like RADIUS, WIDTH, ORIENTATION, LENGTH and only show the parameters that are needed for
   each wall entry */


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
        RCP<NumberVisualDependency<int> > NsurfType_Dep = rcp(
           new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), Surface_List->getEntryRCP("S3: Number of surface types")));

        RCP<NumberVisualDependency<int> > NsurfMacro_Dep = rcp(
           new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), Surface_List->getEntryRCP("S2: Number of macro surfaces")));

        RCP<NumberArrayLengthDependency<int, std::string> > surfTypeLengthDep = rcp(
           new NumberArrayLengthDependency<int, std::string>( 
               Surface_List->getEntryRCP("S3: Number of surface types"), 
               SurfaceGeometry_List->getEntryRCP("SG1: Surface Type Array"))
        );

        RangeValidatorDependency<int>::RangeToValidatorMap dimranges;
        dimranges[std::pair<int,int>(1,1)] = arraySurfTypeVali1;
        dimranges[std::pair<int,int>(2,2)] = arraySurfTypeVali2;
        dimranges[std::pair<int,int>(3,3)] = arraySurfTypeVali3;

        RCP<RangeValidatorDependency<int> >
        surfTypeValiDep = rcp(
                new RangeValidatorDependency<int>( 
                    Mesh_List->getEntryRCP("M1_Ndim"), 
                    SurfaceGeometry_List->getEntryRCP("SG1: Surface Type Array"), dimranges)
        );

        /* need to write a dependency for "WW1: Compute wall-wall interactions?" to show up if Ndim=3D and at least one surface type is atomic */
        RCP<BoolVisualDependency> UWWType_Dep = rcp(new BoolVisualDependency(PotentialsWW_List->getEntryRCP("WW1: Compute wall-wall interactions?"),  
			PotentialsWW_List->getEntryRCP("WW2: Type of wall-wall interactions"), true));



     /* DEPENDENCY SHEET ENTRIES*/
        depSheet_Tramonto->addDependency(NsurfType_Dep);
        depSheet_Tramonto->addDependency(NsurfMacro_Dep);
        depSheet_Tramonto->addDependency(surfTypeLengthDep);
        depSheet_Tramonto->addDependency(surfTypeValiDep);

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

