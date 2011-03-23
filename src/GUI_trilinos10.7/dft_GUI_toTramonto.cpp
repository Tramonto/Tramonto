/*#include "Teuchos_ParameterList.hpp"*/
using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_toTramonto( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                        Teuchos::RCP<Teuchos::ParameterList> Mesh_List, 
                        Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                        Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                        Teuchos::RCP<Teuchos::ParameterList> SurfaceGeometry_List) 
{
  int i;
  /****************************** DIMENSION PARAMETER SECTION **********************************************************/
  /* this routine translates the parameters from the GUI to Tramonto.  */


    /***************************************/
    /* params from mesh section of the GUI */
    /***************************************/
    
    if (Mesh_List->get<bool>("M2_Dimensionless_Distance_Entry")) Length_ref=1.0; 
    else               Length_ref=Mesh_List->get<double>("M3_Reference_Length");

    Ndim=Mesh_List->get<int>("M1_Ndim");

    Array<double> M4_Size_x = Mesh_List->get<Array<double> >("M4_Size_x");
    for (i=0; i<Ndim; i++) Size_x[i]=M4_Size_x[i]; 


    Array<double> M5_Esize_x = Mesh_List->get<Array<double> >("M5_Esize_x");
    for (i=0; i<Ndim; i++) Esize_x[i]=M5_Esize_x[i]; 

    Array<string> BC_Xmin=Mesh_List->get<Array<string> >("M6_Boundary_Conditions_XMIN");
    for (i=0; i<Ndim; i++){
       if (BC_Xmin[i]=="Bulk") Type_bc[i][0]=IN_BULK;
       else if (BC_Xmin[i]=="In_Wall") Type_bc[i][0]=IN_WALL;
       else if (BC_Xmin[i]=="Reflect") Type_bc[i][0]=REFLECT;
       else if (BC_Xmin[i]=="Periodic") Type_bc[i][0]=PERIODIC;
       else if (BC_Xmin[i]=="Last_node") Type_bc[i][0]=LAST_NODE;
       else if (BC_Xmin[i]=="Last_node_restart") Type_bc[i][0]=LAST_NODE_RESTART;
    }

    Array<string> BC_Xmax=Mesh_List->get<Array<string> >("M7_Boundary_Conditions_XMAX");
    for (i=0; i<Ndim; i++){
       if (BC_Xmax[i]=="Bulk") Type_bc[i][1]=IN_BULK;
       else if (BC_Xmin[i]=="In_Wall") Type_bc[i][1]=IN_WALL;
       else if (BC_Xmin[i]=="Reflect") Type_bc[i][1]=REFLECT;
       else if (BC_Xmin[i]=="Periodic") Type_bc[i][1]=PERIODIC;
       else if (BC_Xmin[i]=="Last_node") Type_bc[i][1]=LAST_NODE;
       else if (BC_Xmin[i]=="Last_node_restart") Type_bc[i][1]=LAST_NODE_RESTART;
    }

    /*********************************************/
    /* params from functional section of the GUI */
    /*********************************************/
    if (Functional_List->get<string>("F0_Type_of_Calculation")=="Equilibrium (homogeneous boundary conditions)")          Type_interface=UNIFORM_INTERFACE;
    else if (Functional_List->get<string>("F0_Type_of_Calculation")=="Equilibrium (inhomogeneous boudary conditions)")    Type_interface=PHASE_INTERFACE;
    else if (Functional_List->get<string>("F0_Type_of_Calculation")=="Steady State Diffusion (inhomogeneous boundaries)") Type_interface=DIFFUSIVE_INTERFACE;

    if (Functional_List->get<string>("F1_HS_Functional")=="Ideal Fluid / No Volume Exclusion") Type_func=NONE;
    else if (Functional_List->get<string>("F1_HS_Functional")=="HS_Rosenfeld_Original")        Type_func=FMT1;
    else if (Functional_List->get<string>("F1_HS_Functional")=="HS_Rosenfeld_Zerocrossover")   Type_func=FMT2;
    else if (Functional_List->get<string>("F1_HS_Functional")=="HS_White_Bear")                Type_func=FMT3;
    else if (Functional_List->get<string>("F1_HS_Functional")=="HS_White_Bear2")               Type_func=FMT4;

    if (Functional_List->get<string>("F2_PAIRPOTcore_Functional")=="No Mean Field Functional")                            Type_attr=NONE;
    else if (Functional_List->get<string>("F2_PAIRPOTcore_Functional")=="UATT_core=[Umin; r=0,Rmin]")                     Type_attr=MFPAIR_RMIN_UMIN;
    else if (Functional_List->get<string>("F2_PAIRPOTcore_Functional")=="UATT_core=[0;r=0,Sigma]")                        Type_attr=MFPAIR_SIGMA_ZERO;
    else if (Functional_List->get<string>("F2_PAIRPOTcore_Functional")=="UATT_core=[U(sigma);r=0,Sigma]")                 Type_attr=MFPAIR_SIGMA_USIGMA;
    else if (Functional_List->get<string>("F2_PAIRPOTcore_Functional")=="UATT_core=[(0;r=0,Sigma),(Umin;r=Sigma,Rmin)]")  Type_attr=MFPAIR_SIGTOUMIN_UMIN;
    else if (Functional_List->get<string>("F2_PAIRPOTcore_Functional")=="UATTcore=[0;r=0,Rzero]")                         Type_attr=MFPAIR_RCSZERO_ZERO;

    if (Functional_List->get<string>("F3_CHARGE_Functional")=="No Charge or No Poisson")           Type_coul=NONE;
    else if (Functional_List->get<string>("F3_CHARGE_Functional")=="Charge_Mean_Field")            Type_coul=BARE;
    else if (Functional_List->get<string>("F3_CHARGE_Functional")=="Charge_with_DeltaC_RPM")       Type_coul=DELTAC_RPM;
    else if (Functional_List->get<string>("F3_CHARGE_Functional")=="Charge_with_DeltaC_General")   Type_coul=DELTAC_GENERAL;
    else if (Functional_List->get<string>("F3_CHARGE_Functional")=="Charge(MF)_with_Polarization") Type_coul=POLARIZE;

    if (Functional_List->get<string>("F4_POLYMER_Functional")=="No Bonded Fluids")      Type_poly=NONE;
    else if (Functional_List->get<string>("F4_POLYMER_Functional")=="Polymer_CMS")      Type_poly=CMS;
    else if (Functional_List->get<string>("F4_POLYMER_Functional")=="Polymer_CMS_SCFT") Type_poly=CMS_SCFT;
    else if (Functional_List->get<string>("F4_POLYMER_Functional")=="Polymer_TC_iSAFT")                    Type_poly=WTC;
    else if (Functional_List->get<string>("F4_POLYMER_Functional")=="Polymer_JDC_iSAFT(seg)")              Type_poly=WJDC;
    else if (Functional_List->get<string>("F4_POLYMER_Functional")=="Polymer_JDC_iSAFT(segRho compField)") Type_poly=WJDC2;
    else if (Functional_List->get<string>("F4_POLYMER_Functional")=="Polymer_JDC_iSAFT(comp)")             Type_poly=WJDC3;

    /****************************************************/
    /* params from surface geometry section of the GUI  */
    /****************************************************/

    Nwall=Surface_List->get<int>("S1: Number of Surfaces");
    Nlink=Surface_List->get<int>("S2: Number of macro surfaces");
    Nwall_type=Surface_List->get<int>("S3: Number of surface types");

    Array<string> WallTypes=SurfaceGeometry_List->get<Array<string> >("SG1: Surface Type Array");
    for (i=0; i<Nwall_type; i++){
       if (WallTypes[i]=="PLANE: Infinite in two dimensions")       Surface_type[i]=0;
       else if (WallTypes[i]=="PLANE: Infinite in one dimension")   Surface_type[i]=1;
       else if (WallTypes[i]=="CYLINDERICAL SURFACE: Infinite_L")   Surface_type[i]=2;
       else if (WallTypes[i]=="CYLINDRICAL PORE: Infinite_L")       Surface_type[i]=7;
       else if (WallTypes[i]=="SLIT PORE: Finite_L")                Surface_type[i]=8;
       else if (WallTypes[i]=="SLIT TAPERED PORE: Finite_L")        Surface_type[i]=9;
       else if (WallTypes[i]=="PLANE: Finite three dimensions")     Surface_type[i]=1;
       else if (WallTypes[i]=="SPHERICAL SURFACE: COLLOID")         Surface_type[i]=2;
       else if (WallTypes[i]=="SPHERICAL SURFACE: ATOM")            Surface_type[i]=3;
       else if (WallTypes[i]=="SPHERICAL SURFACE: POINT ATOM")      Surface_type[i]=4;
       else if (WallTypes[i]=="SPHERICAL CAVITY")                   Surface_type[i]=7;
       else if (WallTypes[i]=="CYLINDERICAL SURFACE: Finite_L")     Surface_type[i]=5;
       else if (WallTypes[i]=="CYLINDERICAL SURFACE PERIODIC: Infinite_L") Surface_type[i]=6;
       else if (WallTypes[i]=="CYLINDRICAL PORE: Finite_L")         Surface_type[i]=8;
       else if (WallTypes[i]=="CYLINDRICAL TAPERED PORE: Finite_L") Surface_type[i]=9;
   }
    
       
  return;
}

