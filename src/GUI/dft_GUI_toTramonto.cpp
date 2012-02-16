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
                        Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                        Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List,
                        Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                        Teuchos::RCP<Teuchos::ParameterList> PolymerGraft_List,
                        Teuchos::RCP<Teuchos::ParameterList> PolymerArch_List,
                        Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List,
                        Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                        Teuchos::RCP<Teuchos::ParameterList> SurfaceGeometry_List) 
{
  int i,j,k,counter;
  /****************************** DIMENSION PARAMETER SECTION **********************************************************/
  /* this routine translates the parameters from the GUI to Tramonto.  */


    /***************************************/
    /* params from mesh section of the GUI */
    /***************************************/
    
    if (Mesh_List->get<bool>("M2_Dimensionless_Distance_Entry")) Length_ref=-1.0; 
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
       else if (BC_Xmax[i]=="In_Wall") Type_bc[i][1]=IN_WALL;
       else if (BC_Xmax[i]=="Reflect") Type_bc[i][1]=REFLECT;
       else if (BC_Xmax[i]=="Periodic") Type_bc[i][1]=PERIODIC;
       else if (BC_Xmax[i]=="Last_node") Type_bc[i][1]=LAST_NODE;
       else if (BC_Xmax[i]=="Last_node_restart") Type_bc[i][1]=LAST_NODE_RESTART;
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
    /* params from fluid physics section of the GUI  */
    /****************************************************/
   Ncomp=Fluid_List->get<int>("F1_Ncomp");

    if(Fluid_List->get<string>("F2_HSDiamType")=="HS_diam=Sigma") Type_hsdiam=SIGMA_DIAM;
    else if(Fluid_List->get<string>("F2_HSDiamType")=="Manual Definition") Type_hsdiam=MANUAL_HS_DIAM;
    else if(Fluid_List->get<string>("F2_HSDiamType")=="Barker-Henderson") Type_hsdiam=BH_DIAM;
    
    Array<double> F3_HSDiam = Fluid_List->get<Array<double> >("F3_HSDiam");
    for (i=0; i<Ncomp; i++) HS_diam[i]=F3_HSDiam[i]; 

    if(Fluid_List->get<string>("F4_PairPotType")=="none") Type_pairPot=PAIR_HARD;
    else if(Fluid_List->get<string>("F4_PairPotType")=="LJ 12-6 potential (cut/shift)") Type_pairPot=PAIR_LJ12_6_CS;
    else if(Fluid_List->get<string>("F4_PairPotType")=="Coulomb potential as mean field (cut/shift)") Type_pairPot=PAIR_COULOMB_CS;
    else if(Fluid_List->get<string>("F4_PairPotType")=="Coulomb potential as mean field (cut only)") Type_pairPot=PAIR_COULOMB;
    else if(Fluid_List->get<string>("F4_PairPotType")=="Yukawa potential (cut/shift)") Type_pairPot=PAIR_YUKAWA_CS;
    else if(Fluid_List->get<string>("F4_PairPotType")=="Exponential potential (cut/shift)") Type_pairPot=PAIR_EXP_CS;
    else if(Fluid_List->get<string>("F4_PairPotType")=="Square Well potential") Type_pairPot=PAIR_SW;
    else if(Fluid_List->get<string>("F4_PairPotType")=="LJ 12-6 plus Yukawa potential (cut/shift)") Type_pairPot=PAIR_LJandYUKAWA_CS;
    else if(Fluid_List->get<string>("F4_PairPotType")=="r^12 repulsion plus Yukawa potential (cut/shift)") Type_pairPot=PAIR_r12andYUKAWA_CS;
    else if(Fluid_List->get<string>("F4_PairPotType")=="r^18 repulsion plus Yukawa potential (cut/shift)") Type_pairPot=PAIR_r18andYUKAWA_CS;
    else if(Fluid_List->get<string>("F4_PairPotType")=="r^N repulsion plus Yukawa potential (cut/shift)") Type_pairPot=PAIR_rNandYUKAWA_CS;

    Temp=Fluid_List->get<double>("F6_Temperature");
   
    if (PotentialsFF_List->get<string>("PF0_Off_Diagonal_Definitions")=="Manual Definition") Mix_type=1;
    else if (PotentialsFF_List->get<string>("PF0_Off_Diagonal_Definitions")=="Lorentz-Berthlot Mixing") Mix_type=0;

    if (Mix_type==0){ 
       Array<double> PF1_SigmaF=PotentialsFF_List->get<Array<double> >("PF1_SigmaF");
       for (i=0; i<Ncomp; i++) Sigma_ff[i][i]=PF1_SigmaF[i]; 

       Array<double> PF2_EpsF=PotentialsFF_List->get<Array<double> >("PF2_EpsF");
       for (i=0; i<Ncomp; i++) Eps_ff[i][i]=PF2_EpsF[i]; 

       Array<double> PF3_CutF=PotentialsFF_List->get<Array<double> >("PF3_CutF");
       for (i=0; i<Ncomp; i++) Cut_ff[i][i]=PF3_CutF[i]; 
    }
    else{
       TwoDArray<double> PF1_SigmaFF=PotentialsFF_List->get<TwoDArray<double> >("PF1_SigmaFF");
       for (i=0; i<Ncomp; i++)  
          for (j=0; j<Ncomp; j++){ Sigma_ff[i][j]=PF1_SigmaFF[i][j]; 
       }

       TwoDArray<double> PF2_EpsFF=PotentialsFF_List->get<TwoDArray<double> >("PF2_EpsFF");
       for (i=0; i<Ncomp; i++)  
          for (j=0; j<Ncomp; j++) Eps_ff[i][j]=PF2_EpsFF[i][j]; 

       TwoDArray<double> PF3_CutFF=PotentialsFF_List->get<TwoDArray<double> >("PF3_CutFF");
       for (i=0; i<Ncomp; i++)  
          for (j=0; j<Ncomp; j++) Cut_ff[i][j]=PF3_CutFF[i][j]; 

       TwoDArray<double> PF4_EpsYukawaFF=PotentialsFF_List->get<TwoDArray<double> >("PF4_EpsYukawaFF");
       for (i=0; i<Ncomp; i++)  
          for (j=0; j<Ncomp; j++) EpsYukawa_ff[i][j]=PF4_EpsYukawaFF[i][j]; 

       TwoDArray<double> PF5_YukawaKFF=PotentialsFF_List->get<TwoDArray<double> >("PF5_ExpDecayParamFF");
       for (i=0; i<Ncomp; i++)  
          for (j=0; j<Ncomp; j++) YukawaK_ff[i][j]=PF5_YukawaKFF[i][j]; 

       TwoDArray<double> PF6_NpowFF=PotentialsFF_List->get<TwoDArray<double> >("PF6_NpowFF");
       for (i=0; i<Ncomp; i++)  
          for (j=0; j<Ncomp; j++) Npow_ff[i][j]=PF6_NpowFF[i][j]; 

       TwoDArray<double> PF10_BondFF=PotentialsFF_List->get<TwoDArray<double> >("PF10_BondFF");
       for (i=0; i<Ncomp; i++)  
          for (j=0; j<Ncomp; j++) Bond_ff[i][j]=PF10_BondFF[i][j]; 
    }
 

    Array<double> PF7_Mass=PotentialsFF_List->get<Array<double> >("PF7_Mass");
    for (i=0; i<Ncomp; i++) Mass[i]=PF7_Mass[i]; 

    Array<double> PF8_Charge=PotentialsFF_List->get<Array<double> >("PF8_Charge");
    for (i=0; i<Ncomp; i++) Charge_f[i]=PF8_Charge[i]; 

    Array<double> PF9_Pol=PotentialsFF_List->get<Array<double> >("PF9_Polarization");
    for (i=0; i<Ncomp; i++) Pol[i]=PF9_Pol[i]; 

    /****************************************************/
    /* params from polymer input section of the GUI  */
    /****************************************************/
    Npol_comp=Polymer_List->get<int>("P1: Npoly_comp");

    Array<int> P2_Nblock=Polymer_List->get<Array<int> >("P2: Nblock_per_polymer");
    for (i=0; i<Npol_comp; i++) Nblock[i]=P2_Nblock[i]; 

     /** note that Nseg_per_block is the same as block[][] in dft_input.c...Need to use it to compute some other things,
        but we don't want to do those things here if possible ***/
    TwoDArray<int> P4_SegsPerBlock=Polymer_List->get<TwoDArray<int> >("P4: Nseg_perBlock");
      for (i=0; i<Npol_comp; i++)  
        for (j=0; j<Nblock[i]; j++) Nseg_per_block[i][j]=P4_SegsPerBlock[i][j]; 

     
    TwoDArray<int> P5_SegTypePerBlock=Polymer_List->get<TwoDArray<int> >("P5: SegType_perBlock");
      for (i=0; i<Npol_comp; i++)  
        for (j=0; j<Nblock[i]; j++) SegType_per_block[i][j]=P5_SegTypePerBlock[i][j]; 

    if (Polymer_List->get<string>("P6: Polymer achitecture entry")=="Read From File") Type_poly_arch=POLY_ARCH_FILE;
    else if (Polymer_List->get<string>("P6: Polymer achitecture entry")=="Linear Chains - Automatic set-up") Type_poly_arch=LIN_POLY;
    else if (Polymer_List->get<string>("P6: Polymer achitecture entry")=="Symmetric Linear Chains - Automatic set-up") Type_poly_arch=LIN_POLY_SYM;
    else Type_poly_arch=SET_IN_GUI;
   
    if (Polymer_List->get<string>("P6: Polymer achitecture entry")=="Read From File") {
       Poly_file_name=(char*)Polymer_List->get<string>("P7: Polymer architecture filename").c_str();
    }
    else Poly_file_name=NULL;


           /* Variables specific to Grafted polymers */
    if (PolymerGraft_List->get<bool>("PG1: Grafted Polymers?")){
/**       Array<bool> PG2_GraftPolymerTF=PolymerGraft_List->get<Array<bool> >("PG2: Grafted_polymer_TF");
       for (i=0; i<Npol_comp; i++) if(PG2_GraftPolymerTF[i]) Grafted[i]=TRUE;
       else Grafted[i]=FALSE; ***/

       Array<int> PG3_GraftWallID=PolymerGraft_List->get<Array<int> >("PG3: Grafted_wall_ID[ipol_comp]");
       for (i=0; i<Npol_comp; i++) {
            Graft_wall[i]=PG3_GraftWallID[i]; 
            if (Graft_wall[i]<0) Grafted[i]=FALSE;
            else Grafted[i]=TRUE;
       }

       Array<double> PG4_GraftDensity=PolymerGraft_List->get<Array<double> >("PG4: Grafted_wall_Density[ipol_comp]");
       for (i=0; i<Npol_comp; i++) Rho_g[i]=PG4_GraftDensity[i]; 
    }
    
          /* CMS specific variables */
   if (Type_poly==CMS || Type_poly==CMS_SCFT){
       Ncr_files=PolymerCMS_List->get<int>("CMS1: N_CrFiles");
       Cr_file=(char*)PolymerCMS_List->get<string>("CMS2: Cr_File_1").c_str();
       Cr_file2=(char*)PolymerCMS_List->get<string>("CMS3: Cr_File_2").c_str();
       Crfac=PolymerCMS_List->get<double>("CMS4: CrFac");

       TwoDArray<double> CMS5_CradHSArray=PolymerCMS_List->get<TwoDArray<double> >("CMS5: Cr HSRadius");
       for (i=0; i<Ncomp; i++)  for (j=0; j<Ncomp; j++) Cr_rad_hs[i][j]=CMS5_CradHSArray[i][j]; 
    }
    

    if (Polymer_List->get<string>("P6: Polymer achitecture entry")=="Set up in GUI") {

       Nseg_tot=PolymerArch_List->get<int>("PA1: NSeg_tot");
       Nbond_max=PolymerArch_List->get<int>("PA2: Nbond_max");
        
       Array<int> PA3_NbondPerSeg_Array=PolymerArch_List->get<Array<int> >("PA3: NBondsperSeg"); 
       counter=0;
       for (i=0;i<Npol_comp;i++)
          for (j=0;j<Nmer[i];j++){ 
             Nbond[i][j]=PA3_NbondPerSeg_Array[counter++]; 
          }

       TwoDArray<int> PA4_BondsPerSeg=PolymerArch_List->get<TwoDArray<int> >("PA4: BondAll");
       counter=0;
       for (i=0;i<Npol_comp;i++)
          for (j=0;j<Nmer[i];j++){
              for (k=0;k<Nbond[i][j];k++) { Bonds[i][j][k]=PA4_BondsPerSeg[counter][k]; }
               counter++;
           }

       TwoDArray<int> PA5_PolSymArray=PolymerArch_List->get<TwoDArray<int> >("PA5: PolySym");
       counter=0;
       for (i=0;i<Npol_comp;i++)
          for (j=0;j<Nmer[i];j++){
              for (k=0;k<Nbond[i][j];k++) { pol_sym_tmp[i][j][k]=PA5_PolSymArray[counter][k];}
              counter++;
          }

    }


    /****************************************************/
    /* params from surface geometry section of the GUI  */
    /****************************************************/

 /*   Nwall=Surface_List->get<int>("S1: Number of Surfaces");
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
*/

    /****************************************************/
    /* params from thermodynamics section of the GUI  */
    /****************************************************/

    /****************************************************/
    /* params from solvers section of the GUI  */
    /****************************************************/

    /****************************************************/
    /* params from continuation section of the GUI  */
    /****************************************************/
    
       
  return;
}

