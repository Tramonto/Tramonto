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
                        Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                        Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                        Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List,
                        Teuchos::RCP<Teuchos::ParameterList> Continuation_List,
                        Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                        Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                        Teuchos::RCP<Teuchos::ParameterList> LoadBalance_List,
                        Teuchos::RCP<Teuchos::ParameterList> PhysicsMethod_List,
                        Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List,
                        Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                        Teuchos::RCP<Teuchos::ParameterList> Output_List,
                        Teuchos::RCP<Teuchos::ParameterList> DensProfile_List,
                        Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                        Teuchos::RCP<Teuchos::ParameterList> SurfaceGeometry_List) 
{
  string str_tmp;
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
    /* params from state point section of the GUI  */
    /****************************************************/
    if (StatePoint_List->get<bool>("BF1: Dimensionless Density Entry?")) Density_ref=-1.0; 
    else    Density_ref=StatePoint_List->get<double>("BF2: Density Conversion Factor");

   
    if (Type_interface == UNIFORM_INTERFACE){
        if (Type_poly==NONE){ 
           Array<double> BF3_RhoLBB=StatePoint_List->get<Array<double> >("BF3: Rho_b_0[icomp]"); 
           for (i=0;i<Ncomp;i++){
              Rho_b[i]=BF3_RhoLBB[i]; 
              if (Type_interface==DIFFUSIVE_INTERFACE) D_coef[i]=0.0; 
           }
        }
        else{ 
            Array<double> BF3_RhoLBB=StatePoint_List->get<Array<double> >("BF3: Rho_b_0[ipol_comp]"); 
            for (i=0;i<Npol_comp;i++){
               Rho_b[i]=BF3_RhoLBB[i]; 
               if (Type_interface==DIFFUSIVE_INTERFACE) D_coef[i]=0.0;
            }
        }
    
        Grad_dim=0;
        Velocity=0.0;
        Lconstrain_interface=FALSE;
    }
    else{
       
        if (Type_poly==NONE){ 
           Array<double> BF3_RhoLBB=StatePoint_List->get<Array<double> >("BF3: Rho_b_0[icomp]"); 
           Array<double> BF4_RhoRTF=StatePoint_List->get<Array<double> >("BF4: Rho_b_1[icomp]"); 
           Array<double> D1_DiffCoeff=Diffusion_List->get<Array<double> >("D1: Diff_Coeff[icomp]"); 
           for (i=0;i<Ncomp;i++){ 
             Rho_b_LBB[i]=BF3_RhoLBB[i]; 
             Rho_b_RTF[i]=BF4_RhoRTF[i]; 
             Rho_b[i]=Rho_b_LBB[i];
             if (Type_interface==DIFFUSIVE_INTERFACE) D_coef[i]=D1_DiffCoeff[i];
           }
        }
        else{ 
           Array<double> BF3_RhoLBB=StatePoint_List->get<Array<double> >("BF3: Rho_b_0[ipol_comp]"); 
           Array<double> BF4_RhoRTF=StatePoint_List->get<Array<double> >("BF4: Rho_b_1[ipol_comp]"); 
           Array<double> D1_DiffCoeff=Diffusion_List->get<Array<double> >("D1: Diff_Coeff[ipol_comp]"); 
           for (i=0;i<Npol_comp;i++){ 
             Rho_b_LBB[i]=BF3_RhoLBB[i]; 
             Rho_b_RTF[i]=BF4_RhoRTF[i]; 
             Rho_b[i]=Rho_b_LBB[i];
             if (Type_interface==DIFFUSIVE_INTERFACE) D_coef[i]=D1_DiffCoeff[i];
           }
        }

        Grad_dim=StatePoint_List->get<int>("BF5: Direction of Gradient");
        if (Type_interface==DIFFUSIVE_INTERFACE) Velocity= Diffusion_List->get<double>("D3: Velocity");
        if (StatePoint_List->get<bool>("BF6: Constrained Interface?")) Lconstrain_interface=TRUE;
        else Lconstrain_interface=FALSE;
    }

    X_const_mu=StatePoint_List->get<double>("BF7: X_const");


    if (Type_coul != NONE){
       if (ChargedFluid_List->get<string>("CF1: Type Dielectric Constant(s)")=="Uniform Dielectric Constant") Type_dielec=DIELEC_CONST;
       else if (ChargedFluid_List->get<string>("CF1: Type Dielectric Constant(s)")=="2-Dielec Const: fluid and wall regions") Type_dielec=DIELEC_WF;
       else if (ChargedFluid_List->get<string>("CF1: Type Dielectric Constant(s)")=="3-Dielec Const: bulk fluid, fluid near wall, wall regions") Type_dielec=DIELEC_WF_PORE;

       if (ChargedFluid_List->get<bool>("CF2: Entry of Relative Dielectric Constant(s)?")) Dielec_ref=-1.;
       else Dielec_ref=ChargedFluid_List->get<double>("CF3: Reference Dielectric Constant");

       Dielec_bulk=ChargedFluid_List->get<double>("CF4.0: Dielec Const Bulk Fluid");
       Dielec_pore=ChargedFluid_List->get<double>("CF4.1: Dielec Const Near Wall Fluid");
       Dielec_X=ChargedFluid_List->get<double>("CF4.2: Size of Near Wall region");
  
       Sigma_Angstroms_plasma=ChargedFluid_List->get<double>("CF5.0: Sigma for Plasma Parameter (Angstroms)");
       Temp_K_plasma=ChargedFluid_List->get<double>("CF5.1: Temperature for Plasma Parameter (Kelvin)");
       DielecConst_plasma=ChargedFluid_List->get<double>("CF5.2: Dielectric const for Plasma Parameter");
 
       if (ChargedFluid_List->get<string>("CF6.0: Electrostatic Potential(s) entry type")=="Enter electrostatic potentials in mV units") Flag_mV_elecpot=TRUE;
       else Flag_mV_elecpot=FALSE;


       Elec_pot_LBB=ChargedFluid_List->get<double>("CF6.1: Elec_pot_0");
       Elec_pot_LBB=ChargedFluid_List->get<double>("CF6.2: Elec_pot_1");
    }

    /****************************************************/
    /* params from the continuation section of the GUI  */
    /****************************************************/
    Lbinodal=FALSE;
    if (Continuation_List->get<string>("C1: Continuation Type")=="None") Loca.method=-1;
    else if (Continuation_List->get<string>("C1: Continuation Type")=="LOCA: Simple Parameter Continuation") Loca.method=1;
    else if (Continuation_List->get<string>("C1: Continuation Type")=="LOCA: Arc-Length Parameter Continuation") Loca.method=2;
    else if (Continuation_List->get<string>("C1: Continuation Type")=="LOCA: Spinodal Continuation") Loca.method=3;
    else if (Continuation_List->get<string>("C1: Continuation Type")=="LOCA: Binodal Continuation"){
         Loca.method=4;
         Lbinodal=TRUE;
    }
    else if (Continuation_List->get<string>("C1: Continuation Type")=="Mesh Stepping with LOCA Binodal" ||
             Continuation_List->get<string>("C1: Continuation Type")=="Mesh Continuation"){

         Nsteps=Continuation_List->get<int>("C5: Number of Steps");
         Plane_new_nodes=Continuation_List->get<int>("C9.1: Direction of Mesh Cont.");
         for (i=0;i<Ndim;i++){
             if (i==Plane_new_nodes) Del_1[i]=Continuation_List->get<double>("C4: Continuation Step Size");
             else Del_1[i]=0.0;
         }
         if (Continuation_List->get<string>("C9.2: Location for node count change")=="center") Pos_new_nodes=0;
         else if (Continuation_List->get<string>("C9.2: Location for node count change")=="Left(dim=0), Bottom(dim=1), or Back(dim=2) boundary") Pos_new_nodes=-1;
         else if (Continuation_List->get<string>("C9.2: Location for node count change")=="Right(dim=0), Top(dim=1), or Front(dim=2) boundary") Pos_new_nodes=1;

         if (Continuation_List->get<string>("C1: Continuation Type")=="Mesh Stepping with LOCA Binodal"){
            Loca.method=4;
            Loca.cont_type1=CONT_MESH;
            Lbinodal=TRUE;
            Loca.step_size=Del_1[Plane_new_nodes];
            Loca.num_steps=0;
            Loca.aggr=0.0;
         }
         else if (Continuation_List->get<string>("C1: Continuation Type")=="Mesh Continuation") Loca.method=-1;          

    }

    if (Continuation_List->get<string>("C1: Continuation Type")!="None"){

        if (Continuation_List->get<string>("C2: Continuation Parameter")=="Temperature" ){
           Loca.cont_type1=CONT_TEMP;
        }
        else if (Continuation_List->get<string>("C2: Continuation Parameter")=="Density (1 species)" ||
                 Continuation_List->get<string>("C2: Continuation Parameter")=="Chemical Potential (1 species)"){
           if (Continuation_List->get<string>("C2: Continuation Parameter")=="Density (1 species)" ) Loca.cont_type1=CONT_RHO_I;
           else Loca.cont_type1=CONT_BETAMU_I;
           NID_Cont=1;
           Cont_ID[0][0]=Continuation_List->get<int>("C3: icomp");
        }
        else if (Continuation_List->get<string>("C2: Continuation Parameter")=="Electrostatic parameter (iwall)"){
           Loca.cont_type1=CONT_ELECPARAM_I;
           NID_Cont=1;
           Cont_ID[0][0]=Continuation_List->get<int>("C3: iwall");
        }
        else if (Continuation_List->get<string>("C2: Continuation Parameter")=="Wall-Wall energy param (iwall_type or iwall_type,jwall_type pair)"){
           Loca.cont_type1=CONT_EPSW_I;

           if (PotentialsFF_List->get<string>("PF0_Off_Diagonal_Definitions")=="Lorentz-Berthlot Mixing"){
              NID_Cont=1;
              Cont_ID[0][0]=Continuation_List->get<int>("C3: iwall_type");
           }
           else{
              NID_Cont=2;
              Array<int> C3_WW = Continuation_List->get<Array<int> >("C3: iwall_type,jwall_type");
              for (i=0; i<2; i++) Cont_ID[0][i]=C3_WW[i]; 
           }
        }
        else if (Continuation_List->get<string>("C2: Continuation Parameter")=="Fluid-fluid energy param (i or ij pair)"){
           Loca.cont_type1=CONT_EPSFF_IJ;
           if (PotentialsFF_List->get<string>("PF0_Off_Diagonal_Definitions")=="Lorentz-Berthlot Mixing"){
              NID_Cont=1;
              Cont_ID[0][0]=Continuation_List->get<int>("C3: icomp");
           }
           else{
              NID_Cont=2;
              Array<int> C3_FF = Continuation_List->get<Array<int> >("C3: icomp,jcomp");
              for (i=0; i<2; i++) Cont_ID[0][i]=C3_FF[i]; 
           }
        }
        else if (Continuation_List->get<string>("C2: Continuation Parameter")=="Fluid size: SigmaFF(i or ij pair)"){
           Loca.cont_type1=CONT_SIGMAFF_IJ;
           if (PotentialsFF_List->get<string>("PF0_Off_Diagonal_Definitions")=="Lorentz-Berthlot Mixing"){
              NID_Cont=1;
              Cont_ID[0][0]=Continuation_List->get<int>("C3: icomp");
           }
           else{
              NID_Cont=2;
              Array<int> C3Sig_FF = Continuation_List->get<Array<int> >("C3: icomp,jcomp");
              for (i=0; i<2; i++) Cont_ID[0][i]=C3Sig_FF[i]; 
           }
        }
        else if (Continuation_List->get<string>("C2: Continuation Parameter")=="Wall-Fluid energy param (iwall_type,icomp)"){
              Loca.cont_type1=CONT_EPSWF_IJ;
              NID_Cont=2;
              Array<int> C3_WF = Continuation_List->get<Array<int> >("C3: iwall_type,icomp");
              for (i=0; i<2; i++) Cont_ID[0][i]=C3_WF[i]; 
        }

        Loca.step_size=Continuation_List->get<double>("C4: Continuation Step Size");
        Loca.num_steps=Continuation_List->get<int>("C5: Number of Steps");
        Loca.aggr=Continuation_List->get<double>("C6: Step Aggressivnes");

    }

    if (Continuation_List->get<string>("C1: Continuation Type")=="LOCA: Binodal Continuation" || 
        Continuation_List->get<string>("C1: Continuation Type")=="Mesh Stepping with LOCA Binodal"){

        if (Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Temperature" ){
           Loca.cont_type2=CONT_TEMP;
        }
        else if (Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Density (1 species)" ||
                 Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Chemical Potential (1 species)"){
           if (Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Density (1 species)" ) {
                 Loca.cont_type2=CONT_RHO_I;
           }
           else  Loca.cont_type2=CONT_BETAMU_I;
           NID_Cont=1;
           Cont_ID[1][0]=Continuation_List->get<int>("C8: icomp");
        }
        else if (Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Electrostatic parameter (iwall)"){
           Loca.cont_type2=CONT_ELECPARAM_I;
           NID_Cont=1;
           Cont_ID[1][0]=Continuation_List->get<int>("C8: iwall");
        }
        else if (Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Wall-Wall energy param (iwall_type or iwall_type,jwall_type pair)"){
           Loca.cont_type2=CONT_EPSW_I;

           if (PotentialsFF_List->get<string>("PF0_Off_Diagonal_Definitions")=="Lorentz-Berthlot Mixing"){
              NID_Cont=1;
              Cont_ID[1][0]=Continuation_List->get<int>("C8: iwall_type");
           }
           else{
              NID_Cont=2;
              Array<int> C8_WW = Continuation_List->get<Array<int> >("C8: iwall_type,jwall_type");
              for (i=0; i<2; i++) Cont_ID[1][i]=C8_WW[i]; 
           }
        }
        else if (Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Fluid-fluid energy param (i or ij pair)"){
           Loca.cont_type2=CONT_EPSFF_IJ;
           if (PotentialsFF_List->get<string>("PF0_Off_Diagonal_Definitions")=="Lorentz-Berthlot Mixing"){
              NID_Cont=1;
              Cont_ID[1][0]=Continuation_List->get<int>("C8: icomp");
           }
           else{
              NID_Cont=2;
              Array<int> C8_FF = Continuation_List->get<Array<int> >("C8: icomp,jcomp");
              for (i=0; i<2; i++) Cont_ID[1][i]=C8_FF[i]; 
           }
        }
        else if (Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Fluid size: SigmaFF(i or ij pair)"){
           Loca.cont_type2=CONT_SIGMAFF_IJ;
           if (PotentialsFF_List->get<string>("PF0_Off_Diagonal_Definitions")=="Lorentz-Berthlot Mixing"){
              NID_Cont=1;
              Cont_ID[1][0]=Continuation_List->get<int>("C8: icomp");
           }
           else{
              NID_Cont=2;
              Array<int> C8Sig_FF = Continuation_List->get<Array<int> >("C8: icomp,jcomp");
              for (i=0; i<2; i++) Cont_ID[1][i]=C8Sig_FF[i]; 
           }
        }
        else if (Continuation_List->get<string>("C7: Dependent Cont. Param. (Binodals)")=="Wall-Fluid energy param (iwall_type,icomp)"){
              Loca.cont_type2=CONT_EPSWF_IJ;
              NID_Cont=2;
              Array<int> C8_WF = Continuation_List->get<Array<int> >("C8: iwall_type,icomp");
              for (i=0; i<2; i++) Cont_ID[1][i]=C8_WF[i]; 
        }
     }

    /*****************************************************/
    /* params from numerical methods section of the GUI  */
    /*****************************************************/
      if (LoadBalance_List->get<string>("LB1: Load Balancing Approach")=="Linear Matrix Balance") Load_Bal_Flag=LB_LINEAR;
      else if (LoadBalance_List->get<string>("LB1: Load Balancing Approach")=="Weighted Recursive Bisection") Load_Bal_Flag=LB_WEIGHTS;
      else if (LoadBalance_List->get<string>("LB1: Load Balancing Approach")=="Geometric Recursive Bisection") Load_Bal_Flag=LB_BOX;

     if (Coarsening_List->get<string>("C1: Coarsening Type")=="none") Mesh_coarsening=FALSE;
     else if (Coarsening_List->get<string>("C1: Coarsening Type")=="residual and jacobian coarsening") Mesh_coarsening=TRUE;
     else if (Coarsening_List->get<string>("C1: Coarsening Type")=="Poisson-Boltzmann zone") Mesh_coarsening=PB_ZONE;
     else if (Coarsening_List->get<string>("C1: Coarsening Type")=="Bulk zone") Mesh_coarsening=BULK_ZONE;

     if (Coarsening_List->get<bool>("C5.0: Truncate jacobian integrals?")){
        Lcut_jac=TRUE;
        Jac_threshold=Coarsening_List->get<double>("C5.1: Truncation threshhold");
     }
     else{
        Lcut_jac=FALSE;
        Jac_threshold=0.0;
     }

     if(Coarsening_List->get<string>("C2.0: Type of Jacobian coarsening")=="Jacobian Coarsening identical to Residual Coarsening") Coarser_jac=JAC_RESID_ZONES_SAME;
     else if(Coarsening_List->get<string>("C2.0: Type of Jacobian coarsening")=="Jacobian: factor of 2 more coarse than resid in most refined zone") Coarser_jac=JAC_ZONE0_FAC2LESSTHANRESID;
     else if(Coarsening_List->get<string>("C2.0: Type of Jacobian coarsening")=="Jacobian: factor of 2 more coarse than resid in all but most coarse zone") Coarser_jac=JAC_ZONES_FAC2LESSTHANRESID;
     else if(Coarsening_List->get<string>("C2.0: Type of Jacobian coarsening")=="Jacobian: use most coarse zone to define entire matrix") Coarser_jac=JAC_ZONES_ALLMOSTCOARSE;
     else if(Coarsening_List->get<string>("C2.0: Type of Jacobian coarsening")=="Jacobian: use 2nd most coarse zone in all but most coarse region") Coarser_jac=JAC_ZONES_SECONDMOSTCOARSE;
     else if(Coarsening_List->get<string>("C2.0: Type of Jacobian coarsening")=="Jacobian: set Esize_jac for all matrix calculations") Coarser_jac=JAC_ZONES_SETFIXED_ESIZE;

     Jac_grid=Coarsening_List->get<double>("C2.1: Esize_Jac");

     Nzone = Coarsening_List->get<int>("C3: Nzone");

     Array<double> Rmin_Array = Coarsening_List->get<Array<double> >("C4: Rmin for each zone");
     for (i=0;i<Nzone-1;i++) Rmax_zone[i]=Rmin_Array[i+1];

     if (Coarsening_List->get<bool>("C6.0: 1D boundary zone?")){
       L1D_bc=TRUE;
       Dim_1Dbc=Coarsening_List->get<int>("C6.1: Dim_1D_bc");
       X_1D_bc=Coarsening_List->get<double>("C6.2: X_1D_bc");
     }
     else L1D_bc=FALSE;

     if(PhysicsMethod_List->get<bool>("PM1: Attractions in A22 Block of matrix?")) ATTInA22Block=TRUE;
     else ATTInA22Block=FALSE;

     if(PhysicsMethod_List->get<string>("PM2: Physics Scaling?")=="No Physics Scaling") Physics_scaling=FALSE;
     else if(PhysicsMethod_List->get<string>("PM2: Physics Scaling?")=="Automatic Calculation") Physics_scaling=AUTOMATIC;
     else if(PhysicsMethod_List->get<string>("PM2: Physics Scaling?")=="Manual Input") Physics_scaling=MANUAL_INPUT;

     if(PhysicsMethod_List->get<bool>("PM3: Analytic Jacobian?")) Analyt_WJDC_Jac=TRUE;
     else Analyt_WJDC_Jac=FALSE;

     Array<double> PM4_ScaleFac_Array=PhysicsMethod_List->get<Array<double> >("PM4: ScaleFac[icomp]");
     for (i=0; i<Npol_comp;i++) 
          for (j=0; j<Nblock[i];j++){ 
             Scale_fac_WJDC[i][SegType_per_block[i][j]]=PM4_ScaleFac_Array[SegType_per_block[i][j]]; 
     }

     if(NonlinearSolver_List->get<string>("NLS1: Nonlinear Solver")=="Newton Built-In") NL_Solver=NEWTON_BUILT_IN;
     else if(NonlinearSolver_List->get<string>("NLS1: Nonlinear Solver")=="Newton NOX") NL_Solver=NEWTON_NOX;
     if(NonlinearSolver_List->get<string>("NLS1: Nonlinear Solver")=="Picard Built-In") NL_Solver=PICARD_BUILT_IN;
     else if(NonlinearSolver_List->get<string>("NLS1: Nonlinear Solver")=="Picard NOX") NL_Solver=PICARD_NOX;
     else if(NonlinearSolver_List->get<string>("NLS1: Nonlinear Solver")=="Picard/Newton Built-In") NL_Solver=PICNEWTON_BUILT_IN;
     else if(NonlinearSolver_List->get<string>("NLS1: Nonlinear Solver")=="Picard/Newton NOX") NL_Solver=PICNEWTON_NOX;

     Max_NL_iter=NonlinearSolver_List->get<int>("NLS2: Max Nonlinear Iterations");
     NL_abs_tol=NonlinearSolver_List->get<double>("NLS3: Newton Tolerance Absolute");
     NL_rel_tol=NonlinearSolver_List->get<double>("NLS4: Newton Tolerance Relative");
     NL_update_scalingParam=NonlinearSolver_List->get<double>("NLS5: Minimum update fraction");
     NL_abs_tol_picard=NonlinearSolver_List->get<double>("NLS6: Picard Tolerance Absolute");
     NL_rel_tol_picard=NonlinearSolver_List->get<double>("NLS7: Picard Tolerance Relative");

     if(LinearSolver_List->get<string>("LS4: Linear Solver Approach")=="GMRES") Az_solver=0;
     else if(LinearSolver_List->get<string>("LS4: Linear Solver Approach")=="cg") Az_solver=1;
     else if(LinearSolver_List->get<string>("LS4: Linear Solver Approach")=="tfqmr") Az_solver=2;
     else if(LinearSolver_List->get<string>("LS4: Linear Solver Approach")=="cgs") Az_solver=3;
     else if(LinearSolver_List->get<string>("LS4: Linear Solver Approach")=="bicgstab") Az_solver=4;

     if(LinearSolver_List->get<string>("LS5: Matrix Scaling option")=="none") Az_scaling=-1;
     else if(LinearSolver_List->get<string>("LS5: Matrix Scaling option")=="row_sum") Az_scaling=0;
     else if(LinearSolver_List->get<string>("LS5: Matrix Scaling option")=="jacobi") Az_scaling=1;
     else if(LinearSolver_List->get<string>("LS5: Matrix Scaling option")=="symrow_sum") Az_scaling=2;

     if(LinearSolver_List->get<string>("LS6: Preconditioner option")=="none") Az_preconditioner=-1;
     else if(LinearSolver_List->get<string>("LS6: Preconditioner option")=="ilu") Az_preconditioner=0;
     else if(LinearSolver_List->get<string>("LS6: Preconditioner option")=="jacobi") Az_preconditioner=1;
     else if(LinearSolver_List->get<string>("LS6: Preconditioner option")=="symmetric Gauss-Seidel") Az_preconditioner=2;
     else if(LinearSolver_List->get<string>("LS6: Preconditioner option")=="LSpoly3") Az_preconditioner=3;
     else if(LinearSolver_List->get<string>("LS6: Preconditioner option")=="ilut") Az_preconditioner=4;

     Az_ilut_fill_param=LinearSolver_List->get<double>("LS7: Number of Fill Levels for ILUT");

    /****************************************************/
    /* params from output selection section of the GUI  */
    /****************************************************/
     if(Output_List->get<string>("O1: Screen Output")=="Basic Screen Output") Iwrite_screen=SCREEN_BASIC;
     else if (Output_List->get<string>("O1: Screen Output")=="No Screen Output") Iwrite_screen=SCREEN_NONE;
     else if (Output_List->get<string>("O1: Screen Output")=="Errors Only to Screen") Iwrite_screen=SCREEN_ERRORS_ONLY;
     else if (Output_List->get<string>("O1: Screen Output")=="Verbose Screen Output") Iwrite_screen=SCREEN_VERBOSE;
     else if (Output_List->get<string>("O1: Screen Output")=="Debugging Output: Residual") Iwrite_screen=SCREEN_DEBUG_RESID;

     if(Output_List->get<string>("O2: Files Output")=="Basic Files") Iwrite_files=FILES_BASIC;
     else if(Output_List->get<string>("O2: Files Output")=="Extended Files") Iwrite_files=FILES_EXTENDED;
     else if(Output_List->get<string>("O2: Files Output")=="Debug Files") Iwrite_files=FILES_DEBUG;
     else if(Output_List->get<string>("O2: Files Output")=="Debug Matrix Files") Iwrite_files=FILES_DEBUG_MATRIX;

     Iwrite=DENSITIES;
     if (Iwrite_screen==SCREEN_VERBOSE || Iwrite_files==FILES_DEBUG) Iwrite=VERBOSE;
     if (Iwrite_files==FILES_EXTENDED) Iwrite=EXTENDED;
     else if (Iwrite_files==FILES_DEBUG_MATRIX) Iwrite=VERBOSE_MATRIX;
     if (Iwrite_screen==SCREEN_NONE) Iwrite=NO_SCREEN;

     if(Output_List->get<string>("O3.0: dft_output.dat: State Point Output")=="Density[i]") Print_rho_switch=SWITCH_RHO;
     else if(Output_List->get<string>("O3.0: dft_output.dat: State Point Output")=="Betamu[i]") Print_rho_switch=SWITCH_MU;
     else if(Output_List->get<string>("O3.0: dft_output.dat: State Point Output")=="kappa") Print_rho_switch=SWITCH_ION;
     else if(Output_List->get<string>("O3.0: dft_output.dat: State Point Output")=="Density[i], Betamu[i], kappa, and pressure") Print_rho_switch=SWITCH_ALLTYPES_ICOMP;
     else if(Output_List->get<string>("O3.0: dft_output.dat: State Point Output")=="Density[i], Betamu[i], and pressure") Print_rho_switch=SWITCH_ALLTYPES_ICOMP;
     else if(Output_List->get<string>("O3.0: dft_output.dat: State Point Output")=="Density[all], Betamu[all], kappa, and pressure") Print_rho_switch=SWITCH_ALLTYPES;
     else if(Output_List->get<string>("O3.0: dft_output.dat: State Point Output")=="Density[all], Betamu[all], and pressure") Print_rho_switch=SWITCH_ALLTYPES;
     else if(Output_List->get<string>("O3.0: dft_output.dat: State Point Output")=="No State Point Output") Print_rho_switch=SWITCH_NO_STATEOUT;

     if(Output_List->get<string>("O3.1: dft_output.dat: Adsorption & Energy Output")=="adsorption/volume and energy/volume (bulk density & pressure)"){
         if(Print_rho_switch==SWITCH_ALLTYPES)  Print_rho_switch=SWITCH_BULK_OUTPUT_ALL;
         else Print_rho_switch=SWITCH_BULK_OUTPUT;
     }

     if (Output_List->get<bool>("O3.2: dft_output.dat: per unit area?")) Lper_area=TRUE;
     else Lper_area=FALSE;

     if (Output_List->get<bool>("O3.3: dft_output.dat: Correct for reflections?")) Lcount_reflect=TRUE;
     else Lcount_reflect=FALSE;

     if (Output_List->get<string>("O3.4: dft_output.dat: Mesh Output")=="Position of all surfaces")  Print_mesh_switch=TRUE;
     else if (Output_List->get<string>("O3.4: dft_output.dat: Mesh Output")=="Separations between surfaces")  Print_mesh_switch=SWITCH_SURFACE_SEP;
   
     if (Output_List->get<bool>("O4: Print g(r)?")) Lprint_gofr=TRUE;
     else Lprint_gofr=FALSE;

     if (Output_List->get<bool>("O5: Print surface-surface interactions?")) Lprint_pmf=TRUE;
     else Lprint_pmf=FALSE;

    /***********************************************************/
    /* params from initial guess selection section of the GUI  */
    /***********************************************************/
     str_tmp=DensProfile_List->get<string>("DP1.0: Type of Initial Density Profile");
     if (str_tmp=="Construct a new Density Profile") Restart=NORESTART;
     else if (str_tmp=="Restart from File") Restart=RESTART_BASIC;
     else if (str_tmp=="Restart with Step to constant") Restart=RESTART_STEP;
     else if (str_tmp=="Restart densities only (no other fields)") Restart=RESTART_DENSONLY;
     else if (str_tmp=="Restart incomplete Ncomp") Restart=RESTART_FEWERCOMP;
     else if (str_tmp=="Restart with 1D profile (in 2D or 3D)") Restart=RESTART_1DTOND;

    if (str_tmp !="Construct a new Density Profile") {
       DensityFile=(char*)DensProfile_List->get<string>("DP1.1: Density Restart File").c_str();
       if (Continuation_List->get<string>("C1: Continuation Type")=="LOCA: Binodal Continuation" ||
           Continuation_List->get<string>("C1: Continuation Type")=="Mesh Stepping with LOCA Binodal" ||
           Continuation_List->get<string>("C1: Continuation Type")=="LOCA: Spinodal Continuation" ){
                DensityFile2=(char*)DensProfile_List->get<string>("DP1.2: 2nd Density RestartFile").c_str();
        }
        else DensityFile2=NULL;
    }
    else{
       DensityFile=NULL;
       DensityFile2=NULL;
    }
   
    Nmissing_densities=DensProfile_List->get<int>("DP1.3: Number of missing components");
    Rho_max=DensProfile_List->get<double>("DP1.4: Rho max");

    str_tmp=DensProfile_List->get<string>("DP2.0: Density Profile Construct Type");
    if (str_tmp=="Constant Bulk Density") Iguess=CONST_RHO;
    else if (str_tmp=="Rho_bulk*exp(-Vext/kT)") Iguess=EXP_RHO;
    else if (str_tmp=="Step function profile") Iguess=STEP_PROFILE;
    else if (str_tmp=="Chopped profile to rho_bulk") Iguess=CHOP_RHO;
    else if (str_tmp=="Chopped profile to rho_step") Iguess=CHOP_RHO_STEP;
    else if (str_tmp=="Linear profile (for diffusion)") Iguess=LINEAR;

    Nsteps=DensProfile_List->get<int>("DP2.1: Nsteps");

    Array<int> OrientStep_Array = DensProfile_List->get<Array<int> >("DP2.2: Orient_step[istep]");
    Array<double> XstartStep_Array = DensProfile_List->get<Array<double> >("DP2.3: Xstart_step[istep]");
    Array<double> XendStep_Array = DensProfile_List->get<Array<double> >("DP2.4: Xend_step[istep]");
    for (i=0;i<Nsteps;i++){
        Orientation_step[i]=OrientStep_Array[i];
        Xstart_step[i]=XstartStep_Array[i];
        Xend_step[i]=XendStep_Array[i];
    }
    TwoDArray<double> Rhomax_Array=DensProfile_List->get<TwoDArray<double> >("DP2.5 Rho_step[icomp][istep]");
    for (i=0; i<Ncomp; i++)  for (j=0; j<Nsteps; j++) Rho_step[i][j]=Rhomax_Array[i][j]; 

    str_tmp=DensProfile_List->get<string>("DP3: Dependent Field Construct Type");
    if (str_tmp=="Bulk values for dependent fields") Iguess_fields=BULK;
    else if (str_tmp=="Compute fields based on density profile") Iguess_fields=CALC_ALL_FIELDS;
    else if (str_tmp=="Compute nonlocal densities only - other fields are bulk") Iguess_fields=CALC_RHOBAR_ONLY;
    else if (str_tmp=="Compute chain eq. and nonlocal densities - other fields are bulk") Iguess_fields=CALC_RHOBAR_AND_G;

    str_tmp=DensProfile_List->get<string>("DP4.0: External Field Restart Type");
    if (str_tmp=="Compute new Vext") Restart_Vext=READ_VEXT_FALSE;
    else if (str_tmp=="Restart from File") Restart_Vext=READ_VEXT_TRUE;
    else if (str_tmp=="Sum Two from Files") Restart_Vext=READ_VEXT_SUMTWO;
    else if (str_tmp=="Sum Two with constraints") Restart_Vext=READ_VEXT_STATIC;

    if (str_tmp !="Compute new Vext") {
       Vext_filename=(char*)DensProfile_List->get<string>("DP4.1: External Field Filename").c_str();
       if (str_tmp=="Sum Two from Files" || str_tmp=="Sum Two with constraints") 
             Vext_filename2=(char*)DensProfile_List->get<string>("DP4.2: 2nd External Field Filename").c_str();
       else Vext_filename2=NULL;
    }
    else{
       Vext_filename=NULL;
       Vext_filename2=NULL;
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
    /* params from solvers section of the GUI  */
    /****************************************************/

    /****************************************************/
    /* params from continuation section of the GUI  */
    /****************************************************/
    
       
  return;
}

