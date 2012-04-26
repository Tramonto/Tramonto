using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_potentialsWW_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List)
{

  /****************************************************************************************************************/
  /********************************* SURFACE DEFINITIONS SECTION **************************************************/
  /****************************************************************************************************************/

     /* VALIDATORS*/

   RCP<StringValidator> SurfElecBCType_Validator = rcp(
      new StringValidator(tuple<std::string>( "none - neutral surfaces",
          "constant potential surfaces", "constant surface charge density", "discrete atomic charges")));

   RCP<ArrayStringValidator> SurfElecBCArray_Val=rcp(new ArrayStringValidator(SurfElecBCType_Validator));

   RCP<StringValidator> ChargeType_Validator = rcp(new StringValidator(tuple<std::string>(
	   "point charge", "charged smeared in spherical volume")));

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
           "r^18 repulsion plus Yukawa potential (cut/shift)")));


   /* SET PARAMEBERS */

   PotentialsWW_List->set(" Any Spherical Surfaces?",false,"Set to true if there are any spherical surfaces in this system." );
   PotentialsWW_List->set(" Compute interactions for spherical surfaces?",false,"Indicate whether interactions between spherical surfaces to be output." );
   PotentialsWW_List->set(" Type of wall-wall interactions","none","Identify type of interactions between fixed spherical surfaces. Selection applies to all.",PairPotValidator);

       /* Using mixing rules for cross terms */
   Array<double> SigmaW_Array(Surface_List->get<int>("S3: Number of surface types"),1.0);
   PotentialsWW_List->set("W3: Sigma_w",SigmaW_Array,"Array for characteristic length in wall-wall interactions.");

   Array<double> EpsW_Array(Surface_List->get<int>("S3: Number of surface types"),1.0);
   PotentialsWW_List->set("W4: Eps_w",EpsW_Array,"Array for energy parameters in wall-wall interactions.");

   Array<double> CutW_Array(Surface_List->get<int>("S3: Number of surface types"),3.0);
   PotentialsWW_List->set("W5: Cut_w",CutW_Array,"Array for cutoff distances in wall-wall interactions.");

   Array<double> EpsYukawaW_Array(Surface_List->get<int>("S3: Number of surface types"),1.0);
   PotentialsWW_List->set("W6: EpsYukawa_w",EpsYukawaW_Array,"Array for prefactor of Yukawa term in wall-wall interactions.");

   Array<double> ExpDecayParamW_Array(Surface_List->get<int>("S3: Number of surface types"),1.0);
   PotentialsWW_List->set("W7: ExpDecayParam_w",ExpDecayParamW_Array,"Array for Exponential Decay Parameter in wall-wall interactions.");

       /* Manual entry of all cross terms */
   TwoDArray<double> SigmaWW_Array(Surface_List->get<int>("S3: Number of surface types")*Surface_List->get<int>("S3: Number of surface types"),1.0);
   PotentialsWW_List->set("WW3: Sigma_ww",SigmaWW_Array,"Array for characteristic length in wall-wall interactions.");

   TwoDArray<double> EpsWW_Array(Surface_List->get<int>("S3: Number of surface types")*Surface_List->get<int>("S3: Number of surface types"),1.0);
   PotentialsWW_List->set("WW4: Eps_ww",EpsWW_Array,"Array for energy parameters in wall-wall interactions.");

   TwoDArray<double> CutWW_Array(Surface_List->get<int>("S3: Number of surface types")*Surface_List->get<int>("S3: Number of surface types"),3.0);
   PotentialsWW_List->set("WW5: Cut_ww",CutWW_Array,"Array for cutoff distances in wall-wall interactions.");

   TwoDArray<double> EpsYukawaWW_Array(Surface_List->get<int>("S3: Number of surface types")*Surface_List->get<int>("S3: Number of surface types"),1.0);
   PotentialsWW_List->set("WW6: EpsYukawa_ww",EpsYukawaWW_Array,"Array for prefactor of Yukawa term in wall-wall interactions.");

   TwoDArray<double> ExpDecayParamWW_Array(Surface_List->get<int>("S3: Number of surface types")*Surface_List->get<int>("S3: Number of surface types"),1.0);
   PotentialsWW_List->set("WW7: ExpDecayParam_ww",ExpDecayParamWW_Array,"Array for Exponential Decay Parameter in wall-wall interactions.");

   /* now set up some charged surface parameters */

   SurfaceParamCharge_List->set("SC0: Are any Surfaces Charged?",false,"Set to true if any of the defined surfaces in the problem carry a surface charge.");

   Array<double> DielecW_Array(Surface_List->get<int>("S3: Number of surface types"),1.0);
   SurfaceParamCharge_List->set("SC0.1: DielecConst_Wall Array",DielecW_Array,"Set the dielectric constant in each of the surfaces.  Entries will be reduced by reference value if one was provided.");

   Array<string> TypeBCElec_Array(Surface_List->get<int>("S3: Number of surface types"),"none - neutral surfaces");
   SurfaceParamCharge_List->set("SC1: Type_elec_BC",TypeBCElec_Array,"Select type of boundary condition you would like to apply at the charged surfaces.",SurfElecBCArray_Val);

   SurfaceParamCharge_List->set("SC1.1: Treatment of atomic charges","point charge","Enter method for treating fixed atomic charges.");

   SurfaceParamCharge_List->set("SC2: Number of additional Charge Sources",0,"Enter number of localized fixed source charges in the system that are not associated with a defined surfaces.");
   SurfaceParamCharge_List->set("SC2.1: Treatment of charge sources","point charge","Enter method for treating charge sources at fixed positions.");

   Array<double> ChargeLoc_Array(SurfaceParamCharge_List->get<int>("SC2: Number of additional Charge Sources"),0.0);
   SurfaceParamCharge_List->set("SC2.2: Charge of the source points",ChargeLoc_Array,"Enter the value of the charge that is present at a source point.");

   Array<double> ChargeDiam_Array(SurfaceParamCharge_List->get<int>("SC2: Number of additional Charge Sources"),0.0);
   SurfaceParamCharge_List->set("SC2.3: Charge Source Diameter",ChargeDiam_Array,"Enter the diameter over which the locallized source charges should be spread.");

   TwoDArray<double> ChargePos_Array(Mesh_List->get<int>("M1_Ndim")*SurfaceParamCharge_List->get<int>("SC2: Number of additional Charge Sources"),0.0);
   SurfaceParamCharge_List->set("SC2.4: Charge Source Position",ChargePos_Array,"Enter the location of each charge source.");

   Array<double> SurfaceCharge_Array(Surface_List->get<int>("S1: Number of Surfaces"),0.0);
   SurfaceParamCharge_List->set("SC3: Surface Charge or Potential", SurfaceCharge_Array, "Charge (or surface potential) on each surface in the system.");


   return;
}
/*****************************************************************************************************************************************/
void dft_GUI_potentialsWW_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List)
{
   bool tmp_bool=false;
   string tmp_string;
   int i,j;
   double tmp_1D[NWALL_MAX_TYPE];
   string tmpstring_1D[NWALL_MAX_TYPE];

   RCP<StringValidator> SurfElecBCType_Validator = rcp(
      new StringValidator(tuple<std::string>( "none - neutral surfaces",
           "constant potential surfaces", "constant surface charge density", "discrete atomic charges")));
   RCP<ArrayStringValidator> SurfElecBCArray_Val=rcp(new ArrayStringValidator(SurfElecBCType_Validator));

   RCP<StringValidator> ChargeType_Validator = rcp(new StringValidator(tuple<std::string>(
	   "point charge", "charged smeared in spherical volume")));

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
           "r^18 repulsion plus Yukawa potential (cut/shift)")));

   for (i=0;i<Nwall_type;i++) if (Surface_type[i]==colloids_cyl_sphere || Surface_type[i]==atomic_centers || Surface_type[i]==point_surface) tmp_bool=true;
   PotentialsWW_List->set(" Any Spherical Surfaces?",tmp_bool,"Set to true if there are any spherical surfaces in this system." );

   tmp_bool=false;
   if (Lprint_pmf) tmp_bool=true;
   PotentialsWW_List->set(" Compute interactions for spherical surfaces?",tmp_bool,"Indicate whether interactions between spherical surfaces to be output." );

   if (Type_uwwPot==PAIR_HARD) tmp_string="none";
   else if (Type_uwwPot==PAIR_LJ12_6_CS) tmp_string="LJ 12-6 potential (cut/shift)";
   else if (Type_uwwPot==PAIR_COULOMB_CS) tmp_string="Coulomb potential as mean field (cut/shift)";
   else if (Type_uwwPot==PAIR_COULOMB) tmp_string="Coulomb potential as mean field (cut only)";
   else if (Type_uwwPot==PAIR_YUKAWA_CS) tmp_string="Yukawa potential (cut/shift)";
   else if (Type_uwwPot==PAIR_EXP_CS) tmp_string="Exponential potential (cut/shift)";
   else if (Type_uwwPot==PAIR_SW) tmp_string="Square Well potential";
   else if (Type_uwwPot==PAIR_LJandYUKAWA_CS) tmp_string="LJ 12-6 plus Yukawa potential (cut/shift)";
   else if (Type_uwwPot==PAIR_r12andYUKAWA_CS) tmp_string="r^12 repulsion plus Yukawa potential (cut/shift)";
   else  if (Type_uwwPot==PAIR_r18andYUKAWA_CS) tmp_string="r^18 repulsion plus Yukawa potential (cut/shift)";
   else tmp_string="none";
  
   PotentialsWW_List->set(" Type of wall-wall interactions",tmp_string,"Identify type of interactions between fixed spherical surfaces. Selection applies to all.",PairPotValidator);

       /* Using mixing rules for cross terms */
   for (i=0;i<Nwall_type;i++) tmp_1D[i]=Sigma_ww[i][i];
   Array<double> SigmaW_Array(tmp_1D, tmp_1D+Nwall_type);
   PotentialsWW_List->set("W3: Sigma_w",SigmaW_Array,"Array for characteristic length in wall-wall interactions.");

   for (i=0;i<Nwall_type;i++) tmp_1D[i]=Eps_ww[i][i];
   Array<double> EpsW_Array(tmp_1D, tmp_1D+Nwall_type);
   PotentialsWW_List->set("W4: Eps_w",EpsW_Array,"Array for energy parameters in wall-wall interactions.");

   for (i=0;i<Nwall_type;i++) tmp_1D[i]=Cut_ww[i][i];
   Array<double> CutW_Array(tmp_1D, tmp_1D+Nwall_type);
   PotentialsWW_List->set("W5: Cut_w",CutW_Array,"Array for cutoff distances in wall-wall interactions.");

   for (i=0;i<Nwall_type;i++) tmp_1D[i]=EpsYukawa_ww[i][i];
   Array<double> EpsYukawaW_Array(tmp_1D, tmp_1D+Nwall_type);
   PotentialsWW_List->set("W6: EpsYukawa_w",EpsYukawaW_Array,"Array for prefactor of Yukawa term in wall-wall interactions.");

   for (i=0;i<Nwall_type;i++) tmp_1D[i]=YukawaK_ww[i][i];
   Array<double> ExpDecayParamW_Array(tmp_1D, tmp_1D+Nwall_type);
   PotentialsWW_List->set("W7: ExpDecayParam_w",ExpDecayParamW_Array,"Array for Exponential Decay Parameter in wall-wall interactions.");


       /* Manually set cross terms */
   TwoDArray<double> SigmaWW_Array(Nwall_type,Nwall_type);
   for (i=0; i<Nwall_type;i++) for (j=0; j<Nwall_type;j++) SigmaWW_Array(i,j)=Sigma_ww[i][j];
   PotentialsWW_List->set("WW3: Sigma_ww",SigmaWW_Array,"Array for characteristic length in wall-wall interactions.");

   TwoDArray<double> EpsWW_Array(Nwall_type,Nwall_type);
   for (i=0; i<Nwall_type;i++) for (j=0; j<Nwall_type;j++) EpsWW_Array(i,j)=Eps_ww[i][j];
   PotentialsWW_List->set("WW4: Eps_ww",EpsWW_Array,"Array for energy parameters in wall-wall interactions.");

   TwoDArray<double> CutWW_Array(Nwall_type,Nwall_type);
   for (i=0; i<Nwall_type;i++) for (j=0; j<Nwall_type;j++) CutWW_Array(i,j)=Cut_ww[i][j];
   PotentialsWW_List->set("WW5: Cut_ww",CutWW_Array,"Array for cutoff distances in wall-wall interactions.");
TwoDArray<double> EpsYukawaWW_Array(Nwall_type,Nwall_type);
   for (i=0; i<Nwall_type;i++) for (j=0; j<Nwall_type;j++) EpsYukawaWW_Array(i,j)=EpsYukawa_ww[i][j];
   PotentialsWW_List->set("WW6: EpsYukawa_ww",EpsYukawaWW_Array,"Array for prefactor of Yukawa term in wall-wall interactions.");

   TwoDArray<double> ExpDecayParamWW_Array(Nwall_type,Nwall_type);
   for (i=0; i<Nwall_type;i++) for (j=0; j<Nwall_type;j++) ExpDecayParamWW_Array(i,j)=YukawaK_ww[i][j];
   PotentialsWW_List->set("WW7: ExpDecayParam_ww",ExpDecayParamWW_Array,"Array for Exponential Decay Parameter in wall-wall interactions.");

   /* now set up some charged surface parameters */

   tmp_bool=false;
   for (i=0; i<Nwall_type;i++){
      if (Type_bc_elec[i]==CONST_CHARGE){ tmpstring_1D[i]="constant surface charge density"; tmp_bool=true;}
      else if (Type_bc_elec[i]==CONST_POTENTIAL){ tmpstring_1D[i]="constant potential surfaces"; tmp_bool=true;}
      else if (Type_bc_elec[i]==ATOMIC_CHARGE) {tmpstring_1D[i]="discrete atomic charges"; tmp_bool=true;}
      else tmpstring_1D[i]="none - neutral surfaces";
   }
   SurfaceParamCharge_List->set("SC0: Are any Surfaces Charged?",tmp_bool,"Set to true if any of the defined surfaces in the problem carry a surface charge.");


   Array<string> TypeBCElec_Array(tmpstring_1D, tmpstring_1D+Nwall_type);
   SurfaceParamCharge_List->set("SC1: Type_elec_BC",TypeBCElec_Array,"Select type of boundary condition you would like to apply at the charged surfaces.",SurfElecBCArray_Val);

   if (Type_coul != NONE){
      for (i=0;i<Nwall_type;i++) tmp_1D[i]=Dielec_wall[i];
      Array<double> DielecW_Array(tmp_1D, tmp_1D+Nwall_type);
      SurfaceParamCharge_List->set("SC0.1: DielecConst_Wall Array",DielecW_Array,"Set the dielectric constant in each of the surfaces.  Entries will be reduced by reference value if one was provided.");
   }

   if (Charge_type_atoms==POINT_CHARGE) SurfaceParamCharge_List->set("SC1.1: Treatment of atomic charges","point charge","Enter method for treating fixed atomic charges.");
   else     SurfaceParamCharge_List->set("SC1.1: Treatment of atomic charges","charged smeared in spherical volume","Enter method for treating fixed atomic charges.");



   SurfaceParamCharge_List->set("SC2: Number of additional Charge Sources",Nlocal_charge,"Enter number of localized fixed source charges in the system that are not associated with a defined surfaces.");

   if (Charge_type_local==POINT_CHARGE) SurfaceParamCharge_List->set("SC2.1: Treatment of charge sources","point charge","Enter method for treating charge sources at fixed positions.");
   else     SurfaceParamCharge_List->set("SC2.1: Treatment of charge sources","charged smeared in spherical volume","Enter method for treating charge sources at fixed positions.");

   for (i=0;i<Nlocal_charge;i++) tmp_1D[i]=Charge[i];
   Array<double> ChargeLoc_Array(tmp_1D,tmp_1D+Nlocal_charge);
   SurfaceParamCharge_List->set("SC2.2: Charge of the source points",ChargeLoc_Array,"Enter the value of the charge that is present at a source point.");

   for (i=0;i<Nlocal_charge;i++) tmp_1D[i]=Charge_Diam[i];
   Array<double> ChargeDiam_Array(tmp_1D,tmp_1D+Nlocal_charge);
   SurfaceParamCharge_List->set("SC2.3: Charge Source Diameter",ChargeDiam_Array,"Enter the diameter over which the locallized source charges should be spread.");

   TwoDArray<double> ChargePos_Array(Ndim,Nlocal_charge);
   for (i=0; i<Ndim;i++) for (j=0; j<Nlocal_charge;j++) ChargePos_Array(i,j)=Charge_x[i][j];
   SurfaceParamCharge_List->set("SC2.4: Charge Source Position",ChargePos_Array,"Enter the location of each charge source.");

   for (i=0;i<Nwall;i++) tmp_1D[i]=Elec_param_w[i];
   Array<double> SurfaceCharge_Array(tmp_1D,tmp_1D+Nwall);
   SurfaceParamCharge_List->set("SC3: Surface Charge or Potential", SurfaceCharge_Array, "Charge (or surface potential) on each surface in the system.");

   return;
}
/*****************************************************************************************************************************************/
void dft_GUI_potentialsWW_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List,
                      Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List)
{
     /* DEPENDENCIES */
  RCP<StringCondition> ManualMixCon = rcp(
     new StringCondition(PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"), "Manual Definition"));

  RCP<StringCondition> LBMixCond = rcp(
     new StringCondition(PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"), "Lorentz-Berthlot Mixing"));

  RCP<BoolCondition> SphericalSurfBoolCon = rcp(new BoolCondition(PotentialsWW_List->getEntryRCP(" Any Spherical Surfaces?")));
  RCP<BoolCondition> SphericalSurfBoolCon2 = rcp(new BoolCondition(PotentialsWW_List->getEntryRCP(" Compute interactions for spherical surfaces?")));

  Condition::ConstConditionList Spherical_conList=tuple<RCP<const Condition> >(SphericalSurfBoolCon,SphericalSurfBoolCon2);
  RCP<AndCondition> SphericalSurf_request = rcp(new AndCondition(Spherical_conList));

  Condition::ConstConditionList WWcompute_conList=tuple<RCP<const Condition> >(LBMixCond,SphericalSurfBoolCon2);
  RCP<OrCondition> WWcompute_ORcondition = rcp(new OrCondition(WWcompute_conList));
  RCP<ConditionVisualDependency> WWcompute_Dep = rcp( new ConditionVisualDependency(WWcompute_ORcondition, PotentialsWW_List->getEntryRCP(" Type of wall-wall interactions"), true));

      /* set up SigmaWW conditions */
  Condition::ConstConditionList SigmaWW_conList=tuple<RCP<const Condition> >(SphericalSurfBoolCon,SphericalSurfBoolCon2,ManualMixCon);
  RCP<AndCondition> SigmaWW_allowed = rcp(new AndCondition(SigmaWW_conList));
  RCP<ConditionVisualDependency> SigmaWW_Dep = rcp( new ConditionVisualDependency(SigmaWW_allowed, PotentialsWW_List->getEntryRCP("WW3: Sigma_ww"), true));

      /* now set up SigmaW conditions */
  Condition::ConstConditionList SigmaW_conList=tuple<RCP<const Condition> >(LBMixCond);
  RCP<AndCondition> SigmaW_allowed = rcp(new AndCondition(SigmaW_conList));
  RCP<ConditionVisualDependency> SigmaW_Dep = rcp(new ConditionVisualDependency(SigmaW_allowed, PotentialsWW_List->getEntryRCP("W3: Sigma_w"), true));

  RCP<StringCondition> EpsStringCon1 = rcp(
     new StringCondition(PotentialsWW_List->getEntryRCP(" Type of wall-wall interactions"),
                         tuple<std::string>("LJ 12-6 potential (cut/shift)",
                            "Exponential potential (cut/shift)","Square Well potential",
                            "LJ 12-6 plus Yukawa potential (cut/shift)",
                            "r^12 repulsion plus Yukawa potential (cut/shift)",
                            "r^18 repulsion plus Yukawa potential (cut/shift)")));

      /* set up EpsWW conditions */
  Condition::ConstConditionList EpsWW_conList=tuple<RCP<const Condition> >(SphericalSurfBoolCon,SphericalSurfBoolCon2,ManualMixCon,EpsStringCon1);
  RCP<AndCondition> EpsWW_allowed = rcp(new AndCondition(EpsWW_conList));
  RCP<ConditionVisualDependency> EpsWW_Dep = rcp( new ConditionVisualDependency(EpsWW_allowed, PotentialsWW_List->getEntryRCP("WW4: Eps_ww"), true));

      /* now set up EpsW conditions */
  Condition::ConstConditionList EpsW_conList=tuple<RCP<const Condition> >(EpsStringCon1,LBMixCond);
  RCP<AndCondition> EpsW_allowed = rcp(new AndCondition(EpsW_conList));
  RCP<ConditionVisualDependency> EpsW_Dep = rcp(new ConditionVisualDependency(EpsW_allowed, PotentialsWW_List->getEntryRCP("W4: Eps_w"), true));

  RCP<StringCondition> CutStringCon1 = rcp(
      new StringCondition(PotentialsWW_List->getEntryRCP(" Type of wall-wall interactions"),
                                 tuple<std::string>("LJ 12-6 potential (cut/shift)","Exponential potential (cut/shift)",
                                 "Coulomb potential as mean field (cut/shift)","Coulomb potential as mean field (cut only)",
                                 "Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
                                 "r^12 repulsion plus Yukawa potential (cut/shift)",
                                 "r^18 repulsion plus Yukawa potential (cut/shift)")));

      /* set up CutWW conditions */
  Condition::ConstConditionList CutWW_conList=tuple<RCP<const Condition> >(SphericalSurfBoolCon,SphericalSurfBoolCon2,ManualMixCon,CutStringCon1);
  RCP<AndCondition> CutWW_allowed = rcp(new AndCondition(CutWW_conList));
  RCP<ConditionVisualDependency> CutWW_Dep = rcp( new ConditionVisualDependency(CutWW_allowed, PotentialsWW_List->getEntryRCP("WW5: Cut_ww"), true));

      /* now set up CutW conditions */
  Condition::ConstConditionList CutW_conList=tuple<RCP<const Condition> >(CutStringCon1,LBMixCond);
  RCP<AndCondition> CutW_allowed = rcp(new AndCondition(CutW_conList));
  RCP<ConditionVisualDependency> CutW_Dep = rcp(new ConditionVisualDependency(CutW_allowed, PotentialsWW_List->getEntryRCP("W5: Cut_w"), true));

  RCP<StringCondition> EpsYukawaStringCon1 = rcp(
      new StringCondition(PotentialsWW_List->getEntryRCP(" Type of wall-wall interactions"),
                                 tuple<std::string>("Yukawa potential (cut/shift)",
                                 "LJ 12-6 plus Yukawa potential (cut/shift)",
                                 "r^12 repulsion plus Yukawa potential (cut/shift)",
                                 "r^18 repulsion plus Yukawa potential (cut/shift)")));

      /* set up EpsYukawaWW conditions */
  Condition::ConstConditionList EpsYukawaWW_conList=tuple<RCP<const Condition> >(SphericalSurfBoolCon,SphericalSurfBoolCon2,ManualMixCon,EpsYukawaStringCon1);
  RCP<AndCondition> EpsYukawaWW_allowed = rcp(new AndCondition(EpsYukawaWW_conList));
  RCP<ConditionVisualDependency> EpsYukawaWW_Dep = rcp( new ConditionVisualDependency(EpsYukawaWW_allowed, PotentialsWW_List->getEntryRCP("WW6: EpsYukawa_ww"), true));

      /* now set up EpsYukawaW conditions */
  Condition::ConstConditionList EpsYukawaW_conList=tuple<RCP<const Condition> >(EpsYukawaStringCon1,LBMixCond);
  RCP<AndCondition> EpsYukawaW_allowed = rcp(new AndCondition(EpsYukawaW_conList));
  RCP<ConditionVisualDependency> EpsYukawaW_Dep = rcp(new ConditionVisualDependency(EpsYukawaW_allowed, PotentialsWW_List->getEntryRCP("W6: EpsYukawa_w"), true));

  RCP<StringCondition> ExpDecayParamStringCon1 = rcp(
      new StringCondition(PotentialsWW_List->getEntryRCP(" Type of wall-wall interactions"),
                                 tuple<std::string>("Exponential potential (cut/shift)",
                                 "Yukawa potential (cut/shift)","LJ 12-6 plus Yukawa potential (cut/shift)",
                                 "r^12 repulsion plus Yukawa potential (cut/shift)",
                                 "r^18 repulsion plus Yukawa potential (cut/shift)")));

      /* set up ExpDecayParamWW conditions */
  Condition::ConstConditionList ExpDecayParamWW_conList=tuple<RCP<const Condition> >(SphericalSurfBoolCon,SphericalSurfBoolCon2,ManualMixCon,ExpDecayParamStringCon1);
  RCP<AndCondition> ExpDecayParamWW_allowed = rcp(new AndCondition(ExpDecayParamWW_conList));
  RCP<ConditionVisualDependency> ExpDecayParamWW_Dep = rcp( new ConditionVisualDependency(ExpDecayParamWW_allowed, PotentialsWW_List->getEntryRCP("WW7: ExpDecayParam_ww"), true));

      /* now set up ExpDecayParamW conditions */
  Condition::ConstConditionList ExpDecayParamW_conList=tuple<RCP<const Condition> >(ExpDecayParamStringCon1,LBMixCond);
  RCP<AndCondition> ExpDecayParamW_allowed = rcp(new AndCondition(ExpDecayParamW_conList));
  RCP<ConditionVisualDependency> ExpDecayParamW_Dep = rcp(new ConditionVisualDependency(ExpDecayParamW_allowed, PotentialsWW_List->getEntryRCP("W7: ExpDecayParam_w"), true));


   Dependency::ParameterEntryList PotWW2DRowNumber_Dependents;
   PotWW2DRowNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW3: Sigma_ww"));
   PotWW2DRowNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW4: Eps_ww"));
   PotWW2DRowNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW5: Cut_ww"));
   PotWW2DRowNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW6: EpsYukawa_ww"));
   PotWW2DRowNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW7: ExpDecayParam_ww"));
   RCP<TwoDRowDependency<int,double> > PotWW2DRowNumber_Dep = rcp(
             new TwoDRowDependency<int,double>(Surface_List->getEntryRCP("S3: Number of surface types"),PotWW2DRowNumber_Dependents));

   Dependency::ParameterEntryList PotWW2DColNumber_Dependents;
   PotWW2DColNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW3: Sigma_ww"));
   PotWW2DColNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW4: Eps_ww"));
   PotWW2DColNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW5: Cut_ww"));
   PotWW2DColNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW6: EpsYukawa_ww"));
   PotWW2DColNumber_Dependents.insert(PotentialsWW_List->getEntryRCP("WW7: ExpDecayParam_ww"));
   RCP<TwoDColDependency<int,double> > PotWW2DColNumber_Dep = rcp(
             new TwoDColDependency<int,double>(Surface_List->getEntryRCP("S3: Number of surface types"),PotWW2DColNumber_Dependents));

   Dependency::ParameterEntryList PotWWArrayLength_Dependents;
   PotWWArrayLength_Dependents.insert(PotentialsWW_List->getEntryRCP("W3: Sigma_w"));
   PotWWArrayLength_Dependents.insert(PotentialsWW_List->getEntryRCP("W4: Eps_w"));
   PotWWArrayLength_Dependents.insert(PotentialsWW_List->getEntryRCP("W5: Cut_w"));
   PotWWArrayLength_Dependents.insert(PotentialsWW_List->getEntryRCP("W6: EpsYukawa_w"));
   PotWWArrayLength_Dependents.insert(PotentialsWW_List->getEntryRCP("W7: ExpDecayParam_w"));
   PotWWArrayLength_Dependents.insert(SurfaceParamCharge_List->getEntryRCP("SC0.1: DielecConst_Wall Array"));
   RCP<NumberArrayLengthDependency<int,double> > PotWWArrayLength_Dep = rcp(
             new NumberArrayLengthDependency<int,double>(Surface_List->getEntryRCP("S3: Number of surface types"), PotWWArrayLength_Dependents));

   RCP<NumberArrayLengthDependency<int,string> > TypeElecArrayLength_Dep = rcp(
             new NumberArrayLengthDependency<int,string>(Surface_List->getEntryRCP("S3: Number of surface types"), SurfaceParamCharge_List->getEntryRCP("SC1: Type_elec_BC")));

   RCP<StringVisualDependency> DielecW_Dep = rcp(
        new StringVisualDependency(ChargedFluid_List->getEntryRCP("CF1: Type Dielectric Constant(s)"),
                                   SurfaceParamCharge_List->getEntryRCP("SC0.1: DielecConst_Wall Array"),
                                    tuple<std::string>("2-Dielec Const: fluid and wall regions",
                                      "3-Dielec Const: bulk fluid, fluid near wall, wall regions"), true));

   Dependency::ParameterEntryList SourceChargesLength_Deps;
   SourceChargesLength_Deps.insert(SurfaceParamCharge_List->getEntryRCP("SC2.2: Charge of the source points"));
   SourceChargesLength_Deps.insert(SurfaceParamCharge_List->getEntryRCP("SC2.3: Charge Source Diameter"));
   RCP<NumberArrayLengthDependency<int,double> > SourceChargeLength_Dep = rcp(
             new NumberArrayLengthDependency<int,double>(SurfaceParamCharge_List->getEntryRCP("SC2: Number of additional Charge Sources"), SourceChargesLength_Deps));

   Dependency::ParameterEntryList SurfaceCharge_Deps;
   SurfaceCharge_Deps.insert(SurfaceParamCharge_List->getEntryRCP("SC1: Type_elec_BC"));
   SurfaceCharge_Deps.insert(SurfaceParamCharge_List->getEntryRCP("SC1.1: Treatment of atomic charges"));
   RCP<BoolVisualDependency> SurfaceCharge_Dep = rcp(
        new BoolVisualDependency(SurfaceParamCharge_List->getEntryRCP("SC0: Are any Surfaces Charged?"), SurfaceCharge_Deps, true));

   Dependency::ParameterEntryList SourceCharge_Deps;
   SourceCharge_Deps.insert(SurfaceParamCharge_List->getEntryRCP("SC2.1: Treatment of charge sources"));
   SourceCharge_Deps.insert(SurfaceParamCharge_List->getEntryRCP("SC2.2: Charge of the source points"));
   SourceCharge_Deps.insert(SurfaceParamCharge_List->getEntryRCP("SC2.4: Charge Source Position"));

   RCP<NumberVisualDependency<int> > SourceCharge_Dep = rcp(
        new NumberVisualDependency<int>( SurfaceParamCharge_List->getEntryRCP("SC2: Number of additional Charge Sources"), SourceCharge_Deps));

   RCP<StringVisualDependency> ChargeLocDiam_Dep = rcp(
        new StringVisualDependency(SurfaceParamCharge_List->getEntryRCP("SC2.1: Treatment of charge sources"),
                                   SurfaceParamCharge_List->getEntryRCP("SC2.3: Charge Source Diameter"),
                                    tuple<std::string>("charged smeared in spherical volume"), true));

   RCP<TwoDColDependency<int,double> > SourceChargeColNumber_Dep = rcp(
             new TwoDColDependency<int,double>(SurfaceParamCharge_List->getEntryRCP("SC2: Number of additional Charge Sources"),
                                               SurfaceParamCharge_List->getEntryRCP("SC2.4: Charge Source Position")));

   RCP<TwoDRowDependency<int,double> > SourceChargeRowNumber_Dep = rcp(
             new TwoDRowDependency<int,double>(Mesh_List->getEntryRCP("M1_Ndim"),
                                               SurfaceParamCharge_List->getEntryRCP("SC2.4: Charge Source Position")));

   RCP<StringCondition> SurfCharge_StCon = rcp(
           new StringCondition(SurfacePosition_List->getEntryRCP("SP0: Type of surface position/charge entry"),"Set up in GUI"));

   RCP<BoolCondition> ChargedSurf_TFCond=rcp(new BoolCondition(SurfaceParamCharge_List->getEntryRCP("SC0: Are any Surfaces Charged?")));

   Condition::ConstConditionList SurfChargeNwall_ConList=tuple<RCP<const Condition> >(ChargedSurf_TFCond,SurfCharge_StCon);
   RCP<AndCondition>SurfChargeNwall_andCon=rcp(new AndCondition(SurfChargeNwall_ConList));
   
   RCP<ConditionVisualDependency> SurfChargeNwall_Dep=rcp(new ConditionVisualDependency(SurfChargeNwall_andCon, 
                                  SurfaceParamCharge_List->getEntryRCP("SC3: Surface Charge or Potential"),true));

   RCP<NumberArrayLengthDependency<int,double> > SurfChargeNwallLength_Dep = rcp(
             new NumberArrayLengthDependency<int,double>(Surface_List->getEntryRCP("S1: Number of Surfaces"), SurfaceParamCharge_List->getEntryRCP("SC3: Surface Charge or Potential")));

   /* DEPENDENCY SHEET ENTRIES*/
  depSheet_Tramonto->addDependency(WWcompute_Dep);
  depSheet_Tramonto->addDependency(SigmaWW_Dep);
  depSheet_Tramonto->addDependency(SigmaW_Dep);
  depSheet_Tramonto->addDependency(EpsWW_Dep);
  depSheet_Tramonto->addDependency(EpsW_Dep);
  depSheet_Tramonto->addDependency(CutWW_Dep);
  depSheet_Tramonto->addDependency(CutW_Dep);
  depSheet_Tramonto->addDependency(EpsYukawaWW_Dep);
  depSheet_Tramonto->addDependency(EpsYukawaW_Dep);
  depSheet_Tramonto->addDependency(ExpDecayParamWW_Dep);
  depSheet_Tramonto->addDependency(ExpDecayParamW_Dep);
  depSheet_Tramonto->addDependency(PotWW2DRowNumber_Dep);
  depSheet_Tramonto->addDependency(PotWW2DColNumber_Dep);
  depSheet_Tramonto->addDependency(PotWWArrayLength_Dep);
  depSheet_Tramonto->addDependency(SourceCharge_Dep);
  depSheet_Tramonto->addDependency(SurfaceCharge_Dep);
  depSheet_Tramonto->addDependency(ChargeLocDiam_Dep);
  depSheet_Tramonto->addDependency(SourceChargeLength_Dep);
  depSheet_Tramonto->addDependency(SourceChargeColNumber_Dep);
  depSheet_Tramonto->addDependency(SourceChargeRowNumber_Dep);
  depSheet_Tramonto->addDependency(DielecW_Dep);
  depSheet_Tramonto->addDependency(TypeElecArrayLength_Dep);
  depSheet_Tramonto->addDependency(SurfChargeNwall_Dep);
  depSheet_Tramonto->addDependency(SurfChargeNwallLength_Dep);
   return;
}

