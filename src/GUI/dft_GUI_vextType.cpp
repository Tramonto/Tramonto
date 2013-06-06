//using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_vextType_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List)
{
     /* VALIDATORS*/

        RCP<StringValidator> MethodType_Valid = rcp(new StringValidator(tuple<std::string>(
                       "none",
                       "Pure exclusions / hard surface",
                       "V(r or x=dist. to surface)",
                       "V(r or x=dist. to center)",
                       "V=integral(U(r))")));

        RCP<ArrayStringValidator> MethodType_ValidArray = rcp(new ArrayStringValidator(MethodType_Valid));

        RCP<StringValidator> VextType_Valid = rcp(new StringValidator(tuple<std::string>(
                       "Vext(x)", "U(r) for Vext(U(r))")));
        RCP<ArrayStringValidator> VextType_ValidArray = rcp(new ArrayStringValidator(VextType_Valid));

        RCP<StringValidator> U3D_Valid = rcp(new StringValidator(tuple<std::string>(
                       "none",
                       "LJ 12-6 (cut/shift)",
                       "Coulomb(1/r) as mean field (c/s)",
                       "Coulomb(1/r) as mean field (cut only)",
                       "Yukawa(exp(-alpha r)/r (c/s)",
                       "Exponential(exp(-alpha r) (c/s)",
                       "Square Well",
                       "LJ 12-6 plus Yukawa(c/s)",
                       "r12 repulsion + Yukawa(c/s)",
                       "r18 repulsion + Yukawa(c/s)",
                       "rN repulsion plus Yukawa(c/s)")));
        RCP<ArrayStringValidator> U3D_ValidArray = rcp(new ArrayStringValidator(U3D_Valid));

        RCP<StringValidator> Vext1D_Valid = rcp(new StringValidator(tuple<std::string>(
                       "none",
                       "LJ 9-3 (cut/shift)",
                       "LJ 9-3 (v2) (c/s)",
                       "LJ 9-3 (no c/s)",
                       "LJ 9-3 shifted x",
                       "r9 repulsive (no c/s)",
                       "Exponential (no c/s)",
                       "Linear field (no c/s)",
                       "R7 + Yukawa (c/s)")));
        RCP<ArrayStringValidator> Vext1D_ValidArray = rcp(new ArrayStringValidator(Vext1D_Valid));

        RCP<StringValidator> SemiPerm_Valid = rcp(new StringValidator(tuple<std::string>(
                       "No semi-permeable surfaces", "All surfaces are semi-permeable", "Some surfaces are semi-permeable")));

        RCP<StringValidator> VextOption_Valid = rcp(new StringValidator(tuple<std::string>("none",
                       "1D Vext(x) potentials only", 
                       "Vext=U(r) only (fixed atoms or colloids)",
                       "Vext=U(r) and/or Vext=integral(U(r))dr", 
                       "Mixed 1D Vext(x) and 3D Vext(U(r)) types")));

     /* SET PARAMETERS */
        SurfaceInteraction_List->set("SI0.0: Define Wall-Fluid Interactions?",true, "Indicate if wall-fluid interactions are needed"); 

        Array<string> WFType_Array(Surface_List->get<int>("S3: Number of surface types"),"none");
        SurfaceInteraction_List->set("SI0.1: Wall-Fluid Computation Method",WFType_Array, "Method of computation for wall-fluid interactions", MethodType_ValidArray); 

        SurfaceInteraction_List->set("SI1: VEXT_MAX",20.0, "Set point for Vext inside surfaces. (VEXT_MAX typically between 10. and 20.)");
        SurfaceInteraction_List->set("SI2.0: Vext selections","none", "Indicate the general form for external fields, Vext", VextOption_Valid); 

        Array<string> VextType_Array(Surface_List->get<int>("S3: Number of surface types"),"Vext(x)");
        SurfaceInteraction_List->set("SI2.1: Type_vext",VextType_Array, "Specify type of external field for each surface", VextType_ValidArray); 

        Array<string> Vext1D_Array(Surface_List->get<int>("S3: Number of surface types"),"none");
        SurfaceInteraction_List->set("SI3: V(x) 1D",Vext1D_Array, "Specify analytic external field to apply for surfaces where 1D analytical potentials will be used", Vext1D_ValidArray);

        Array<string> Upair_Array(Surface_List->get<int>("S3: Number of surface types"),"none");
        SurfaceInteraction_List->set("SI4: Upair(r) for Vext(Upair(r))",Upair_Array, "Specify 3D pair potential to use for surfaces with atomic, colloidal, or numerically integrated external fields", U3D_ValidArray);

        SurfaceInteraction_List->set("SI5: Careful Boundaries?",true,"Set to true for careful integration at surface boundaries. This will allow for correct calculation of discontinuities in density profiles that arise from discontinuities in external fields. This applies for any hard core walls.");
        SurfaceInteraction_List->set("SI6: Semipermeable walls","No semi-permeable surfaces","set to true if any of the surface types are semipermeable membranes.",SemiPerm_Valid);

        TwoDArray<double> Vext_semiperm_Array( (Surface_List->get<int>("S3: Number of surface types"))*(Fluid_List->get<int>("F1_Ncomp")),SurfaceInteraction_List->get<double>("SI1: VEXT_MAX"));
        SurfaceInteraction_List->set("SI7: Vext_semiperm Array",Vext_semiperm_Array,"Set value of external field inside the semipermeable surface. A value of VEXT_MAX is used to indicate that the surface will not be permeable.");

  return;
}
/*************************************************************************************************************************************************/
void dft_GUI_vextType_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List)
{
  string tmpSt_1D[NWALL_MAX_TYPE],tmp_st;
  int Nvext1D=0, NvextPairPot=0,i,j,count_tmp=0;;
  bool TF_tmp=false;

     /* VALIDATORS*/

        RCP<StringValidator> MethodType_Valid = rcp(new StringValidator(tuple<std::string>(
                       "none",
                       "Pure exclusions / hard surface",
                       "V(r or x=dist. to surface)",
                       "V(r or x=dist. to center)",
                       "V=integral(U(r))")));

        RCP<ArrayStringValidator> MethodType_ValidArray = rcp(new ArrayStringValidator(MethodType_Valid));

        RCP<StringValidator> VextType_Valid = rcp(new StringValidator(tuple<std::string>(
                       "Vext(x)", "U(r) for Vext(U(r))")));
        RCP<ArrayStringValidator> VextType_ValidArray = rcp(new ArrayStringValidator(VextType_Valid));

        RCP<StringValidator> U3D_Valid = rcp(new StringValidator(tuple<std::string>(
                       "none",
                       "LJ 12-6 (cut/shift)",
                       "Coulomb(1/r) as mean field (c/s)",
                       "Coulomb(1/r) as mean field (cut only)",
                       "Yukawa(exp(-alpha r)/r (c/s)",
                       "Exponential(exp(-alpha r) (c/s)",
                       "Square Well",
                       "LJ 12-6 plus Yukawa(c/s)",
                       "r12 repulsion + Yukawa(c/s)",
                       "r18 repulsion + Yukawa(c/s)")));
        RCP<ArrayStringValidator> U3D_ValidArray = rcp(new ArrayStringValidator(U3D_Valid));

        RCP<StringValidator> Vext1D_Valid = rcp(new StringValidator(tuple<std::string>(
                       "none",
                       "LJ 9-3 (cut/shift)",
                       "LJ 9-3 (v2) (c/s)",
                       "LJ 9-3 (no c/s)",
                       "LJ 9-3 shifted x",
                       "r9 repulsive (no c/s)",
                       "Exponential (no c/s)",
                       "Linear field (no c/s)",
                       "R7 + Yukawa (c/s)")));
        RCP<ArrayStringValidator> Vext1D_ValidArray = rcp(new ArrayStringValidator(Vext1D_Valid));

        RCP<StringValidator> SemiPerm_Valid = rcp(new StringValidator(tuple<std::string>(
                       "No semi-permeable surfaces", "All surfaces are semi-permeable", "Some surfaces are semi-permeable")));

        RCP<StringValidator> VextOption_Valid = rcp(new StringValidator(tuple<std::string>("none",
                       "1D Vext(x) potentials only", 
                       "Vext=U(r) only (fixed atoms or colloids)",
                       "Vext=U(r) and/or Vext=integral(U(r))dr", 
                       "Mixed 1D Vext(x) and 3D Vext(U(r)) types")));


  if (Nwall >0) SurfaceInteraction_List->set("SI0.0: Define Wall-Fluid Interactions?",true, "Indicate if wall-fluid interactions are needed");
  else          SurfaceInteraction_List->set("SI0.0: Define Wall-Fluid Interactions?",false, "Indicate if wall-fluid interactions are needed");
  
  for (i=0;i<Nwall_type;i++){
     if (Ipot_wf_n[i]==VEXT_NONE) tmpSt_1D[i]="none";
     else if (Ipot_wf_n[i]==VEXT_HARD) tmpSt_1D[i]="Pure exclusions / hard surface";
     else if (Ipot_wf_n[i]==VEXT_DIST_TO_SURF) tmpSt_1D[i]="V(r or x=dist. to surface)";
     else if (Ipot_wf_n[i]==VEXT_DIST_TO_CENTER) tmpSt_1D[i]="V(r or x=dist. to center)";
     else if (Ipot_wf_n[i]==VEXT_3D_INTEGRATED){
         tmpSt_1D[i]="V=integral(U(r))";
         TF_tmp=true;
     }
  }
  Array<string> WFType_Array(tmpSt_1D,tmpSt_1D+Nwall_type);
  SurfaceInteraction_List->set("SI0.1: Wall-Fluid Computation Method",WFType_Array, "Method of computation for wall-fluid interactions", MethodType_ValidArray); 

  SurfaceInteraction_List->set("SI1: VEXT_MAX",VEXT_MAX, "Set point for Vext inside surfaces. (VEXT_MAX typically between 10. and 20.)");

  for (i=0;i<Nwall_type;i++){
     if (Type_vext[i]==VEXT_DEFINED){
         tmpSt_1D[i]="Vext(x)";
         Nvext1D++;
     }
     else if (Type_vext[i]==VEXT_PAIR_POTENTIAL){
          tmpSt_1D[i]="U(r) for Vext(U(r))";
          NvextPairPot++;
     }
  }
  if (Nvext1D==0 && NvextPairPot==0) tmp_st="none";
  else if (Nvext1D>0 && NvextPairPot==0) tmp_st="1D Vext(x) potentials only";
  else if (Nvext1D==0 && NvextPairPot>0){
      if (TF_tmp)    tmp_st="Vext=U(r) and/or Vext=integral(U(r))dr";
      else           tmp_st="Vext=U(r) only (fixed atoms or colloids)";
  } 
  else if (Nvext1D>0 && NvextPairPot>0) tmp_st="Mixed 1D Vext(x) and 3D Vext(U(r)) types";
  SurfaceInteraction_List->set("SI2.0: Vext selections",tmp_st, "Indicate the general form for external fields, Vext", VextOption_Valid); 

  Array<string> VextType_Array(tmpSt_1D,tmpSt_1D+Nwall_type);
  SurfaceInteraction_List->set("SI2.1: Type_vext",VextType_Array, "Specify type of external field for each surface", VextType_ValidArray); 

  for (i=0;i<Nwall_type;i++){
     if (Vext_PotentialID[i]==LJ9_3_CS) tmpSt_1D[i]="LJ 9-3 (cut/shift)";
     else if (Vext_PotentialID[i]==LJ9_3_v2_CS) tmpSt_1D[i]="LJ 9-3 (v2) (c/s)";
     else if (Vext_PotentialID[i]==LJ9_3_noCS) tmpSt_1D[i]="LJ 9-3 (no c/s)";
     else if (Vext_PotentialID[i]==LJ9_3_shiftX_CS) tmpSt_1D[i]="LJ 9-3 shifted x";
     else if (Vext_PotentialID[i]==REPULSIVE9_noCS) tmpSt_1D[i]="r9 repulsive (no c/s)";
     else if (Vext_PotentialID[i]==EXP_ATT_noCS) tmpSt_1D[i]="Exponential (no c/s)";
     else if (Vext_PotentialID[i]==LINEAR_noCS) tmpSt_1D[i]="Linear field (no c/s)";
     else if (Vext_PotentialID[i]==R7_YUKAWA_SUM_CS) tmpSt_1D[i]="R7 + Yukawa (c/s)";
     else tmpSt_1D[i]="none";
  }
  
  Array<string> Vext1D_Array(tmpSt_1D,tmpSt_1D+Nwall_type);
  SurfaceInteraction_List->set("SI3: V(x) 1D",Vext1D_Array, "Specify analytic external field to apply for surfaces where 1D analytical potentials will be used", Vext1D_ValidArray);

  for (i=0;i<Nwall_type;i++){
     if (Vext_PotentialID[i]==PAIR_LJ12_6_CS) tmpSt_1D[i]="LJ 12-6 (cut/shift)";
     else if (Vext_PotentialID[i]==PAIR_COULOMB_CS) tmpSt_1D[i]="Coulomb(1/r) as mean field (c/s)";
     else if (Vext_PotentialID[i]==PAIR_COULOMB) tmpSt_1D[i]="Coulomb(1/r) as mean field (cut only)";
     else if (Vext_PotentialID[i]==PAIR_YUKAWA_CS) tmpSt_1D[i]="Yukawa(exp(-alpha r)/r (c/s)";
     else if (Vext_PotentialID[i]==PAIR_EXP_CS) tmpSt_1D[i]="Exponential(exp(-alpha r) (c/s)";
     else if (Vext_PotentialID[i]==PAIR_SW) tmpSt_1D[i]="Square Well";
     else if (Vext_PotentialID[i]==PAIR_LJandYUKAWA_CS) tmpSt_1D[i]="LJ 12-6 plus Yukawa(c/s)";
     else if (Vext_PotentialID[i]==PAIR_r12andYUKAWA_CS) tmpSt_1D[i]="r12 repulsion + Yukawa(c/s)";
     else if (Vext_PotentialID[i]==PAIR_r18andYUKAWA_CS) tmpSt_1D[i]="r18 repulsion + Yukawa(c/s)";
     else tmpSt_1D[i]="none";
  }
  Array<string> Upair_Array(tmpSt_1D,tmpSt_1D+Nwall_type);
  SurfaceInteraction_List->set("SI4: Upair(r) for Vext(Upair(r))",Upair_Array, "Specify 3D pair potential to use for surfaces with atomic, colloidal, or numerically integrated external fields", U3D_ValidArray);

  if (Lhard_surf==TRUE) SurfaceInteraction_List->set("SI5: Careful Boundaries?",true,"Set to true for careful integration at surface boundaries. This will allow for correct calculation of discontinuities in density profiles that arise from discontinuities in external fields. This applies for any hard core walls.");
  else                  SurfaceInteraction_List->set("SI5: Careful Boundaries?",false,"Set to true for careful integration at surface boundaries. This will allow for correct calculation of discontinuities in density profiles that arise from discontinuities in external fields. This applies for any hard core walls.");

  TF_tmp=false; 
  count_tmp=0;
  for (i=0;i<Nwall_type;i++){
     for (j=0;j<Ncomp;j++){
        if (Lsemiperm[i][j]==TRUE){
          TF_tmp=true;
          count_tmp++;
        }
     }
  }

  if (!TF_tmp || count_tmp==0)
       SurfaceInteraction_List->set("SI6: Semipermeable walls","No semi-permeable surfaces","set to true if any of the surface types are semipermeable membranes.",SemiPerm_Valid);
  else if (TF_tmp && count_tmp==Nwall_type*Ncomp) 
       SurfaceInteraction_List->set("SI6: Semipermeable walls","All surfaces are semi-permeable","set to true if any of the surface types are semipermeable membranes.",SemiPerm_Valid);
  else if (TF_tmp && count_tmp!=Nwall_type*Ncomp) 
       SurfaceInteraction_List->set("SI6: Semipermeable walls","Some surfaces are semi-permeable","set to true if any of the surface types are semipermeable membranes.",SemiPerm_Valid);

  TwoDArray<double> Vext_semiperm_Array(Nwall_type,Ncomp);
  for (i=0; i<Nwall_type;i++) 
      for (j=0; j<Ncomp;j++){
         if (Lsemiperm[i][j]==TRUE) Vext_semiperm_Array(i,j)=Vext_membrane[i][j];
         else Vext_semiperm_Array(i,j)=VEXT_MAX;
      } 
  SurfaceInteraction_List->set("SI7: Vext_semiperm Array",Vext_semiperm_Array,"Set value of external field inside the semipermeable surface. A value of VEXT_MAX is used to indicate that the surface will not be permeable.");

  return;
}
/*************************************************************************************************************************************************/
void dft_GUI_vextType_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List)
{
  /* DEPENDENCIES */

  Dependency::ParameterEntryList NSurfaces_Deps;
  NSurfaces_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI0.1: Wall-Fluid Computation Method"));
  NSurfaces_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI1: VEXT_MAX"));
  NSurfaces_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI5: Careful Boundaries?"));
  NSurfaces_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI6: Semipermeable walls"));
  NSurfaces_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI2.0: Vext selections"));
  NSurfaces_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI2.1: Type_vext"));

  RCP<BoolVisualDependency> NsurfDeps = rcp(
     new BoolVisualDependency( SurfaceInteraction_List->getEntryRCP("SI0.0: Define Wall-Fluid Interactions?"), NSurfaces_Deps),true);

  RCP<StringVisualDependency> VextTypeArray_Dep = rcp(
     new StringVisualDependency(SurfaceInteraction_List->getEntryRCP("SI2.0: Vext selections"),
                                SurfaceInteraction_List->getEntryRCP("SI2.1: Type_vext"), 
                                tuple<std::string>("Mixed 1D Vext(x) and 3D Vext(U(r)) types"),true));

  RCP<StringVisualDependency> Vext1D_Dep = rcp(
     new StringVisualDependency(SurfaceInteraction_List->getEntryRCP("SI2.0: Vext selections"),
                                SurfaceInteraction_List->getEntryRCP("SI3: V(x) 1D"), 
                                tuple<std::string>( "1D Vext(x) potentials only","Mixed 1D Vext(x) and 3D Vext(U(r)) types"),true));

  RCP<StringVisualDependency> VextU3D_Dep = rcp(
     new StringVisualDependency( SurfaceInteraction_List->getEntryRCP("SI2.0: Vext selections"),
                                 SurfaceInteraction_List->getEntryRCP("SI4: Upair(r) for Vext(Upair(r))"),
                                 tuple<std::string>( "Vext=U(r) only (fixed atoms or colloids)", 
                                  "Vext=U(r) and/or Vext=integral(U(r))dr", 
                                  "Mixed 1D Vext(x) and 3D Vext(U(r)) types"),true));

   RCP<StringVisualDependency> SemiPerm_Dep = rcp( new StringVisualDependency( 
                                      SurfaceInteraction_List->getEntryRCP("SI6: Semipermeable walls"), 
                                      SurfaceInteraction_List->getEntryRCP("SI7: Vext_semiperm Array"),
                                      tuple<std::string>("Some surfaces are semi-permeable","All surfaces are semi-permeable")));

   Dependency::ParameterEntryList ArrayLengthSTR_Deps;
   ArrayLengthSTR_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI0.1: Wall-Fluid Computation Method"));
   ArrayLengthSTR_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI2.1: Type_vext"));
   ArrayLengthSTR_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI3: V(x) 1D"));
   ArrayLengthSTR_Deps.insert(SurfaceInteraction_List->getEntryRCP("SI4: Upair(r) for Vext(Upair(r))"));
   RCP<NumberArrayLengthDependency<int,string> > ArrayLengthSTR_Dep = rcp(
             new NumberArrayLengthDependency<int,string>(Surface_List->getEntryRCP("S3: Number of surface types"), ArrayLengthSTR_Deps));

   RCP<TwoDColDependency<int,double> > VextMembraneColNumber_Dep = rcp(
             new TwoDColDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"),SurfaceInteraction_List->getEntryRCP("SI7: Vext_semiperm Array")));
   RCP<TwoDRowDependency<int,double> > VextMembraneRowNumber_Dep = rcp(
             new TwoDRowDependency<int,double>(Surface_List->getEntryRCP("S3: Number of surface types"),SurfaceInteraction_List->getEntryRCP("SI7: Vext_semiperm Array")));



  /* DEPENDENCY SHEET ENTRIES*/
  depSheet_Tramonto->addDependency(NsurfDeps);
  depSheet_Tramonto->addDependency(VextTypeArray_Dep);
  depSheet_Tramonto->addDependency(Vext1D_Dep);
  depSheet_Tramonto->addDependency(VextU3D_Dep);
  depSheet_Tramonto->addDependency(SemiPerm_Dep);
  depSheet_Tramonto->addDependency(ArrayLengthSTR_Dep);
  depSheet_Tramonto->addDependency(VextMembraneColNumber_Dep);
  depSheet_Tramonto->addDependency(VextMembraneRowNumber_Dep);

  return;
}
/*************************************************************************************************************************************************/
