using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_Continuation(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                     Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                     Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                     Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List,
                     Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                     Teuchos::RCP<Teuchos::ParameterList> Continuation_List)
{
  bool set_defaults_from_old_format_file=true;
  int i;
  int tmp_1D[2];
  string tmp_string, cont_type_str,cont_param_str, cont_paramDep_str, cont_mesh_str;

        /**************************************/
        /* Define help strings for this section.*/
        /**************************************/
        cont_type_str="Select the continuation type of interest.";
        cont_param_str="Select the independent parameter for continuation studies.";
        cont_paramDep_str="Select the continuation parameter type for the dependent parameter in binodal studies";
        cont_mesh_str="Indicate the direction where the number of nodes will be varied with continuation";

        /**************************************/
        /* Define validators for this section.*/
        /**************************************/

    RCP<EnhancedNumberValidator<int> > NIDValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));
    RCP<EnhancedNumberValidator<int> > Dimension_Validator = rcp(new EnhancedNumberValidator<int>(0,2,1));

    RCP<StringValidator> NodeChange_Validator = rcp(
           new StringValidator(tuple<std::string>("center","Left(dim=0), Bottom(dim=1), or Back(dim=2) boundary",
                                     "Right(dim=0), Top(dim=1), or Front(dim=2) boundary")));

    RCP<StringValidator> ContType_Validator = rcp(
           new StringValidator(tuple<std::string>("None","LOCA: Simple Parameter Continuation",
                                     "LOCA: Arc-Length Parameter Continuation",
                                     "LOCA: Spinodal Continuation",
                                     "LOCA: Binodal Continuation",
                                     "Mesh Continuation",
                                     "Mesh Stepping with LOCA Binodal")));

    RCP<StringValidator> ContParam_Validator = rcp(
           new StringValidator(tuple<std::string>("None","Temperature","Density (1 species)","Chemical Potential (1 species)",
                               "Wall-Wall energy param (iwall_type or iwall_type,jwall_type pair)", "Wall-Fluid energy param (iwall_type,icomp)", 
                               "Fluid-fluid energy param (i or ij pair)","Electrostatic parameter (iwall)","All Electrostatic parameters (Nwalls)",
                               "Vext_membrane(iwall_type,icomp)","Fluid size: SigmaFF(i or ij pair)")));


        /*********************/
        /* set up parameters */
        /*********************/

        Continuation_List->set("C1: Continuation Type","None", cont_type_str,ContType_Validator);
        Continuation_List->set("C2: Continuation Parameter","None",cont_param_str,ContParam_Validator);

        Continuation_List->set("C3: icomp",0,"ID of Fluid Species whose parameter will be varied in continuation studies.");
        Continuation_List->set("C3: iwall_type",0,"ID for Wall Type whose parameter will be varied in continuation studies.\n Note that all walls of type iwall_type will be affected.");
        Continuation_List->set("C3: iwall",0,"ID for Wall whose parameter will be varied in continuation studies.");
        Array<int> IJcomp_Array(2,0);
        Continuation_List->set("C3: icomp,jcomp",IJcomp_Array,"i,j IDs of the particular pair interaction to be varied in continuation studies.");
        Array<int> IwallIcomp_Array(2,0);
        Continuation_List->set("C3: iwall_type,icomp",IwallIcomp_Array,"iwall_type,icomp IDs of the wall-fluid interaction pair to be varied in continuation studies.");
        Array<int> IJwall_Array(2,0);
        Continuation_List->set("C3: iwall_type,jwall_type",IJwall_Array,"iwall_type,jwall_type IDs of the wall-wall particular pair interaction to be varied in continuation studies.");

        Continuation_List->set("C4: Continuation Step Size",0.01,"Set the initial step size for independent parameter in continuation studies.");
        Continuation_List->set("C5: Number of Steps",1,"Enter the number of continuation steps desired in the study.");
        Continuation_List->set("C6: Step Aggressivnes",0.01,"Set to 0.0 for a constant step.  Increase for more rapid variations in the step size.");

        Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","None",cont_paramDep_str,ContParam_Validator);

        Continuation_List->set("C9.1: Direction of Mesh Cont.",0,"Indicate the direction where the number of nodes will be varied with continuation",Dimension_Validator);
        Continuation_List->set("C9.2: Location for node count change","center",cont_mesh_str,NodeChange_Validator);

        Continuation_List->set("C8: icomp",0,"ID of Fluid species of dependent parameter in binodal studies.");
        Continuation_List->set("C8: iwall_type",0,"ID for Wall Type of the dependent parameter in binodal.\n Note that all walls of type iwall_type will be affected.");
        Continuation_List->set("C8: iwall",0,"ID for Wall property to be used as dependent parameter in binodal studies.");
        Array<int> IJcomp2_Array(2,0);
        Continuation_List->set("C8: icomp,jcomp",IJcomp2_Array,"i,j IDs of the particular pair interaction to be varied as the dependent parameter in binodal studies.");
        Array<int> IwallIcomp2_Array(2,0);
        Continuation_List->set("C8: iwall_type,icomp",IwallIcomp2_Array,"iwall_type,icomp IDs of the wall-fluid interaction pair to be varied as the dependent parameter in binodal studies.");
        Array<int> IJwall2_Array(2,0);
        Continuation_List->set("C8: iwall_type,jwall_type",IJwall2_Array,"iwall_type,jwall_type IDs of the wall-wall pair interaction to be varied as the dependent parameter in continuation studies.");

    if (set_defaults_from_old_format_file){
        if (Loca.method==-1){
          if(Nruns==1) Continuation_List->set("C1: Continuation Type","None", cont_type_str,ContType_Validator);
          else          Continuation_List->set("C1: Continuation Type","Mesh Continuation", cont_type_str,ContType_Validator);
        }
        else if (Loca.method==0 || Loca.method==1){
          Continuation_List->set("C1: Continuation Type","LOCA: Simple Parameter Continuation", cont_type_str,ContType_Validator);
        }
        else if (Loca.method==2){
          Continuation_List->set("C1: Continuation Type","LOCA: Arc-Length Parameter Continuation", cont_type_str,ContType_Validator);
        }
        else if (Loca.method==3){
          Continuation_List->set("C1: Continuation Type","LOCA: Spinodal Continuation", cont_type_str,ContType_Validator);
        }
        else if (Loca.method==4){
          if (Nruns >1) Continuation_List->set("C1: Continuation Type","Mesh Stepping with LOCA Binodal", cont_type_str,ContType_Validator);
          else           Continuation_List->set("C1: Continuation Type","LOCA: Binodal Continuation", cont_type_str,ContType_Validator);
        }

        if (Loca.method==-1 && Nruns==1)Continuation_List->set("C2: Continuation Parameter","None",cont_param_str,ContParam_Validator);
        else{
           if (Loca.cont_type1==CONT_MESH || Nruns>1)  Continuation_List->set("C2: Continuation Parameter","None",cont_param_str,ContParam_Validator);
           else if (Loca.cont_type1==CONT_TEMP) Continuation_List->set("C2: Continuation Parameter","Temperature", cont_param_str,ContParam_Validator);
           else if (Loca.cont_type1==CONT_RHO_I){
                Continuation_List->set("C2: Continuation Parameter","Density (1 species)", cont_param_str,ContParam_Validator);
                Continuation_List->set("C3: icomp",Cont_ID[0][0],"ID of Fluid Species whose parameter will be varied in continuation studies.");
           }
           else if (Loca.cont_type1==CONT_BETAMU_I){
                Continuation_List->set("C2: Continuation Parameter","Chemical Potential (1 species)", cont_param_str,ContParam_Validator);
                Continuation_List->set("C3: icomp",Cont_ID[0][0],"ID of Fluid Species whose parameter will be varied in continuation studies.");
           }
           else if (Loca.cont_type1==CONT_EPSW_I){
                Continuation_List->set("C2: Continuation Parameter","Fluid-fluid energy param (i or ij pair)", cont_param_str,ContParam_Validator);
                Continuation_List->set("C3: iwall_type",Cont_ID[0][0],
                       "ID for Wall Type whose parameter will be varied in continuation studies.\n Note that all walls of type iwall_type will be affected.");
           }
           else if (Loca.cont_type1==CONT_EPSWF_IJ || Loca.cont_type1==CONT_SEMIPERM_IJ){
                if (Loca.cont_type1==CONT_EPSWF_IJ) Continuation_List->set("C2: Continuation Parameter","Wall-Fluid energy param (iwall_type,icomp)", cont_param_str,ContParam_Validator);
                else                               Continuation_List->set("C2: Continuation Parameter","Vext_membrane(iwall_type,icomp)",cont_param_str,ContParam_Validator);
                for(i=0;i<2;i++) tmp_1D[i]=Cont_ID[0][i];
                Array<int> IwallIcomp_Array(tmp_1D,tmp_1D+2);
                Continuation_List->set("C3: iwall_type,icomp",IwallIcomp_Array,"iwall_type,icomp IDs of the wall-fluid interaction pair to be varied in continuation studies.");
           }
           else if (Loca.cont_type1==CONT_EPSFF_IJ || Loca.cont_type1==CONT_SIGMAFF_IJ){
                if (Loca.cont_type1==CONT_EPSFF_IJ) Continuation_List->set("C2: Continuation Parameter","Fluid-fluid energy param (i or ij pair)", cont_param_str,ContParam_Validator);
                else   Continuation_List->set("C2: Continuation Parameter","Fluid size: SigmaFF(i or ij pair)",cont_param_str,ContParam_Validator);
                for(i=0;i<2;i++) tmp_1D[i]=Cont_ID[0][i];
                Array<int> IJcomp_Array(tmp_1D,tmp_1D+2);
                Continuation_List->set("C3: icomp,jcomp",IJcomp_Array,"i,j IDs of the particular pair interaction to be varied in continuation studies.");
           }
 
           else if (Loca.cont_type1==CONT_ELECPARAM_I){
                Continuation_List->set("C2: Continuation Parameter","Electrostatic parameter (iwall)", cont_param_str,ContParam_Validator);
                Continuation_List->set("C3: iwall",Cont_ID[0][0],"ID for Wall whose parameter will be varied in continuation studies.");
           }
           else if (Loca.cont_type1==CONT_ELECPARAM_ALL){
                Continuation_List->set("C2: Continuation Parameter","All Electrostatic parameters (Nwalls)", cont_param_str,ContParam_Validator);
           }
           else Continuation_List->set("C2: Continuation Parameter","None",cont_param_str,ContParam_Validator);
        }

        if (Nruns>1){
           Continuation_List->set("C9.1: Direction of Mesh Cont.",Plane_new_nodes,"Indicate the direction where the number of nodes will be varied with continuation",Dimension_Validator);
           if (Pos_new_nodes==0)       Continuation_List->set("C9.2: Location for node count change","center",cont_mesh_str,NodeChange_Validator);
           else if (Pos_new_nodes==-1) Continuation_List->set("C9.2: Location for node count change","Left(dim=0), Bottom(dim=1), or Back(dim=2) boundary",cont_mesh_str,NodeChange_Validator);
           else if (Pos_new_nodes==1)  Continuation_List->set("C9.2: Location for node count change","Right(dim=0), Top(dim=1), or Front(dim=2) boundary",cont_mesh_str,NodeChange_Validator);

           Continuation_List->set("C4: Continuation Step Size",Del_1[Plane_new_nodes],"Set the initial step size for independent parameter in continuation studies.");
           Continuation_List->set("C5: Number of Steps",Nruns,"Enter the number of continuation steps desired in the study.");
        }
        else{
           Continuation_List->set("C4: Continuation Step Size",Loca.step_size,"Set the initial step size for independent parameter in continuation studies.");
           Continuation_List->set("C5: Number of Steps",Loca.num_steps,"Enter the number of continuation steps desired in the study.");
           Continuation_List->set("C6: Step Aggressivnes",Loca.aggr,"Set to 0.0 for a constant step.  Increase for more rapid variations in the step size.");
        }

        if (Lbinodal){
           Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","None",cont_paramDep_str,ContParam_Validator);

           if (Loca.cont_type2==CONT_MESH)  Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","None",cont_paramDep_str,ContParam_Validator);
           else if (Loca.cont_type2==CONT_TEMP) Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Temperature", cont_paramDep_str,ContParam_Validator);
           else if (Loca.cont_type2==CONT_RHO_I){
                Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Density (1 species)", cont_paramDep_str,ContParam_Validator);
                Continuation_List->set("C8: icomp",Cont_ID[1][0],"ID of Fluid Species whose parameter will be varied in continuation studies.");
           }
           else if (Loca.cont_type2==CONT_BETAMU_I){
                Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Chemical Potential (1 species)", cont_paramDep_str,ContParam_Validator);
                Continuation_List->set("C8: icomp",Cont_ID[1][0],"ID of Fluid Species whose parameter will be varied in continuation studies.");
           }
           else if (Loca.cont_type2==CONT_EPSW_I){
                Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Fluid-fluid energy param (i or ij pair)", cont_paramDep_str,ContParam_Validator);
                Continuation_List->set("C8: iwall_type",Cont_ID[1][0],
                       "ID for Wall Type whose parameter will be varied in continuation studies.\n Note that all walls of type iwall_type will be affected.");
           }
           else if (Loca.cont_type2==CONT_EPSWF_IJ || Loca.cont_type2==CONT_SEMIPERM_IJ){
                if (Loca.cont_type2==CONT_EPSWF_IJ) Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Wall-Fluid energy param (iwall_type,icomp)", cont_paramDep_str,ContParam_Validator);
                else                               Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Vext_membrane(iwall_type,icomp)",cont_paramDep_str,ContParam_Validator);
                for(i=0;i<2;i++) tmp_1D[i]=Cont_ID[1][i];
                Array<int> IwallIcomp_Array(tmp_1D,tmp_1D+2);
                Continuation_List->set("C8: iwall_type,icomp",IwallIcomp_Array,"iwall_type,icomp IDs of the wall-fluid interaction pair to be varied in continuation studies.");
           }
           else if (Loca.cont_type2==CONT_EPSFF_IJ || Loca.cont_type2==CONT_SIGMAFF_IJ){
                if (Loca.cont_type2==CONT_EPSFF_IJ) Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Fluid-fluid energy param (i or ij pair)", cont_paramDep_str,ContParam_Validator);
                else   Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Fluid size: SigmaFF(i or ij pair)",cont_paramDep_str,ContParam_Validator);
                for(i=0;i<2;i++) tmp_1D[i]=Cont_ID[1][i];
                Array<int> IJcomp_Array(tmp_1D,tmp_1D+2);
                Continuation_List->set("C8: icomp,jcomp",IJcomp_Array,"i,j IDs of the particular pair interaction to be varied in continuation studies.");
           }
 
           else if (Loca.cont_type2==CONT_ELECPARAM_I){
                Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","Electrostatic parameter (iwall)", cont_paramDep_str,ContParam_Validator);
                Continuation_List->set("C8: iwall",Cont_ID[1][0],"ID for Wall whose parameter will be varied in continuation studies.");
           }
           else if (Loca.cont_type2==CONT_ELECPARAM_ALL){
                Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","All Electrostatic parameters (Nwalls)", cont_paramDep_str,ContParam_Validator);
           }
           else Continuation_List->set("C7: Dependent Cont. Param. (Binodals)","None",cont_paramDep_str,ContParam_Validator);

        }
  
    }


      /************************/
      /* set up dependencies */
      /************************/

    RCP<StringCondition> ContinuationNoMesh2Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C1: Continuation Type"),
                        tuple<std::string>("LOCA: Simple Parameter Continuation",
                             "LOCA: Arc-Length Parameter Continuation", "LOCA: Spinodal Continuation",
                             "LOCA: Binodal Continuation","Mesh Stepping with LOCA Binodal")));

    RCP<StringCondition> ContinuationNoMeshCon = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C1: Continuation Type"),
                        tuple<std::string>("LOCA: Simple Parameter Continuation",
                             "LOCA: Arc-Length Parameter Continuation", "LOCA: Spinodal Continuation",
                             "LOCA: Binodal Continuation")));

    RCP<StringCondition> Cont_I_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C2: Continuation Parameter"),
                         tuple<std::string>("Density (1 species)","Chemical Potential (1 species)")));

    RCP<StringCondition> Cont_IW_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C2: Continuation Parameter"),
                         tuple<std::string>("Electrostatic parameter (iwall)")));

    RCP<StringCondition> Cont_EpsFF_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C2: Continuation Parameter"),
                         tuple<std::string>("Fluid-fluid energy param (i or ij pair)")));

    RCP<StringCondition> Cont_EpsWW_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C2: Continuation Parameter"),
                         tuple<std::string>("Wall-Wall energy param (iwall_type or iwall_type,jwall_type pair)")));

    RCP<StringCondition> Cont_EpsWF_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C2: Continuation Parameter"),
                         tuple<std::string>("Wall-Fluid energy param (iwall_type,icomp)","Vext_membrane(iwall_type,icomp)")));

    RCP<StringCondition> Cont_IMix_Con = rcp(
           new StringCondition(PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"),"Lorentz-Berthlot Mixing"),true);

    RCP<StringCondition> Cont_IMix_MANCon = rcp(
           new StringCondition(PotentialsFF_List->getEntryRCP("PF0_Off_Diagonal_Definitions"),"Manual Definition"));

    RCP<StringCondition> Binodal_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C1: Continuation Type"),
                         tuple<std::string>("LOCA: Binodal Continuation","Mesh Stepping with LOCA Binodal")));

    Condition::ConstConditionList Icomp_ANDFconList=tuple<RCP<const Condition> >(Cont_IMix_Con,Cont_EpsFF_Con);
       RCP<AndCondition> C3Icomp_EpsFCon = rcp(new AndCondition(Icomp_ANDFconList));
    Condition::ConstConditionList Icomp_ORconList=tuple<RCP<const Condition> >(Cont_I_Con,C3Icomp_EpsFCon);
       RCP<OrCondition> C3Icomp_ORCon = rcp(new OrCondition(Icomp_ORconList));
    Condition::ConstConditionList Icomp_conList=tuple<RCP<const Condition> >(ContinuationNoMeshCon,C3Icomp_ORCon);
       RCP<AndCondition> C3Icomp_Con = rcp(new AndCondition(Icomp_conList));
       RCP<ConditionVisualDependency> C3IcompVis_Dep = rcp(
          new ConditionVisualDependency(C3Icomp_Con, Continuation_List->getEntryRCP("C3: icomp")));


    Condition::ConstConditionList IW_conList=tuple<RCP<const Condition> >(ContinuationNoMeshCon,Cont_IW_Con);
    RCP<AndCondition> C3IW_Con = rcp(new AndCondition(IW_conList));
      RCP<ConditionVisualDependency> C3IWVis_Dep = rcp(
          new ConditionVisualDependency(C3IW_Con, Continuation_List->getEntryRCP("C3: iwall")));


    Condition::ConstConditionList EpsFF_conList=tuple<RCP<const Condition> >(ContinuationNoMeshCon,Cont_EpsFF_Con,Cont_IMix_MANCon);
    RCP<AndCondition> EpsFF_Con = rcp(new AndCondition(EpsFF_conList));
      RCP<ConditionVisualDependency> C3EpsFF_Dep = rcp(
          new ConditionVisualDependency(EpsFF_Con, Continuation_List->getEntryRCP("C3: icomp,jcomp")));

    Condition::ConstConditionList EpsWW_conList=tuple<RCP<const Condition> >(ContinuationNoMeshCon,Cont_EpsWW_Con,Cont_IMix_MANCon);
    RCP<AndCondition> EpsWW_Con = rcp(new AndCondition(EpsWW_conList));
      RCP<ConditionVisualDependency> C3EpsWW_Dep = rcp(
          new ConditionVisualDependency(EpsWW_Con, Continuation_List->getEntryRCP("C3: iwall_type,jwall_type")));

    Condition::ConstConditionList EpsW_conList=tuple<RCP<const Condition> >(ContinuationNoMeshCon,Cont_EpsWW_Con,Cont_IMix_Con);
    RCP<AndCondition> EpsW_Con = rcp(new AndCondition(EpsW_conList));
      RCP<ConditionVisualDependency> C3EpsW_Dep = rcp(
          new ConditionVisualDependency(EpsW_Con, Continuation_List->getEntryRCP("C3: iwall_type")));

    Condition::ConstConditionList EpsWF_conList=tuple<RCP<const Condition> >(ContinuationNoMeshCon,Cont_EpsWF_Con);
    RCP<AndCondition> EpsWF_Con = rcp(new AndCondition(EpsWF_conList));
      RCP<ConditionVisualDependency> C3EpsWF_Dep = rcp(
          new ConditionVisualDependency(EpsWF_Con, Continuation_List->getEntryRCP("C3: iwall_type,icomp")));

   Dependency::ParameterEntryList ContinuationParamDeps;
   ContinuationParamDeps.insert(Continuation_List->getEntryRCP("C4: Continuation Step Size"));
   ContinuationParamDeps.insert(Continuation_List->getEntryRCP("C5: Number of Steps"));

   RCP<StringVisualDependency> ContParams_Dep = rcp(
       new StringVisualDependency(Continuation_List->getEntryRCP("C1: Continuation Type"),ContinuationParamDeps,
           tuple<std::string>("LOCA: Simple Parameter Continuation",
                             "LOCA: Arc-Length Parameter Continuation", "LOCA: Spinodal Continuation",
                             "LOCA: Binodal Continuation", "Mesh Continuation",
                             "Mesh Stepping with LOCA Binodal")));

   RCP<StringVisualDependency> ContParamMesh_Dep = rcp(
       new StringVisualDependency(Continuation_List->getEntryRCP("C1: Continuation Type"),Continuation_List->getEntryRCP("C2: Continuation Parameter"),
           tuple<std::string>("LOCA: Arc-Length Parameter Continuation", "LOCA: Spinodal Continuation","LOCA: Binodal Continuation","LOCA: Simple Parameter Continuation")));

   RCP<StringVisualDependency> Aggr_Dep = rcp(
       new StringVisualDependency(Continuation_List->getEntryRCP("C1: Continuation Type"),Continuation_List->getEntryRCP("C6: Step Aggressivnes"),
           tuple<std::string>("LOCA: Arc-Length Parameter Continuation", "LOCA: Spinodal Continuation","LOCA: Binodal Continuation","LOCA: Simple Parameter Continuation")));

   Dependency::ParameterEntryList MeshContDeps;
   MeshContDeps.insert(Continuation_List->getEntryRCP("C9.1: Direction of Mesh Cont."));
   MeshContDeps.insert(Continuation_List->getEntryRCP("C9.2: Location for node count change"));
   RCP<StringVisualDependency> MeshCont_Dep = rcp(
       new StringVisualDependency(Continuation_List->getEntryRCP("C1: Continuation Type"),MeshContDeps,
           tuple<std::string>("Mesh Continuation","Mesh Stepping with LOCA Binodal")));

   RCP<ConditionVisualDependency> Binodal_Dep = rcp(
       new ConditionVisualDependency(Binodal_Con, Continuation_List->getEntryRCP("C7: Dependent Cont. Param. (Binodals)")));

    RCP<StringCondition> Cont2_I_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C7: Dependent Cont. Param. (Binodals)"),
                         tuple<std::string>("Density (1 species)","Chemical Potential (1 species)")));

    RCP<StringCondition> Cont2_IW_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C7: Dependent Cont. Param. (Binodals)"),
                         tuple<std::string>("Electrostatic parameter (iwall)")));

    RCP<StringCondition> Cont2_EpsFF_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C7: Dependent Cont. Param. (Binodals)"),
                         tuple<std::string>("Fluid-fluid energy param (i or ij pair)")));

    RCP<StringCondition> Cont2_EpsWW_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C7: Dependent Cont. Param. (Binodals)"),
                         tuple<std::string>("Wall-Wall energy param (iwall_type or iwall_type,jwall_type pair)")));

    RCP<StringCondition> Cont2_EpsWF_Con = rcp(
           new StringCondition(Continuation_List->getEntryRCP("C7: Dependent Cont. Param. (Binodals)"),
                         tuple<std::string>("Wall-Fluid energy param (iwall_type,icomp)","Vext_membrane(iwall_type,icomp)")));


    Condition::ConstConditionList Bin_Icomp_ANDFconList=tuple<RCP<const Condition> >(Cont_IMix_Con,Cont2_EpsFF_Con);
       RCP<AndCondition> C8Icomp_EpsFCon = rcp(new AndCondition(Bin_Icomp_ANDFconList));
    Condition::ConstConditionList BinIcomp_ORconList=tuple<RCP<const Condition> >(Cont2_I_Con,C8Icomp_EpsFCon);
       RCP<OrCondition> C8Icomp_ORCon = rcp(new OrCondition(BinIcomp_ORconList));
    Condition::ConstConditionList Icomp_BinconList=tuple<RCP<const Condition> >(ContinuationNoMesh2Con,C8Icomp_ORCon,Binodal_Con);
       RCP<AndCondition> C8Icomp_Con = rcp(new AndCondition(Icomp_BinconList));
       RCP<ConditionVisualDependency> C8IcompVis_Dep = rcp(
          new ConditionVisualDependency(C8Icomp_Con, Continuation_List->getEntryRCP("C8: icomp")));


    Condition::ConstConditionList IW_BinconList=tuple<RCP<const Condition> >(ContinuationNoMesh2Con,Cont2_IW_Con,Binodal_Con);
    RCP<AndCondition> C8IW_Con = rcp(new AndCondition(IW_BinconList));
      RCP<ConditionVisualDependency> C8IWVis_Dep = rcp(
          new ConditionVisualDependency(C8IW_Con, Continuation_List->getEntryRCP("C8: iwall")));


    Condition::ConstConditionList EpsFF_BinconList=tuple<RCP<const Condition> >(ContinuationNoMesh2Con,Cont2_EpsFF_Con,Cont_IMix_MANCon,Binodal_Con);
    RCP<AndCondition> EpsFF_BinCon = rcp(new AndCondition(EpsFF_BinconList));
      RCP<ConditionVisualDependency> C8EpsFF_Dep = rcp(
          new ConditionVisualDependency(EpsFF_BinCon, Continuation_List->getEntryRCP("C8: icomp,jcomp")));


    Condition::ConstConditionList EpsWW_BinconList=tuple<RCP<const Condition> >(ContinuationNoMesh2Con,Cont2_EpsWW_Con,Cont_IMix_MANCon,Binodal_Con);
    RCP<AndCondition> EpsWW_BinCon = rcp(new AndCondition(EpsWW_BinconList));
      RCP<ConditionVisualDependency> C8EpsWW_Dep = rcp(
          new ConditionVisualDependency(EpsWW_BinCon, Continuation_List->getEntryRCP("C8: iwall_type,jwall_type")));

    Condition::ConstConditionList EpsW_BinconList=tuple<RCP<const Condition> >(ContinuationNoMesh2Con,Cont2_EpsWW_Con,Cont_IMix_Con,Binodal_Con);
    RCP<AndCondition> EpsW_BinCon = rcp(new AndCondition(EpsW_BinconList));
      RCP<ConditionVisualDependency> C8EpsW_Dep = rcp(
          new ConditionVisualDependency(EpsW_BinCon, Continuation_List->getEntryRCP("C8: iwall_type")));

    Condition::ConstConditionList EpsWF_BinconList=tuple<RCP<const Condition> >(ContinuationNoMesh2Con,Cont2_EpsWF_Con,Binodal_Con);
    RCP<AndCondition> EpsWF_BinCon = rcp(new AndCondition(EpsWF_BinconList));
      RCP<ConditionVisualDependency> C8EpsWF_Dep = rcp(
          new ConditionVisualDependency(EpsWF_BinCon, Continuation_List->getEntryRCP("C8: iwall_type,icomp")));


      /*****************************************/
      /* add the dependencies for this section.*/
      /*****************************************/
      depSheet_Tramonto->addDependency(ContParams_Dep);
      depSheet_Tramonto->addDependency(Aggr_Dep);
      depSheet_Tramonto->addDependency(ContParamMesh_Dep);
      depSheet_Tramonto->addDependency(C3IcompVis_Dep);
      depSheet_Tramonto->addDependency(C3IWVis_Dep);
      depSheet_Tramonto->addDependency(C3EpsW_Dep);
      depSheet_Tramonto->addDependency(C3EpsFF_Dep);
      depSheet_Tramonto->addDependency(C3EpsWW_Dep);
      depSheet_Tramonto->addDependency(C3EpsWF_Dep);

      depSheet_Tramonto->addDependency(Binodal_Dep);
      depSheet_Tramonto->addDependency(C8IcompVis_Dep);
      depSheet_Tramonto->addDependency(C8IWVis_Dep);
      depSheet_Tramonto->addDependency(C8EpsFF_Dep);
      depSheet_Tramonto->addDependency(C8EpsWW_Dep);
      depSheet_Tramonto->addDependency(C8EpsW_Dep);
      depSheet_Tramonto->addDependency(C8EpsWF_Dep);
      depSheet_Tramonto->addDependency(MeshCont_Dep);

  return;
}

