using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_NumericalMethods(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                         Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                         Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                         Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List)
{
    /****************************************************************************************************************/
  /****************************** FUNCTIONAL CONTROL PARAMETER SECTION ********************************************/
  /****************************************************************************************************************/
    /**************************************/
    /* Define validators for this section.*/
    /**************************************/
    RCP<StringValidator> NLSolverTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Newton Built-In","Newton LOCA", 
						  "Picard Built-In","Picard LOCA",
						  "Picard/Newton Built-In","Picard/Newton LOCA")));

    RCP<EnhancedNumberValidator<int> > Niter_Validator = rcp(new EnhancedNumberValidator<int>());

    RCP<EnhancedNumberValidator<double> > NLtolValidator = rcp(new EnhancedNumberValidator<double>(1.e-10,1.e-2,1.e-6,8));

    RCP<EnhancedNumberValidator<double> > NLupdate_Validator = rcp(new EnhancedNumberValidator<double>(0.001,1.0,0.01));

    RCP<StringValidator> LoadBalValidator = rcp(
           new StringValidator(tuple<std::string>("Linear Matrix Balance","Geometric Recursive Bisection","Weighted Recursive Bisection")));

    RCP<StringValidator> LSolverTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Schur Solver","General/AztecOO solver")));

    RCP<StringValidator> LSolverApproachValidator = rcp(
           new StringValidator(tuple<std::string>("GMRES","cg","tfqmr","cg2","bicgstab")));

    RCP<StringValidator> LScalingValidator = rcp(
           new StringValidator(tuple<std::string>("none","row_sum","jacobi","symrow_sum")));

    RCP<StringValidator> LPreconditionerValidator = rcp(
           new StringValidator(tuple<std::string>("none","ilu","ilut","Jacobi","symmetric Gauss-Seidel","LSpoly3")));

    RCP<EnhancedNumberValidator<int> > NlevelILUT_Validator = rcp(new EnhancedNumberValidator<int>());

    RCP<StringValidator> CoarsenTypeValidator = rcp(
         new StringValidator(tuple<std::string>("none",
        "Jacobian Coarsening consistent with Residual Coarsening",
	"Jacobian: factor of 2 in most refined zone",
	"Jacobian: factor of 2 in all but most coarse zone",
	"Jacobian: use most coarse zone to define entire matrix",
	"Jacobian: use 2nd most coarse zone in all but most coarse region",
	"Jacobian: use a specified mesh spacing (Esize_jac) for matrix calculations.")));


    /***************************************************************/
    /* set up parameters */
    /***************************************************************/
    Solver_List->set("S1: Load Balancing Approach", "Weighted Recursive Bisection", "Select a method for load balancing a parallel processing job", LoadBalValidator);

    Coarsening_List->set("C1: residual coarsening?",false,"Set to true for coarsening of residual equations in zones away from the surfaces.\n  The approach increases the mesh size by a factor of 2 in all dimensions in each successive\n user defined zone based on distance from surfaces.");
    Coarsening_List->set("C2: jacobian coarsening?",false,"Set to true for coarsening of jacobian integrals based on the\n distance of a given node point from the surfaces in the problem.");
    Coarsening_List->set("C3: Type of Jacobian coarsening","none","Select a method for coarsening the mesh or jacobian",CoarsenTypeValidator);
    Coarsening_List->set("C4: Number of coarsening zones",1,"Set number of coarsening zones in problem.");
    Array<double> Rmax_zone_Array( (Coarsening_List->get<int>("C4: Number of coarsening zones"))-1,0.0);
    Coarsening_List->set("C5: Rmax_zone", Rmax_zone_Array, "define the maximum distance from surface in each zone. \n The distances should be arranged from nearest to furthest where\n the nearest zone to the surface is the most refined zone,\n and the furthest zone from the surfaces is the most coarse zone.");

    

    NonlinearSolver_List->set("NLS1: Nonlinear Solver", "Newton Built-In", "Select nonlinear solver type you would like to use.", NLSolverTypeValidator);
    NonlinearSolver_List->set("NLS2: Max Number of Nonlinear Iterations", 10, "Set the Maximum number of nonlinear iterations", Niter_Validator);
    NonlinearSolver_List->set("NLS3: Logical_Physics_Scaling", true, "Set to true to turn on physics scaling.\n Physics scaling attempts to moderate the effect of terms like exp(mu) in the JDC polymer options");
    NonlinearSolver_List->set("NLS4: Nonlinear Tolerance: Absolute", 0.00000001, "Set absolute convergence tolerance for nonlinear solver",NLtolValidator);
    NonlinearSolver_List->set("NLS5: Nonlinear Tolerance: Relative", 0.0001, "Set relative convergence tolerance for nonlinear solver",NLtolValidator);
    NonlinearSolver_List->set("NLS6: Minimum update step", 1.0, "Set the minimum update fraction if using the built in Newton solver.\n Or set the Mixing fraction if using a Picard Solver.\n If using mixed Picard-Newton Solver enter the minimum Newton update fraction.\n The Picard mixing parameter will be set to the newton update parameter divided by 100.", NLupdate_Validator);

    LinearSolver_List->set("LS1: Linear Solver", "Schur Solver", "Select linear solver type you would like to use.",LSolverTypeValidator);
    LinearSolver_List->set("LS2: Max Number of Linear Iterations", 100, "Set the Maximum number of nonlinear iterations", Niter_Validator);
    LinearSolver_List->set("LS3: Convergence Tolerance for linear solver", 0.001, "Set convergence tolerance for linear solver");
    LinearSolver_List->set("LS4: Linear Solver Approach", "GMRES", "Select linear solver approach you would like to use.",LSolverApproachValidator);
    LinearSolver_List->set("LS5: Matrix Scaling option", "none", "Select linear solver approach you would like to use.",LScalingValidator);
    LinearSolver_List->set("LS6: Preconditioner option", "none", "Select preconditioner approach you would like to use.",LPreconditionerValidator);
    LinearSolver_List->set("LS7: Number of Levels for ILUT", 4, "Set the levels for ilut preconditioner", NlevelILUT_Validator);

    /*******************************/
    /* define dependent parameters */
    /*******************************/

    /**********************************************************************************************/
    /* show the dependent parameters only if the independent parameters has a particular setting. */
    /**********************************************************************************************/
      RCP<StringVisualDependency> PhysScale_Dep = rcp(
           new StringVisualDependency( "F4_POLYMER_Functional",Functional_List,"NLS3: Logical_Physics_Scaling", NonlinearSolver_List, 
               tuple<std::string>("Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

      RCP<StringVisualDependency> LevelILUT_Dep = rcp(
           new StringVisualDependency( "LS6: Preconditioner option",LinearSolver_List,"LS7: Number of Levels for ILUT", LinearSolver_List, 
               tuple<std::string>("ilut")));

      RCP<BoolVisualDependency> JacCoarseType_Dep = rcp(
          new BoolVisualDependency( "C2: jacobian coarsening?",Coarsening_List,"C3: Type of Jacobian coarsening", Coarsening_List));

      RCP<NumberArrayLengthDependency> RmaxZoneLength_Dep = rcp(
          new NumberArrayLengthDependency( "C4: Number of coarsening zones", Coarsening_List, "C5: Rmax_zone",Coarsening_List));

      RCP<NumberVisualDependency<int> > RmaxZoneView_Dep = rcp(
           new NumberVisualDependency<int>("C4: Number of coarsening zones",Coarsening_List,"C5: Rmax_zone", Coarsening_List)); 

    /*****************************************/
    /* add the dependencies for this section.*/
    /*****************************************/
      depSheet_Tramonto->addDependency(PhysScale_Dep);
      depSheet_Tramonto->addDependency(LevelILUT_Dep);
      depSheet_Tramonto->addDependency(JacCoarseType_Dep);
      depSheet_Tramonto->addDependency(RmaxZoneLength_Dep);
      depSheet_Tramonto->addDependency(RmaxZoneView_Dep);
  /****************************************************************************************************************/
  /****************************** END FUNCTIONAL CONTROL PARAMETER SECTION ****************************************/
  /****************************************************************************************************************/
 
  return;
}

