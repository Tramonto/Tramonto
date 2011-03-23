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

    RCP<StringValidator> CoarsenTypeValidator = rcp(
           new StringValidator(tuple<std::string>("none","Bulk zone","Poisson-Boltzmann zone", "1D zone in 2D or 3D problem",
                                  "residual and jacobian coarsening","jacobian only coarsening")));

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
    RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));

    RCP<StringValidator> JacCoarsenTypeValidator = rcp(
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


    Coarsening_List->set("C1: Type of coarsening", "none", "Select a method for coarsening of the problem.\n options: bulk_zone: Set to true to replace DFT Euler-Lagrange equation with rho(r)=rho_b some distance from teh surfaces.\n \t Poisson-Bolzmann zone: Set to true to replace DFT Euler-Lagrange equation with rho(r)=rho_b some distance from teh surfaces.\n \t 1D zone: Set to true if the 2D or 3D problem will converge to a 1D solution at some distance from the surface.\n\t Residual coarsening: Coarsening of residual equations and Jacobian integrals in zones away from the surfaces.\n\t Jacobian coarsening: Coarsening of jacobian integrals only in zones away from surfaces.", CoarsenTypeValidator);

    Coarsening_List->set("C2: truncate integrals in jacobian?" ,false,"Set to true for truncation of jacobian integrals below some threshhold.");
    Coarsening_List->set("C3: Type of Jacobian coarsening","none","Select a method for coarsening the mesh or jacobian",JacCoarsenTypeValidator);
    Coarsening_List->set("C4: Number of coarsened zones",0,"Set number of coarsened zones in problem (don't count most refined zone).");

    Array<double> Rmax_zone_Array( (Coarsening_List->get<int>("C4: Number of coarsened zones")),0.0);
    Coarsening_List->set("C5: Rmax_zone", Rmax_zone_Array, "define the maximum distance from surface in each zone. \n The distances should be arranged from nearest to furthest where\n the nearest zone to the surface is the most refined zone,\n and the furthest zone from the surfaces is the most coarse zone.");
    Coarsening_List->set("C6: Dimension 1D_BC",0,"This is direction that still has density variations in the 1D zone (0=x,1=y,2=z).",DimValidator);
    Coarsening_List->set("C7: X 1D_BC",0.0,"Distance from domain boundary where 1D zone should be applied. Note that the boundary is applied on both sides of the domain.");

    

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
         new StringVisualDependency( 
             Functional_List->getEntryRCP("F4_POLYMER_Functional"),
             NonlinearSolver_List->getEntryRCP("NLS3: Logical_Physics_Scaling"),
             tuple<std::string>("Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)",
                                "Polymer_JDC_iSAFT(comp)"))
     );

     RCP<StringVisualDependency> LevelILUT_Dep = rcp(
         new StringVisualDependency( 
             LinearSolver_List->getEntryRCP("LS6: Preconditioner option"),
             LinearSolver_List->getEntryRCP("LS7: Number of Levels for ILUT"), 
             tuple<std::string>("ilut"))
     );

     RCP<StringVisualDependency> TypeJacCoarse_Dep = rcp(
         new StringVisualDependency( 
             Coarsening_List->getEntryRCP("C1: Type of coarsening"),
             Coarsening_List->getEntryRCP("C3: Type of Jacobian coarsening"), 
             tuple<std::string>("residual and jacobian coarsening","jacobian only coarsening"))
     );

     RCP<StringVisualDependency> Nzone_Dep = rcp(
         new StringVisualDependency( 
             Coarsening_List->getEntryRCP("C1: Type of coarsening"),
             Coarsening_List->getEntryRCP("C4: Number of coarsened zones"), 
             tuple<std::string>("residual and jacobian coarsening",
                  "jacobian only coarsening","Bulk zone","Poisson-Boltzmann zone",
                  "1D zone in 2D or 3D problem"))
     );

     RCP<StringVisualDependency> RmaxZone_Dep = rcp(
         new StringVisualDependency( 
             Coarsening_List->getEntryRCP("C1: Type of coarsening"),
             Coarsening_List->getEntryRCP("C5: Rmax_zone"),  
             tuple<std::string>("residual and jacobian coarsening",
              "jacobian only coarsening","Bulk zone","Poisson-Boltzmann zone"))
     );

     RCP<StringVisualDependency> Dim1DZone_Dep = rcp(
         new StringVisualDependency( 
             Coarsening_List->getEntryRCP("C1: Type of coarsening"),
             Coarsening_List->getEntryRCP("C6: Dimension 1D_BC"),  
             tuple<std::string>("1D zone in 2D or 3D problem"))
     );

     RCP<StringVisualDependency> X1DZone_Dep = rcp(
         new StringVisualDependency( 
             Coarsening_List->getEntryRCP("C1: Type of coarsening"),
             Coarsening_List->getEntryRCP("C7: X 1D_BC"),
             tuple<std::string>("1D zone in 2D or 3D problem"))
     );

     RCP<NumberArrayLengthDependency<int, double> > RmaxZoneLength_Dep = rcp(
         new NumberArrayLengthDependency<int, double> ( 
             Coarsening_List->getEntryRCP("C4: Number of coarsened zones"), 
             Coarsening_List->getEntryRCP("C5: Rmax_zone"))
     );


    /*****************************************/
    /* add the dependencies for this section.*/
    /*****************************************/
      depSheet_Tramonto->addDependency(PhysScale_Dep);
      depSheet_Tramonto->addDependency(LevelILUT_Dep);

      depSheet_Tramonto->addDependency(TypeJacCoarse_Dep);
      depSheet_Tramonto->addDependency(Nzone_Dep);
      depSheet_Tramonto->addDependency(RmaxZone_Dep);
      depSheet_Tramonto->addDependency(Dim1DZone_Dep);
      depSheet_Tramonto->addDependency(X1DZone_Dep);
      depSheet_Tramonto->addDependency(RmaxZoneLength_Dep);
  /****************************************************************************************************************/
  /****************************** END FUNCTIONAL CONTROL PARAMETER SECTION ****************************************/
  /****************************************************************************************************************/
 
  return;
}

