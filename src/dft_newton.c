/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
 
/*
 *  FILE: dft_newton.c
 *
 *  This file solves the dft problem using Newton's method.
 */
#include "mpi.h" 


#ifdef HAVE_DFT_SCHUR_SOLVER
#include "dft_schur_solver.h"
#endif

#include "dft_globals_const.h" 
#include "rf_allo.h"

static void print_resid_norm(int iter);
void fill_test(double **x, int flag);

/* NEWTON SOLVER using dft_SolverManager */
int solve_problem(double **x, double **x2)
/*
 * x and x2 in dft_main are allocated [Nunk_per_node][Nnodes_box]
 * x2 only relevent for Lbinodal calculations, so can mostly be ignored
 */

{
  int iter;
  double **xOwned;

  /* Construct dft_SolverManager with information on number of unknowns*/
  Solver_manager = dft_solvermanager_create(Nunk_per_node, Unk2Phys, NULL, NULL, MPI_COMM_WORLD);

  /* Give Nodal Row and Column maps */
  (void) dft_solvermanager_setnodalrowmap(Solver_manager, Nnodes_per_proc, L2G_node);
  (void) dft_solvermanager_setnodalcolmap(Solver_manager, Nnodes_box     , B2G_node);

  /* SolverManager can now set up its own numbering scheme, set up unknown-based Maps */
  (void) dft_solvermanager_finalizeblockstructure(Solver_manager);

  /* Set initial guess on owned nodes and reconcile ghost nodes using importr2c */
  xOwned = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
  set_initial_guess(Iguess1, xOwned);
  (void) dft_solvermanager_importr2c(Solver_manager, xOwned, x);

  /* If requested, write out initial guess */
   if (Iwrite == VERBOSE) print_profile_box(x,"rho_init.dat");

  /* Do same for second solution vector when Lbinodal is true */
  if (Lbinodal) {
    set_initial_guess(BINODAL_FLAG, xOwned);
    (void) dft_solvermanager_importr2c(Solver_manager, xOwned, x2);
    if (Iwrite == VERBOSE) print_profile_box(x2,"rho_init2.dat");
  }
  safe_free((void **) &xOwned);

#ifdef LOCA_REIMPLEMENTED
  if (Loca.method != -1)
    iter = solve_continuation(x, x2, Sten_Type[POLYMER_CR], aux_info);
  else
#endif
     iter = newton_solver(x, NULL);

  /* Call the destructor for the dft_SolverManager */
  dft_solvermanager_destruct(Solver_manager);

  return(iter);
}

/****************************************************************************/
/****************************************************************************/
void box2owned(double** xBox, double** xOwned) {
  /* Routine to restrict 2d box array to 2d node array */
  /* This is the opposite of importr2c, and strips ghost nodes */
  int iunk, inode, ibox;

  for (ibox=0; ibox<Nnodes_box; ibox++) {
    inode = B2L_node[ibox];
    if (inode!=-1) {
      for (iunk=0; iunk<Nunk_per_node; iunk++)
        xOwned[iunk][inode] = xBox[iunk][ibox];
    }
  }
}

/****************************************************************************/
/****************************************************************************/
int newton_solver(double** x, void* con_ptr) {

  int iter=0;
  int converged=FALSE, converged2=TRUE;
  double** delta_x;
  delta_x = (double **) array_alloc(2, Nunk_per_node, Nnodes_box, sizeof(double));

  do {
    iter++;

    (void) dft_solvermanager_initializeproblemvalues(Solver_manager);

    /* Call Matrix and Residual Fill routine, resid_only_flag=FALSE)*/
    /*fill_resid_and_matrix_control(x, FALSE); */
    fill_test(x, FALSE);

    (void) dft_solvermanager_finalizeproblemvalues(Solver_manager);
    if (Iwrite != NO_SCREEN) print_resid_norm(iter);
    (void) dft_solvermanager_setupsolver(Solver_manager);
    (void) dft_solvermanager_solve(Solver_manager);
    
    /* I am assuming getLhs returns box coordinates (e.g. Column Map)!! */
    (void) dft_solvermanager_getlhs(Solver_manager, delta_x);
    
    /* Do: x += delta_x, and check for convergence */
    converged = update_solution(x, delta_x, iter);

  } while (iter < Max_Newton_iter && (!converged || !converged2));

  if (!converged || !converged2) {
    printf("\tNewton Solver: Failed to converge in %d iterations\n",iter);
    iter = -iter;
  }
  else
    printf("\tNewton Solver: Successful convergence in %d iterations\n",iter);

  safe_free((void *) &delta_x);
  return iter;
}

/****************************************************************************/
/****************************************************************************/
int update_solution(double** x, double** delta_x, int iter) {
/* Routine to update solution vector x using delta_x and 
 * to test for convergence of Newton's method.
 * 
 * Note: Modifications to Newton's method, including damping, have not yet
 *       been translated form previous update_solution.
 *       iter  value may be used in some damping methods
 */

  int iunk, ibox, inode;
  double updateNorm=0.0, temp;
  char *yo = "newupdate solution";

  for (ibox=0; ibox<Nnodes_box; ibox++) {

    /* Increment updateNorm only for owned nodes (inode=-1 for ghosts) */
    inode = B2L_node[ibox];
    if (inode != -1) {
      for (iunk=0; iunk<Nunk_per_node; iunk++) {
        temp = delta_x[iunk][ibox]/(Newton_rel_tol*x[iunk][ibox] + Newton_abs_tol);
        updateNorm +=  temp*temp;
      }
    }

    /* Update all solution componenets */
    for (iunk=0; iunk<Nunk_per_node; iunk++)
      x[iunk][ibox] += delta_x[iunk][ibox];
  }

  updateNorm = sqrt(gsum_double(updateNorm));

  if (Proc==0 && Iwrite != NO_SCREEN)
    printf("\n\t\t%s: Weighted norm of update vector =  %g\n", yo, updateNorm);

  if (updateNorm > 1.0) return(FALSE);
  else                  return(TRUE);

}

/*****************************************************************************************************/
static void print_resid_norm(int iter)
{
  int iunk, j;
  double norm=0.0;
  double **f;
  f = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
  dft_solvermanager_getrhs(Solver_manager, f);

  for (iunk=0; iunk<Nunk_per_node; iunk++) {
    for (j=0; j< Nnodes_per_proc; j++) {
       norm += f[iunk][j] * f[iunk][j];
    }
  }
  norm = gsum_int(norm);
  if (Proc==0) printf("\t\tResidual norm at iteration %d = %g\n",iter, sqrt(norm));
}
/*****************************************************************************************************/

void fill_test(double **x, int flag) 
/* Quick test problem for matrix loading, solve, and Newton's method */
/* For the first unknown at a node, x[0]=inode 
 * for subsequent unknowns at that node, x[iunk]^2=x[iunk-1]
 * So, when run at 3 unknowns per node, global node 16 should have
 * solutions 16, +-4, +-2, for the three unknowns at that node
 */
{
  int loc_inode, inode_box, iunk, inode;
  double f;

  for (loc_inode=0; loc_inode< Nnodes_per_proc; loc_inode++) {
    inode = L2G_node[loc_inode];
    inode_box = L2B_node[loc_inode];
    iunk = 0;

    f = x[iunk][inode_box] - inode - 1;
    dft_solvermanager_insertrhsvalue(Solver_manager, iunk, loc_inode, -f);
    dft_solvermanager_insertonematrixvalue(Solver_manager,iunk,loc_inode,iunk,inode_box,1.0);

    /* For iunk 1 and higher, the square of the unknow is the previous value */
    for (iunk=1; iunk<Nunk_per_node; iunk++) {
      f = x[iunk][inode_box]*x[iunk][inode_box] - x[iunk-1][inode_box];
      dft_solvermanager_insertrhsvalue(Solver_manager, iunk, loc_inode, -f);
      dft_solvermanager_insertonematrixvalue(Solver_manager,iunk,loc_inode,iunk,inode_box, 2*x[iunk][inode_box]);
      dft_solvermanager_insertonematrixvalue(Solver_manager,iunk,loc_inode,iunk-1,inode_box, -1.0);
    }
  } 
}

