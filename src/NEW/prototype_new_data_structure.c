
/*
Pseudo code for Tramonto interacting with modified dft_Solvermanager
Prototype code for [node][unknown#] data structures in Tramonto

NOTE: The term "Unknown" now means exclusively the undiscretized
quantity (e.g. density, rho_bar_0) described by the physics being
solved and not the discretized number of unknowns. This term did 
have two meanings, but in the data structure I am proposing here
there is no need to name the discretized degrees-of-freedom.
*/

#define MAX_UNKNOWNS 100
///Tramonto block physics variables:
int Num_Blocks;
int Num_block_unk[Num_Blocks];  //Number of unknowns per physics block
                                //Is currently called Num_total_eq in dft_main
int Num_Unknowns; //sum of previous array over all blocks
int UnknownIndices[Num_Blocks][MAX_UNKNOWNS]; //Indices of unknwons
int Phys[NUM_BLOCKS]; //Identifier of Physics type, e.g DENSITY, POISSON, RHO_BARS,...

// Give the Block Physics Variable values for 2x2 block problem with 6 rho-bars and 1 rho
Num_Blocks = 2;
Num_block_unk[0] = 6;
Num_block_unk[1] = 1;
Num_Unknowns=7;

UnknownIndicies[0][0] = 0;
UnknownIndicies[0][1] = 1;
UnknownIndicies[0][2] = 2;
UnknownIndicies[0][3] = 3;
UnknownIndicies[0][4] = 4;
UnknownIndicies[0][5] = 5;

UnknownIndicies[1][0] = 6;

Phys[0]=PHO_BARS;
Phys[1]=DENSITY;

//Tramonto Mesh Variables:
int num_owned_nodes, num_box_nodes;
  // Maps: Owned-to-global, box-to-global, local-to-box
int OwnedNodes[num_owned_nodes], BoxNodes[num_box_nodes], L2B[num_owned_nodes];

//Construct dft_Solvermanager: DIFFERENT ARGUMENTS THAN CURRENT CODE!
//Physics block variables now passed in here, for consistancy
dft_solvermanager_create(Num_Blocks, Num_block_unk, UnknownIndices, MPI_COMM_WORLD);

//Added "Nodal" to the name to highlight the change from previous version
dft_solvermanager_setNodalrowmap(num_owned_nodes, OwnedNodes);
dft_solvermanager_setNodalcolumnmap(num_box_nodes, BoxNodes);

// The following routine creates the real Maps for the matrix blocks
// but Tramonto doesn't need to know them. The choice of node ordering
// or physics Unknown numbering can be made internally here
dft_solvermanager_finalizeblockstructure();


// Note: No need for graph methods in dft_SolutionManager, we
// would love to fill the matrix directly wiithout a graph-building
// stage.

// Pseudo code for matrix and residual fill:
dft_solvermanager_initializeproblemvalues();

int i,j; // Physics block indices
for (i=0; i< Num_Blocks; i++) {
for (j=0; j< Num_Blocks; j++) {

 //POISSON example implemented...
 if (Phys[i]==POISSON && Phys[j]==POISSON) {
  for (inode=0; inod<num_owned_nodes; inode++) {
  ibox = L2B(inode);

  // residual vector f aka Rhs is accessed with local nodes
  // solution vector x aka Lhs accessed with box nodes
  ipsi = UnknownIndices[i][0]; //Unknown number for Psi variable

  if (FILL_RESIDUAL) {
    // Residual fill for 1D diffusion operator f[inode][ipsi]
    f = x[ibox-1][ipsi] -2*x[ibox][ipsi] + x[ibox+1][ipsi];

    //Row number is now expanded to 2 integers [node][unk#]
    //dft_solutionmanager_insertRhsValue(i, nodeNumber, rowUnkNumber, value);
    dft_solutionmanager_insertRhsValue(i, inode, ipsi, -f);
  }

  if (FILL_JACOBIAN) {
    //Jacobian fill for same: this loads three entries in 1 call
    //nodeIndices is an integer temp array, Values a double temp array
    numEntries=3;
    nodeIndicies[0]=ibox-1; nodeIndicies[1]=ibox; nodeIndicies[2]=ibox+1; 
    Values[0]=1.0; Values[1]=-2.0; Values[2]=1.0; 

    // Row and Column are now each expressed with 2 integers: node,unknown
    // so an extra 2 arguments added to this routine: rowUnkNumber and colUnkNumber
    // I think the colUnkNumber can be an integer and not a more general array of 
    // integers.
    //dft_solutionmanager_insertMatrixValues(i, j, nodeNumber, rowUnkNumber,
    //                                       numEntries, Values, nodeIndices, colUnkNumber);
    dft_solutionmanager_insertMatrixValues(i, j, inode, ipsi, numEntries, Values, nodeIndices, ipsi);
  }

  }//end inode loop
 }//end IF POISSON 

 // Some code fragments for rho_bar_0 fill
 // Looking at A12 block, derivative of rho_bar_0 w.r.t density
 
 if (Phys[i]==RHO_BAR && Phys[j]==DENSITY) {
 irb0 =  UnknownIndices[i][0]; //Unknown number for rho_bar_0
 jden =  UnknownIndices[j][0]; //Unknown number for density
   for (inode=0; inod<num_owned_nodes; inode++) {
     ibox = L2B(inode);
     f=0;

     // Uses RB0 structure holding Stencil info for this integral
     for (isten=0; isten<RB0.stenPoints; isten++) {
       jbox = RB0.StencilOffset(isten, ibox);
       weight = RB0.StencilWeight(isten);
       f -= weight*x[jbox][jden];

       // Add 1 entry into Jacobian
       Values[0]=-weight; nodeIndices[0] = jbox;
       dft_solutionmanager_insertMatrixValues(i, j, inode, irb0, 1, Values, &jbox, iden);
     }
     // LOAD Residual contribution from this physics block,node,unknown
     dft_solutionmanager_insertRhsValue(i, inode, irb0, -f);
   }
 }
 
}//end j 
}//end i

// Mike: It would be great to have an internal flag e.g. "MatrixIsFinalized" 
// in dft_SolutionManager that allows us to fill the matrix graph/values
// from scratch the first fill, but then keeps the graph around after that.
// The first time this finalizeproblemvalues method gets called, OptimizeStorage
// gets called and then this new flag gets switched to true. I am not sure
// whether the  insertMatrixValues  routine will be different based on
// this flag (Insert for first fill vs SumInto for subsequent). A resetMatrix
// method can be added in case we want to change the sparsity pattern during 
// a run.
dft_solvermanager_finalizeproblemvalues();



//Solve/Newton
dft_solvermanager_solve();

// Newton update: Lhs is 2D array
dft_solvermanager_getLhs(delta_x);
for (i=0; i<num_box_nodes; i++)
  for (j=0; j<Num_Unknowns; j++)
    x[i][j] += delta_x[i][j];

printf("Total Time = %g (sec)\n", old_time*0.1);
