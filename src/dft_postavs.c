#include <stdio.h>

#define MAXNODES 500000
#define MAXCELLS 125000
#define MAXVARPERNODE 3
#define MAXVARPERCELL 3

int NumProc;
int Num_Node_Values;
int Num_Cell_Values;
int Total_AVS_cells = 0;
int Total_AVS_nodes = 0;
int AVS_Call;

int CellNodes[MAXCELLS][4];
double NodeX[MAXNODES], NodeY[MAXNODES];
double NodeValues[MAXNODES][MAXVARPERNODE];
double CellValues[MAXCELLS][MAXVARPERCELL];

FILE *OutFP;

/************************************************************************/
int
find_node (double x, double y)
{
  int i, node_num = -1;

  for (i = 0; i < Total_AVS_nodes; i++)
    {
      if ((NodeX[i] == x) && (NodeY[i] == y))
	{
	  node_num = i;
	  break;
	}
    }

  return (node_num);
}

/************************************************************************/
int
find_cell (double x, double y)
{
  int i, cell_num = -1;
  int node;

  for (i = 0; i < Total_AVS_cells; i++)
    {
      if ((x >= NodeX[CellNodes[i][0]]) && (x <= NodeX[CellNodes[i][2]]))
	{
	  if ((y >= NodeY[CellNodes[i][1]]) && (y <= NodeY[CellNodes[i][3]]))
	    {
	      cell_num = i;
	      break;
	    }
	}
    }

  return (cell_num);
}

/************************************************************************/
void
read_node_coords ()
{
  FILE *fp;
  char filename[16];
  int proc, i, j;
  double x, y, junk;

  for (proc = 0; proc < NumProc; proc++)
    {
      sprintf (filename, "avs_mesh%03d%04d", AVS_Call, proc);
      fp = fopen (filename, "r");
      if (fp == NULL)
	{
	  printf ("Error in read_node_coords:  file %s not found\n",
		  filename);
	  exit (0);
	}

      i = 0;
      while (fscanf (fp, "%le %le", &x, &y) != EOF)
	{

	  if ((CellNodes[Total_AVS_cells][i] = find_node (x, y)) == -1)
	    {
	      NodeX[Total_AVS_nodes] = x;
	      NodeY[Total_AVS_nodes] = y;
	      for (j = 0; j < Num_Node_Values; j++)
		fscanf (fp, "%le ", &(NodeValues[Total_AVS_nodes][j]));
	      CellNodes[Total_AVS_cells][i] = Total_AVS_nodes;
	      Total_AVS_nodes++;
	    }
	  else
	    {
	      /* skip over the data values */
	      for (j = 0; j < Num_Node_Values; j++)
		fscanf (fp, "%le ", &junk);
	    }

	  i++;
	  if (i == 4)
	    {
	      Total_AVS_cells++;
	      i = 0;
	    }
	}
      fclose (fp);
    }
}


/************************************************************************/
void
read_cell_values ()
{
  FILE *fp;
  char filename[16];
  int proc, cell, j;
  int elid, norder, limit;
  double x, y;
  double rho, u, v, e, pres, errest;

  for (proc = 0; proc < NumProc; proc++)
    {
      sprintf (filename, "avs_cells%03d%04d", AVS_Call, proc);
      fp = fopen (filename, "r");
      if (fp == NULL)
	{
	  printf ("Error in read_cell_values:  file %s not found\n",
		  filename);
	  exit (0);
	}

      while (fscanf (fp, "%le %le ", &x, &y) != EOF)
	{
	  if ((cell = find_cell (x, y)) == -1)
	    printf ("ERROR Cell not found: (%e,%e)\n", x, y);
	  else
	    {
	      for (j = 0; j < Num_Cell_Values; j++)
		fscanf (fp, "%le ", &(CellValues[cell][j]));
	    }
	}
      fclose (fp);
    }
}

/************************************************************************/
void
print_node_coords ()
{
  int i;

  for (i = 0; i < Total_AVS_nodes; i++)
    {
      fprintf (OutFP, "%d %e %e %e\n", i, NodeX[i], NodeY[i], 0.0);
    }
}

/************************************************************************/
void
print_cell_nodes ()
{
  int i, j;

  for (i = 0; i < Total_AVS_cells; i++)
    {
      fprintf (OutFP, "%d 1 quad ", i);
      for (j = 0; j < 4; j++)
	{
	  fprintf (OutFP, "%d ", CellNodes[i][j]);
	}
      fprintf (OutFP, "\n");
    }
}

/************************************************************************/
void
print_node_data ()
{
  int i, j;

  for (i = 0; i < Total_AVS_nodes; i++)
    {
      fprintf (OutFP, "%d ", i);
      for (j = 0; j < Num_Node_Values; j++)
	fprintf (OutFP, "%le ", NodeValues[i][j]);
      fprintf (OutFP, "\n");
    }
}

/************************************************************************/
void
print_cell_data ()
{
  int i, j;

  for (i = 0; i < Total_AVS_cells; i++)
    {
      fprintf (OutFP, "%d ", i);
      for (j = 0; j < Num_Cell_Values; j++)
	fprintf (OutFP, "%le ", CellValues[i][j]);
      fprintf (OutFP, "\n");

    }
}


/************************************************************************/
void
main (int argc, char *argv[])
{
  char filename[16];
  int j;

  NumProc = atoi (argv[1]);

  AVS_Call = atoi (argv[2]);

  Num_Node_Values = atoi (argv[3]);

  Num_Cell_Values = atoi (argv[4]);

  read_node_coords ();

  read_cell_values ();

  sprintf (filename, "avs%02d.inp", AVS_Call);
  OutFP = fopen (filename, "w");

  fprintf (OutFP, "%d %d %d %d %d\n",
	   Total_AVS_nodes, Total_AVS_cells, Num_Node_Values,
	   Num_Cell_Values, 0);

  print_node_coords ();

  print_cell_nodes ();

  /* components of nodal solutions, but AVS wants one. */
  fprintf (OutFP, "%d ", Num_Node_Values);
  for (j = 0; j < Num_Node_Values; j++)
    fprintf (OutFP, "1 ");
  fprintf (OutFP, "\n");
  for (j = 0; j < Num_Node_Values; j++)
    fprintf (OutFP, "Y%d, no_units \n", j);


  print_node_data ();

  /* components of cell solutions  */
  fprintf (OutFP, "%d ", Num_Cell_Values);
  for (j = 0; j < Num_Cell_Values; j++)
    fprintf (OutFP, "1 ");
  fprintf (OutFP, "\n");
  for (j = 0; j < Num_Cell_Values; j++)
    fprintf (OutFP, "u%d, no_units\n", j);

  print_cell_data ();

  fclose (OutFP);
}
