/*
  Prototypes from dft_ldbal.c
*/
struct rcb_dot;
struct rcb_box;
struct rcb_tree;

void rcb (struct rcb_dot **dots,	/* array of dots, may be extended */
	  int *pdotnum,		/* # of dots - decomposition changes it */
	  int *pdotmax,		/* max # of dots array can hold, may extend */
	  int wtflag,		/* (0) do not (1) do use weights */
	  double overalloc,	/* amount to overallocate by when realloc
				   of dot array must be done
				   1.0 = no extra, 1.5 = 50% extra space */
	  int reuse,		/* (0) don't use (1) use previous cuts
				   stored in treept as initial guesses */
	  int *pdottop,		/* dots >= this index are new */
	  struct rcb_box *rcbbox,	/* bounding box of final RCB sub-domain */
	  struct rcb_tree *treept	/* tree of RCB cuts - only single cut on exit */
  );

void rcb_error (int size);
void rcb_box_merge (void *in, void *inout, int *len, MPI_Datatype * dptr);
void rcb_median_merge (void *in, void *inout, int *len, MPI_Datatype * dptr);
void rcb_check (struct rcb_dot *dotpt, int dotnum, int dotorig,
		struct rcb_box *rcbbox);
void rcb_stats (double timetotal, struct rcb_dot *dotpt, int dotnum,
		double *timers, int *counters, struct rcb_box *rcbbox,
		int reuse);
