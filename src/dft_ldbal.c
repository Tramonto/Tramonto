/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/* Recursive coordinate bisectioning (RCB) routine
   operates on "dots" as defined in rcb.h */

/* Notes:
   dots are balanced across procs by weight (if used)
   on return, proc owns dotnum "dots" in dense array of max-length dotmax
   all dots will be inside (or on surface of) 3-d box defined by rcbbox
   input weights (if used) are real numbers > 0.0
   can extend "rcb_dot" data structure in calling program, see rcb.h
   returned RCB tree only contains one cut on each proc,
     need to do MPI_Allgather if wish to collect it on all procs
   set RCB_CHECK and RCB_STATS at top of rcb.c as desired
*/

#include "dft_ldbal.h"

#define MYHUGE 1.0e30
#define TINY   1.0e-6

/*  RCB_CHECK = 0  No consistency check on input or results */
/*  RCB_CHECK = 1  Check input weights and final results for consistency */
int RCB_CHECK = 1;

/*  RCB_STATS = 0  No statistics logging */
/*  RCB_STATS = 1  Log times and counts, print summary */
/*  RCB_STATS = 2  Log times and counts, print for each proc */
int RCB_STATS = 1;

void rcb(

struct rcb_dot  **dots,       /* array of dots, may be extended */
int              *pdotnum,    /* # of dots - decomposition changes it */
int              *pdotmax,    /* max # of dots array can hold, may extend */
int               wtflag,     /* (0) do not (1) do use weights */
double            overalloc,  /* amount to overallocate by when realloc
                                 of dot array must be done
                                 1.0 = no extra, 1.5 = 50% extra space */
int               reuse,      /* (0) don't use (1) use previous cuts
                                  stored in treept as initial guesses */
int              *pdottop,    /* dots >= this index are new */
struct rcb_box   *rcbbox,     /* bounding box of final RCB sub-domain */
struct rcb_tree  *treept      /* tree of RCB cuts - only single cut on exit */

)

{
  int    proc,nprocs;              /* my proc id, total # of procs */
  struct rcb_dot *dotbuf, *dotpt;  /* local dot arrays */
  struct rcb_median med, medme;    /* median data */
  struct rcb_box boxtmp;           /* tmp rcb box */
  int    keep, outgoing;           /* message exchange counters */
  int    incoming, incoming2;      /* message exchange counters */
  int   *dotmark;                  /* which side of median for each dot */
  int   *dotlist;                  /* list of active dots */
  int    dotnum;                   /* number of dots */
  int    dotmax;                   /* max # of dots arrays can hold */
  int    dottop;                   /* dots >= this index are new */
  int    dotnew;                   /* # of new dots after send/recv */
  int    numlist;                  /* number of active dots I own */
  int    proclower, procupper;     /* lower/upper proc in partition */
  int    procmid;                  /* 1st proc in upper half of part */
  int    procpartner, procpartner2=0; /* proc(s) to exchange with */
  int    readnumber;               /* # of proc partner(s) to read from */
  int    markactive;               /* which side of cut is active = 0/1 */
  int    dim;                      /* which of 3 axes median cut is on */
  int    indexlo=0, indexhi=0;         /* indices of dot closest to median */
  double valuemin, valuemax;       /* lower/upper bounds of active region */
  double valuehalf;                /* median cut position */
  double targetlo, targethi;       /* desired wt in lower/upper half */
  double weight;                   /* wt of entire partition */
  double weightlo, weighthi;       /* wt in lower/upper half of non-active */
  double wtsum,wtok,wtupto,wtmax;  /* temporary wts */
  double tolerance;                /* largest single weight of a dot */
  int    first_iteration;          /* flag for 1st time thru median search */
  int    allocflag;                /* have to re-allocate space */
  int    breakflag;                /* for breaking out of median iteration */
  int    length;                   /* message length */
  double time1=0,time2=0,time3=0,time4;  /* timers */
  double timestart=0,timestop=0;       /* timers */
  double timers[4];                /* diagnostic timers 
			              0 = start-up time before recursion
				      1 = time before median iterations
				      2 = time in median iterations
				      3 = communication time */
  int    counters[7];              /* diagnostic counts
			              0 = # of median iterations
				      1 = # of dots sent
				      2 = # of dots received
				      3 = most dots this proc ever owns
				      4 = most dot memory this proc ever allocs
				      5 = # of times a previous cut is re-used
				      6 = # of reallocs of dot array */
  int    i,j,k;                    /* local variables */

  /* function prototypes */

  void rcb_error(int);
  void rcb_check(struct rcb_dot *, int, int, struct rcb_box *);
  void rcb_stats(double, struct rcb_dot *,int, double *, 
		 int *, struct rcb_box *, int);

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;
  MPI_Request request, request2;
  MPI_Status status;
  MPI_Op box_op, med_op;
  MPI_Datatype box_type, med_type;
  MPI_User_function rcb_box_merge, rcb_median_merge;

  if (RCB_STATS) {
    MPI_Barrier(MPI_COMM_WORLD);
    timestart = time1 = MPI_Wtime();
  }
  /* setup for parallel */

  MPI_Comm_rank(MPI_COMM_WORLD,&proc);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  /* local copies of calling parameters */

  dotpt = *dots;
  dottop = dotnum = *pdotnum;
  dotmax = *pdotmax;

  /* initialize timers and counters */

  if (RCB_STATS) {
    counters[0] = 0;
    counters[1] = 0;
    counters[2] = 0;
    counters[3] = dotnum;
    counters[4] = dotmax;
    counters[5] = 0;
    counters[6] = 0;
    timers[1] = 0.0;
    timers[2] = 0.0;
    timers[3] = 0.0;
  }

  /* create mark and list arrays for dots */

  allocflag = 0;
  if (dotmax > 0) {
    dotmark = (int *) malloc((unsigned) dotmax * sizeof(int));
    dotlist = (int *) malloc((unsigned) dotmax * sizeof(int));
    if (dotmark == NULL || dotlist == NULL) rcb_error(dotmax*sizeof(int));
  }
  else {
    dotmark = NULL;
    dotlist = NULL;
  }

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(6,MPI_DOUBLE,&box_type);
  MPI_Type_commit(&box_type);
  MPI_Type_contiguous(sizeof(struct rcb_median),MPI_CHAR,&med_type);
  MPI_Type_commit(&med_type);

  MPI_Op_create(&rcb_box_merge,1,&box_op);
  MPI_Op_create(&rcb_median_merge,1,&med_op);

  /* set dot weighs = 1.0 if user didn't */

  if (!wtflag)
    for (i = 0; i < dotnum; i++) dotpt[i].weight = 1.0;

  /* check that all weights > 0 */

  if (RCB_CHECK) {
    j = 0;
    for (i = 0; i < dotnum; i++) if (dotpt[i].weight <= 0.0) j++;
    MPI_Allreduce(&j,&k,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (k > 0) {
      if (proc == 0) printf("RCB ERROR: %d dot weights are <= 0\n",k);
      return;
    }
  }

  /* initialize sub-domain bounding box to entire domain */

  boxtmp.lo[0] = boxtmp.lo[1] = boxtmp.lo[2] = MYHUGE;
  boxtmp.hi[0] = boxtmp.hi[1] = boxtmp.hi[2] = -MYHUGE;
  
  for (i = 0; i < dotnum; i++) {
    for (j = 0; j < 3; j++) {
      if (dotpt[i].x[j] < boxtmp.lo[j])
	boxtmp.lo[j] = dotpt[i].x[j];
      if (dotpt[i].x[j] > boxtmp.hi[j])
	boxtmp.hi[j] = dotpt[i].x[j];
    }
  }

  MPI_Allreduce(&boxtmp,rcbbox,1,box_type,box_op,MPI_COMM_WORLD);

  /* create local communicator for use in recursion */

  MPI_Comm_dup(MPI_COMM_WORLD,&local_comm);

  if (RCB_STATS) {
    time2 = MPI_Wtime();
    timers[0] = time2 - time1;
  }

  /* recursively halve until just one proc in partition */
  
  proclower = 0;
  procupper = nprocs - 1;
  

  while (proclower != procupper) {

    if (RCB_STATS) time1 = MPI_Wtime();
    
    /* procmid = 1st proc in upper half of partition */
    /* if odd # of procs, lower partition gets extra one */

    procmid = proclower + (procupper - proclower) / 2 + 1;

    /* determine communication partner(s) */

    if (proc < procmid)
      procpartner = proc + (procmid - proclower);
    else
      procpartner = proc - (procmid - proclower);
    
    readnumber = 1;
    if (procpartner > procupper) {
      readnumber = 0;
      procpartner--;
    }
    if (proc == procupper && procpartner != procmid - 1) {
      readnumber = 2;
      procpartner2 = procpartner + 1;
    }
    
    /* weight = summed weight of entire partition */
    /* search tolerance = largest single weight (plus epsilon) */
    /* targetlo = desired weight in lower half of partition */
    /* targethi = desired weight in upper half of partition */

    wtmax = wtsum = 0.0;
    for (i = 0; i < dotnum; i++) {
      wtsum += dotpt[i].weight;
      if (dotpt[i].weight > wtmax) wtmax = dotpt[i].weight;
    }

    MPI_Allreduce(&wtsum,&weight,1,MPI_DOUBLE,MPI_SUM,local_comm);
    MPI_Allreduce(&wtmax,&tolerance,1,MPI_DOUBLE,MPI_MAX,local_comm);

    tolerance *= 1.0 + TINY;
    targetlo = ((double) (procmid - proclower)) / 
      (procupper + 1 - proclower) * weight;
    targethi = weight - targetlo;

    /* dim = dimension (xyz = 012) to bisect on */

    dim = 0;
    if (rcbbox->hi[1] - rcbbox->lo[1] >
	rcbbox->hi[0] - rcbbox->lo[0])
      dim = 1;
    if (dim == 0 && rcbbox->hi[2] - rcbbox->lo[2] >
	rcbbox->hi[0] - rcbbox->lo[0])
      dim = 2;
    if (dim == 1 && rcbbox->hi[2] - rcbbox->lo[2] >
	rcbbox->hi[1] - rcbbox->lo[1])
      dim = 2;
    
    /* create mark array and active list for dots */

    if (allocflag) {
      allocflag = 0;
      free(dotlist);
      free(dotmark);
      dotmark = (int *) malloc((unsigned) dotmax * sizeof(int));
      dotlist = (int *) malloc((unsigned) dotmax * sizeof(int));
      if (dotmark == NULL || dotlist == NULL) rcb_error(dotmax*sizeof(int));
    }

    /* initialize active list to all dots */

    for (i = 0; i < dotnum; i++) dotlist[i] = i;
    numlist = dotnum;

    /* weightlo/hi = total weight in non-active parts of partition */

    weightlo = weighthi = 0.0;
    valuemin = rcbbox->lo[dim];
    valuemax = rcbbox->hi[dim];
    first_iteration = 1;

    if (RCB_STATS) time2 = MPI_Wtime();

    /* median iteration */
    /* zoom in on bisector until correct # of dots in each half of partition */
    /* as each iteration of median-loop begins, require:
            all non-active dots are marked with 0/1 in dotmark
	    valuemin <= every active dot <= valuemax
            weightlo, weighthi = total wt of non-active dots */
    /* when leave median-loop, require only:
            valuehalf = correct cut position
            all dots <= valuehalf are marked with 0 in dotmark
            all dots >= valuehalf are marked with 1 in dotmark */

    while (1) {

      /* choose bisector value */
      /* use old value on 1st iteration if old cut dimension is the same */
      /* on 2nd option: could push valuehalf towards geometric center 
	 with "1.0-factor" to force overshoot */

      if (first_iteration && reuse && dim == treept[procmid].dim) {
	if (RCB_STATS) counters[5]++;
	valuehalf = treept[procmid].cut;
	if (valuehalf < valuemin || valuehalf > valuemax)
	  valuehalf = 0.5 * (valuemin + valuemax);	  
      }
      else if (weight)
	valuehalf = valuemin + (targetlo - weightlo) /
	  (weight - weightlo - weighthi) * (valuemax - valuemin);
      else
	valuehalf = 0.5 * (valuemin + valuemax);

      first_iteration = 0;
      
      /* initialize local median data structure */

      medme.lototal = medme.hitotal = 0.0;
      medme.valuelo = -MYHUGE;
      medme.valuehi = MYHUGE;
      medme.wtlo = medme.wthi = 0.0;
      medme.countlo = medme.counthi = 0;
      medme.proclo = medme.prochi = proc;

      /* mark all active dots on one side or other of bisector */
      /* also set all fields in median data struct */
      /* save indices of closest dots on either side */

      for (j = 0; j < numlist; j++) {
	i = dotlist[j];
	if (dotpt[i].x[dim] <= valuehalf) {            /* in lower part */
	  medme.lototal += dotpt[i].weight;
	  dotmark[i] = 0;
	  if (dotpt[i].x[dim] > medme.valuelo) {       /* my closest dot */
	    medme.valuelo = dotpt[i].x[dim];
	    medme.wtlo = dotpt[i].weight;
	    medme.countlo = 1;
	    indexlo = i;
	  }                                            /* tied for closest */
	  else if (dotpt[i].x[dim] == medme.valuelo) {
	    medme.wtlo += dotpt[i].weight;
	    medme.countlo++;
	  }
	}
	else {                                         /* in upper part */
	  medme.hitotal += dotpt[i].weight;
	  dotmark[i] = 1;
	  if (dotpt[i].x[dim] < medme.valuehi) {       /* my closest dot */
	    medme.valuehi = dotpt[i].x[dim];
	    medme.wthi = dotpt[i].weight;
	    medme.counthi = 1;
	    indexhi = i;
	  }                                            /* tied for closest */
	  else if (dotpt[i].x[dim] == medme.valuehi) {
	    medme.wthi += dotpt[i].weight;
	    medme.counthi++;
	  }
	}
      }

      /* combine median data struct across current subset of procs */

      if (RCB_STATS) counters[0]++;
      MPI_Allreduce(&medme,&med,1,med_type,med_op,local_comm);

      /* test median guess for convergence */
      /* move additional dots that are next to cut across it */
                                              
      if (weightlo + med.lototal < targetlo) {    /* lower half TOO SMALL */

	weightlo += med.lototal;
	valuehalf = med.valuehi;

	if (med.counthi == 1) {                  /* only one dot to move */
	  if (weightlo + med.wthi < targetlo) {  /* move it, keep iterating */
	    if (proc == med.prochi) dotmark[indexhi] = 0;
	  }
	  else {                                 /* only move if beneficial */
	    if (weightlo + med.wthi - targetlo < targetlo - weightlo)
	      if (proc == med.prochi) dotmark[indexhi] = 0;
	    break;                               /* all done */
	  }
	}
	else {                                   /* multiple dots to move */
	  breakflag = 0;
	  wtok = 0.0;
	  if (medme.valuehi == med.valuehi) wtok = medme.wthi;   
	  if (weightlo + med.wthi >= targetlo) { /* can't move all of them */
	    MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,local_comm);
	    wtmax = targetlo - weightlo;
	    if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
	    breakflag = 1;
	  }                                      /* wtok = most I can move */
	  for (j = 0, k = 0, wtsum = 0.0; j < numlist && wtsum <= wtok; j++) {
	    i = dotlist[j];
	    if (dotpt[i].x[dim] == med.valuehi) { /* only move if better */
	      if (wtsum + dotpt[i].weight - wtok < wtok - wtsum) {
		dotmark[i] = 0;
		k++;
	      }
	      wtsum += dotpt[i].weight;
	    }
	  }
	  if (breakflag) break;                   /* done if not moved all */
	}

	weightlo += med.wthi;                     /* as close as needed */
	if (fabs(weight - 2.0*weightlo) <= tolerance) break;

	valuemin = med.valuehi;                   /* iterate again */
	markactive = 1;
      }

      else if (weighthi + med.hitotal < targethi) {
                                                  /* upper half TOO SMALL */
	weighthi += med.hitotal;
	valuehalf = med.valuelo;

	if (med.countlo == 1) {                  /* only one dot to move */
	  if (weighthi + med.wtlo < targethi) {  /* move it, keep iterating */
	    if (proc == med.proclo) dotmark[indexlo] = 1;
	  }
	  else {                                 /* only move if beneficial */
	    if (weighthi + med.wtlo - targethi < targethi - weighthi)
	      if (proc == med.proclo) dotmark[indexlo] = 1;
	    break;                               /* all done */
	  }
	}
	else {                                   /* multiple dots to move */
	  breakflag = 0;
	  wtok = 0.0;
	  if (medme.valuelo == med.valuelo) wtok = medme.wtlo;   
	  if (weighthi + med.wtlo >= targethi) { /* can't move all of them */
	    MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,local_comm);
	    wtmax = targethi - weighthi;
	    if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
	    breakflag = 1;
	  }                                      /* wtok = most I can move */
	  for (j = 0, wtsum = 0.0; j < numlist && wtsum <= wtok; j++) {
	    i = dotlist[j];
	    if (dotpt[i].x[dim] == med.valuelo) { /* only move if better */
	      if (wtsum + dotpt[i].weight - wtok < wtok - wtsum) 
		dotmark[i] = 1;
	      wtsum += dotpt[i].weight;
	    }
	  }
	  if (breakflag) break;                   /* done if not moved all */
	}

	weighthi += med.wtlo;                     /* as close as needed */
	if (fabs(weight - 2.0*weighthi) <= tolerance) break;

	valuemax = med.valuelo;                   /* iterate again */
	markactive = 0;
      }

      else                  /* Goldilocks result: both partitions JUST RIGHT */
	break;

      /* shrink the active list */
      
      k = 0;
      for (j = 0; j < numlist; j++) {
	i = dotlist[j];
	if (dotmark[i] == markactive) dotlist[k++] = i;
      }
      numlist = k;

    }

    /* found median */

    if (RCB_STATS) time3 = MPI_Wtime();

    /* store cut info in tree only if I am procmid */

    if (proc == procmid) {
      treept[proc].dim = dim;
      treept[proc].cut = valuehalf;
    }
    
    /* use cut to shrink RCB domain bounding box */

    if (proc < procmid)
      rcbbox->hi[dim] = valuehalf;
    else
      rcbbox->lo[dim] = valuehalf;

    /* outgoing = number of dots to ship to partner */
    /* dottop = number of dots that have never migrated */

    markactive = (proc < procpartner);
    for (i = 0, keep = 0, outgoing = 0; i < dotnum; i++)
      if (dotmark[i] == markactive)
	outgoing++;
      else if (i < dottop)
	keep++;
    dottop = keep;
    
    /* alert partner how many dots I'll send, read how many I'll recv */

    MPI_Send(&outgoing,1,MPI_INT,procpartner,0,MPI_COMM_WORLD);
    incoming = 0;
    if (readnumber) {
      MPI_Recv(&incoming,1,MPI_INT,procpartner,0,MPI_COMM_WORLD,&status);
      if (readnumber == 2) {
	MPI_Recv(&incoming2,1,MPI_INT,procpartner2,0,MPI_COMM_WORLD,&status);
	incoming += incoming2;
      }
    }

    /* check if need to malloc more space */

    dotnew = dotnum - outgoing + incoming;

    if (dotnew > dotmax) {
      allocflag = 1;
      dotmax = overalloc * dotnew;
      if (dotmax < dotnew) dotmax = dotnew;
      dotpt = (struct rcb_dot *) 
	realloc(dotpt,(unsigned) dotmax * sizeof(struct rcb_dot));
      if (dotpt == NULL) rcb_error(dotmax*sizeof(struct rcb_dot));
      *dots = dotpt;
      if (RCB_STATS) counters[6]++;
    }

    if (RCB_STATS) {
      counters[1] += outgoing;
      counters[2] += incoming;
      if (dotnew > counters[3]) counters[3] = dotnew;
      if (dotmax > counters[4]) counters[4] = dotmax;
    }
    
    /* malloc comm send buffer */

    if (outgoing > 0) {
      dotbuf = (struct rcb_dot *)
        malloc((unsigned) outgoing * sizeof(struct rcb_dot));
      if (dotbuf == NULL) rcb_error(outgoing*sizeof(struct rcb_dot));
    }
    else 
      dotbuf = NULL;

    /* fill buffer with dots that are marked for sending */
    /* pack down the unmarked ones */
    
    keep = outgoing = 0;
    for (i = 0; i < dotnum; i++) {
      if (dotmark[i] == markactive)
	memcpy(&dotbuf[outgoing++], &dotpt[i], sizeof(struct rcb_dot));
      else
	memcpy(&dotpt[keep++], &dotpt[i], sizeof(struct rcb_dot));
    }

    /* post receives for dot data */

    if (readnumber > 0) {
      length = incoming * sizeof(struct rcb_dot);
      MPI_Irecv(&dotpt[keep],length,MPI_CHAR,
			   procpartner,1,MPI_COMM_WORLD,&request);
      if (readnumber == 2) {
	keep += incoming - incoming2;
	length = incoming2 * sizeof(struct rcb_dot);
	MPI_Irecv(&dotpt[keep],length,MPI_CHAR,
			      procpartner2,1,MPI_COMM_WORLD,&request2);
      }
    }
    
    /* handshake before sending data to insure recvs have been posted */
    
    if (readnumber > 0) {
      MPI_Send(NULL,0,MPI_INT,procpartner,0,MPI_COMM_WORLD);
      if (readnumber == 2)
	MPI_Send(NULL,0,MPI_INT,procpartner2,0,MPI_COMM_WORLD);
    }
    MPI_Recv(NULL,0,MPI_INT,procpartner,0,MPI_COMM_WORLD,&status);

    /* send dot data to partner */

    length = outgoing * sizeof(struct rcb_dot);
    MPI_Rsend(dotbuf,length,MPI_CHAR,procpartner,1,MPI_COMM_WORLD);
    free(dotbuf);
    
    dotnum = dotnew;

    /* wait until all dots are received */

    if (readnumber > 0) {
      MPI_Wait(&request,&status);
      if (readnumber == 2) MPI_Wait(&request2,&status);
    }

    /* cut partition in half, create new communicators of 1/2 size */

    if (proc < procmid) {
      procupper = procmid - 1;
      i = 0;
    }
    else {
      proclower = procmid;
      i = 1;
    }

    MPI_Comm_split(local_comm,i,proc,&tmp_comm);
    MPI_Comm_free(&local_comm);
    local_comm = tmp_comm;

    if (RCB_STATS) {
      time4 = MPI_Wtime();
      timers[1] += time2 - time1;
      timers[2] += time3 - time2;
      timers[3] += time4 - time3;
    }

  }

  /* have recursed all the way to final single sub-domain */

  /* free all memory used by RCB and MPI */

  MPI_Comm_free(&local_comm);
  MPI_Type_free(&med_type);
  MPI_Type_free(&box_type);
  MPI_Op_free(&box_op);
  MPI_Op_free(&med_op);

  free(dotlist);
  free(dotmark);

  if (RCB_STATS) {
    MPI_Barrier(MPI_COMM_WORLD);
    timestop = time1 = MPI_Wtime();
  }

  /* error checking and statistics */

  if (RCB_CHECK) rcb_check(dotpt,dotnum,*pdotnum,rcbbox);
  if (RCB_STATS) rcb_stats(timestop-timestart,dotpt,dotnum,
			   timers,counters,rcbbox,reuse);

  /* update calling routine parameters */
  
  *pdotnum = dotnum;
  *pdotmax = dotmax;
  *pdottop = dottop;

}



/* ----------------------------------------------------------------------- */

/* error message for malloc/realloc overflow */

void rcb_error(int size)

{
  int proc;

  MPI_Comm_rank(MPI_COMM_WORLD,&proc);
  printf("RCB ERROR: proc = %d could not malloc/realloc %d bytes",proc,size);
  exit(1);
}


/* MPI user-defined reduce operations */

/* min/max merge of each component of a rcb_box */

void rcb_box_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)

{
  int i;
  struct rcb_box *box1,*box2;

  box1 = (struct rcb_box *) in;
  box2 = (struct rcb_box *) inout;

  for (i = 0; i < 3; i++) {
    if (box1->lo[i] < box2->lo[i])
      box2->lo[i] = box1->lo[i];
    if (box1->hi[i] > box2->hi[i])
      box2->hi[i] = box1->hi[i];
  }
}


/* merge median data structure */
/* on input: */
/* in,inout->wtlo, wthi     = # of active dots in both partitions on this proc
             valuelo, valuehi = pos of nearest dot(s) to cut on this proc
             countlo, counthi = # of dot(s) nearest to cut on this proc
             proclo, prochi = not used
   on exit:
   inout->   wtlo, wthi     = total # of active dots in both partitions
             valuelo, valuehi = pos of nearest dot(s) to cut
             countlo, counthi = total # of dot(s) nearest to cut
             proclo, prochi = one unique proc who owns a nearest dot
	                      all procs must get same proclo,prochi
*/

void rcb_median_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)

{
  struct rcb_median *med1,*med2;

  med1 = (struct rcb_median *) in;
  med2 = (struct rcb_median *) inout;
  
  med2->lototal += med1->lototal;
  if (med1->valuelo > med2->valuelo) {
    med2->valuelo = med1->valuelo;
    med2->wtlo = med1->wtlo;
    med2->countlo = med1->countlo;
    med2->proclo = med1->proclo;
  }
  else if (med1->valuelo == med2->valuelo) {
    med2->wtlo += med1->wtlo;
    med2->countlo += med1->countlo;
    if (med1->proclo < med2->proclo)
      med2->proclo = med1->proclo;
  }

  med2->hitotal += med1->hitotal;
  if (med1->valuehi < med2->valuehi) {
    med2->valuehi = med1->valuehi;
    med2->wthi = med1->wthi;
    med2->counthi = med1->counthi;
    med2->prochi = med1->prochi;
  }
  else if (med1->valuehi == med2->valuehi) {
    med2->wthi += med1->wthi;
    med2->counthi += med1->counthi;
    if (med1->prochi < med2->prochi)
      med2->prochi = med1->prochi;
  }
}


/* consistency checks on RCB results */

void rcb_check(struct rcb_dot *dotpt, int dotnum, int dotorig,
	       struct rcb_box *rcbbox)

{
  int i,iflag,proc,nprocs,total1,total2;
  double weight,wtmax,wtmin,wtone,tolerance;

  MPI_Comm_rank(MPI_COMM_WORLD,&proc);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  /* check that total # of dots remained the same */

  MPI_Allreduce(&dotorig,&total1,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&dotnum,&total2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (total1 != total2) {
    if (proc == 0) 
      printf("ERROR: Points before RCB = %d, Points after RCB = %d\n",
	     total1,total2);
  }
  
  /* check that result is load-balanced within log2(P)*max-wt */

  weight = wtone = 0.0;
  for (i = 0; i < dotnum; i++) {
    weight += dotpt[i].weight;
    if (dotpt[i].weight > wtone) wtone = dotpt[i].weight;
  }

  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&wtone,&tolerance,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  /* i = smallest power-of-2 >= nprocs */
  /* tolerance = largest-single-weight*log2(nprocs) */

  for (i = 0; (nprocs >> i) != 0; i++);
  tolerance = tolerance * i * (1.0 + TINY);

  if (wtmax - wtmin > tolerance) {
    if (proc == 0) 
      printf("ERROR: Load-imbalance > tolerance of %g\n",tolerance);
    MPI_Barrier(MPI_COMM_WORLD);
    if (weight == wtmin) printf("  Proc %d has weight = %g\n",proc,weight);
    if (weight == wtmax) printf("  Proc %d has weight = %g\n",proc,weight);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  /* check that final set of points is inside RCB box of each proc */
  
  iflag = 0;
  for (i = 0; i < dotnum; i++) {
    if (dotpt[i].x[0] < rcbbox->lo[0] || dotpt[i].x[0] > rcbbox->hi[0] ||
	dotpt[i].x[1] < rcbbox->lo[1] || dotpt[i].x[1] > rcbbox->hi[1] ||
	dotpt[i].x[2] < rcbbox->lo[2] || dotpt[i].x[2] > rcbbox->hi[2])
      iflag++;
  }
  if (iflag > 0) 
    printf("ERROR: %d points are out-of-box on proc %d\n",iflag,proc);
  
  MPI_Barrier(MPI_COMM_WORLD);
    
}


/* RCB statistics */

void rcb_stats(double timetotal, struct rcb_dot *dotpt,
	       int dotnum, double *timers, int *counters,
	       struct rcb_box *rcbbox, int reuse)

{
  int i,proc,nprocs,sum,min,max;
  double ave,rsum,rmin,rmax;
  double weight,wttot,wtmin,wtmax;

  MPI_Comm_rank(MPI_COMM_WORLD,&proc);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  
  if (proc == 0 &&Iwrite==VERBOSE) printf("RCB total time: %g (secs)\n",timetotal);

  if (proc == 0 &&Iwrite==VERBOSE) printf("RCB Statistics:\n");

  MPI_Barrier(MPI_COMM_WORLD);

  /* distribution info */
  
  for (i = 0, weight = 0.0; i < dotnum; i++) weight += dotpt[i].weight;
  MPI_Allreduce(&weight,&wttot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (proc == 0 &&Iwrite==VERBOSE) {
    printf(" Total weight of dots = %g\n",wttot);
    printf(" Weight on each proc: ave = %g, max = %g, min = %g\n",
	   wttot/nprocs,wtmax,wtmin);
  }

  for (i = 0, weight = 0.0; i < dotnum; i++) 
    if (dotpt[i].weight > weight) weight = dotpt[i].weight;
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  
  if (proc == 0 &&Iwrite==VERBOSE) printf(" Maximum weight of single dot = %g\n",wtmax);

  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2 &&Iwrite==VERBOSE) printf("    Proc %d has weight = %g\n",proc,weight);

  /* counter info */

  MPI_Allreduce(&counters[0],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[0],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[0],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0 &&Iwrite==VERBOSE) 
    printf(" Median iter: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2 &&Iwrite==VERBOSE) 
    printf("    Proc %d median count = %d\n",proc,counters[0]);

  MPI_Allreduce(&counters[1],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[1],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[1],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0 &&Iwrite==VERBOSE) 
    printf(" Send count: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2 &&Iwrite==VERBOSE)
    printf("    Proc %d send count = %d\n",proc,counters[1]);
  
  MPI_Allreduce(&counters[2],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[2],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[2],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0 &&Iwrite==VERBOSE) 
    printf(" Recv count: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2 &&Iwrite==VERBOSE)
    printf("    Proc %d recv count = %d\n",proc,counters[2]);
  
  MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0 &&Iwrite==VERBOSE) 
    printf(" Max dots: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2 &&Iwrite==VERBOSE)
    printf("    Proc %d max dots = %d\n",proc,counters[3]);
  
  MPI_Allreduce(&counters[4],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[4],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[4],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0 &&Iwrite==VERBOSE) 
    printf(" Max memory: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2 &&Iwrite==VERBOSE)
    printf("    Proc %d max memory = %d\n",proc,counters[4]);
  
  if (reuse) {
    MPI_Allreduce(&counters[5],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&counters[5],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&counters[5],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    ave = ((double) sum)/nprocs;
    if (proc == 0 &&Iwrite==VERBOSE) 
      printf(" # of Reuse: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(MPI_COMM_WORLD);
    if (RCB_STATS == 2 &&Iwrite==VERBOSE)
      printf("    Proc %d # of Reuse = %d\n",proc,counters[5]);
  }
  
  MPI_Allreduce(&counters[6],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[6],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[6],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0 &&Iwrite==VERBOSE) 
    printf(" # of OverAlloc: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2 &&Iwrite==VERBOSE)
    printf("    Proc %d # of OverAlloc = %d\n",proc,counters[6]);

  /* timer info */
  
  MPI_Allreduce(&timers[0],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[0],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[0],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0 &&Iwrite==VERBOSE) 
    printf(" Start-up time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2&&Iwrite==VERBOSE)
    printf("    Proc %d start-up time = %g\n",proc,timers[0]);
  
  MPI_Allreduce(&timers[1],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[1],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[1],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0&&Iwrite==VERBOSE) 
    printf(" Pre-median time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2&&Iwrite==VERBOSE)
    printf("    Proc %d pre-median time = %g\n",proc,timers[1]);
  
  MPI_Allreduce(&timers[2],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[2],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[2],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0&&Iwrite==VERBOSE) 
    printf(" Median time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2&&Iwrite==VERBOSE)
    printf("    Proc %d median time = %g\n",proc,timers[2]);
  
  MPI_Allreduce(&timers[3],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[3],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[3],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0&&Iwrite==VERBOSE) 
    printf(" Comm time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2&&Iwrite==VERBOSE)
    printf("    Proc %d comm time = %g\n",proc,timers[3]);
  
  /* RCB boxes for each proc */
  
  if (RCB_STATS == 2) {
    if (proc == 0&&Iwrite==VERBOSE) printf(" RCB sub-domain boxes:\n");
    for (i = 0; i < 3; i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (proc == 0&&Iwrite==VERBOSE) printf("    Dimension %d\n",i+1);
      MPI_Barrier(MPI_COMM_WORLD);
      if (Iwrite==VERBOSE)printf("      Proc = %d: Box = %g %g\n",
	     proc,rcbbox->lo[i],rcbbox->hi[i]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sort_int_array(int n, int ra[])
/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*       modified to include int array and to number from 0
*
*       Sorts the array ra[1,..,n] in ascending numerical order using heapsort
*       algorithym.
*
*/

{
  int l,j,ir,i;
  int rra;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
        ra[1]=rra;
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir)     {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        j += (i=j);
      }
      else j=ir+1;
    }
    ra[i]=rra;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double set_weight_for_node(int inode)
/* This routine returns a weight for a given node, which should be           */
/* proportional to the amount of work that needs to be done in the fill and  */
/* solve for this node. This is currently a function of Zero_density_TF      */
/* Weights must be greater than zero! (weight > 0)                           */

{
  int icomp,idim,ijk[3],inode_box;
  double weight = 1.0, max, fac;
  
  /* Add extra weight near periodic/Reflective boundaries */

  node_to_ijk(inode,ijk);
  inode_box=node_to_node_box(inode);

  for (icomp=0; icomp<Ncomp; icomp++){

     if (Type_attr != NONE) max = Cut_ff[icomp][icomp];
     else                  max = 0.5*Sigma_ff[icomp][icomp];

     for (idim=0; idim<Ndim; idim++){
       if (  ( Type_bc[idim][0] == REFLECT && ijk[idim]*Esize_x[idim] <= max )   ||   
             ( Type_bc[idim][1] == REFLECT && 
               (Nodes_x[idim] - ijk[idim])*Esize_x[idim] <= max ) ) weight += 0.2;

       if (  ( Type_bc[idim][0] == PERIODIC && ijk[idim]*Esize_x[idim] <= max )   ||   
             ( Type_bc[idim][1] == PERIODIC && 
               (Nodes_x[idim] - ijk[idim])*Esize_x[idim] <= max ) ) weight += 0.0;

        if (  ( Type_bc[idim][0] == IN_WALL && ijk[idim]*Esize_x[idim] <= 2.0*max ) ||
              ( Type_bc[idim][1] == IN_WALL && 
               (Nodes_x[idim] - ijk[idim])*Esize_x[idim] <= 2.0*max ) ) weight += 0.5;
    }
  }
      
  /* decrement weight by 0.95/Ncomp for each component which has Zero_density*/

  fac = 1.0;
  for (icomp = 0; icomp < Ncomp; icomp++) 
    if (Zero_density_TF[inode_box][icomp])  fac -= 0.95/((float) Ncomp);   
  weight *= fac;

  /* Adjust weights for mesh coarsening -- trivial and reduced nsten nodes */

  if (Mesh_coarsening || L1D_bc) {
    if (Mesh_coarsen_flag[inode_box] < 0 && Mesh_coarsen_flag[inode_box]>= -3) {
        weight = 0.02;  /* trivial residual */
    }
    else if (Mesh_coarsen_flag[inode_box] == FLAG_1DBC || Mesh_coarsen_flag[inode_box] == FLAG_BULK)
        weight = 0.001;
    else if (Mesh_coarsen_flag[inode_box] > 0) {
 
       /* divide by coarsening of residual stencil */ 
       weight /= 0.8 * POW_INT(2,Ndim*Mesh_coarsen_flag[inode_box]);

       /* conservatively assume that Jacobian stencils are 1 finer */ 
       weight /= POW_INT(2,(int)(Ndim*(Mesh_coarsen_flag[inode_box]-1.0)));
    }
  }

  return weight;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


void load_balance(int flag, double *fill_time, int *N_update, int **update)
/* This routine re-loadbalances the nodes, as stored by unknown numbers in   */
/* *update, by calling the recursive coordinate bisection (rcb) routine */
/* supplied by Steve Plimpton. The rest of this routine is just preparing    */
/* the data to send to the rcb routine, and stuffing the results back into   */
/* our format (*N_update, *update, and Nnodes_per_proc).           */
/* If fill_time != NULL, then this array of timings is used as nodal weights;*/
/* otherwise, the routine set_weight_for_node is used.                       */ 

{
  int loc_inode, inode, iunk, ijk[3];
  struct rcb_dot *dots; /* Array of the structure that holds x[3], the       */ 
                        /* coordinates of nodes owned by this proc, and      */
                        /* weights, the weight of that node to balance       */
  int              pdotmax;    /* max # of dots array can hold, may extend */
  int              pdottop;    /* dots >= this index are new */
  struct rcb_box   rcbbox;     /* bounding box of final RCB sub-domain */
  struct rcb_tree  *treept;    /* tree of RCB cuts - only single cut on exit */
  double navg;
  int    nmin, nmax;
  int idim;

  nmax = gmax_int(*N_update);
  nmin = gmin_int(*N_update);
  navg = gsum_double((double)*N_update) / (double) Num_Proc;
  if (Proc == 0 &&Iwrite==VERBOSE) {
    printf("\n+++++++++++++++++++++++++++++++++++++");
    printf("+++++++++++++++++++++++++++++++++++++++\n");
    printf("Calling RCB routine for load balancing\n\n");
    printf("\tBefore load_balance: Max=%d  Min=%d  Avg=%g (Unknowns per Proc)\n",
             nmax, nmin, navg);
  }

  /* Start MPI if this is the first time in this routine */

  /* Allocate and set up the array of struct of coordinates and weights */

  dots = (struct rcb_dot*) array_alloc
                                (1, Nnodes_per_proc, sizeof(struct rcb_dot));

  for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++) {

    if (MATRIX_FILL_NODAL) inode = (*update)[loc_inode * Nunk_per_node] / Nunk_per_node;
    else                   inode = (*update)[loc_inode];

    node_to_position(inode, dots[loc_inode].x);

    if (fill_time == NULL)
      if (flag == 0) dots[loc_inode].weight = 1.;
      else  dots[loc_inode].weight = set_weight_for_node(inode);
    else
      dots[loc_inode].weight = fill_time[loc_inode];
  }

  /* rcb routine assumes 3D, must pad the unused coords with zeroes */

  if (Ndim == 1) {
    for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++) {
      dots[loc_inode].x[1] = 0.0;
      dots[loc_inode].x[2] = 0.0;
    }
  }
  else if (Ndim == 2) {
    for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++) {
      dots[loc_inode].x[2] = 0.0;
    }
  }

  /* set other parameters */

  pdotmax = Nnodes_per_proc;
  treept = (struct rcb_tree *) array_alloc(1, Num_Proc, sizeof(struct rcb_tree));

  /* call rcb routine, which returns new values of dots, Nnodes_per_proc, etc */

 if (Load_Bal_Flag != LB_LINEAR) rcb(&dots, &(Nnodes_per_proc), &pdotmax, 1, 1.5, 0, &pdottop, &rcbbox, treept);
  /* Print out statistics before load balancing */

  nmax = gmax_int(*N_update);
  nmin = gmin_int(*N_update);
  navg = gsum_double((double)*N_update) / (double)Num_Proc;
  if (Proc == 0&&Iwrite==VERBOSE)
    printf("\tBefore load_balance: Max %d  Min %d  Ave %g Unknowns per Proc\n",
             nmax, nmin, navg);

  /* Reset  arrays for unknowns owned by this proc and print new stats*/

  safe_free((void *) update);
  *N_update = Nunk_per_node * Nnodes_per_proc;
  *update = (int *) array_alloc(1, *N_update, sizeof(int));
  /*printf("in dft_ldbal: *update allocated to %d length\n",*N_update);*/

  nmax = gmax_int(*N_update);
  nmin = gmin_int(*N_update);
  navg = gsum_double((double)*N_update) / (double)Num_Proc;
  if (Proc == 0 && Iwrite == VERBOSE)
    printf("\tAfter load_balance: Max=%d  Min=%d  Avg=%g (Unknowns per Proc)\n",
             nmax, nmin, navg);


  /* Initialize max and min dimensions for each dimension */
 
  for (idim=0; idim<Ndim; idim++){
    Min_IJK[idim] = Nnodes;
    Max_IJK[idim] = 0;
  }
  for (idim=Ndim; idim<3; idim++) {
    Max_IJK[idim] = Min_IJK[2] = 0;
  }

  /* Translate new dots array into node and unknown numbers, put in *update */
  for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++) {

     inode = position_to_node(dots[loc_inode].x);

     for (iunk = 0; iunk < Nunk_per_node; iunk++){
        if (MATRIX_FILL_NODAL) (*update)[loc_find(iunk,loc_inode,LOCAL_N)] = inode * Nunk_per_node + iunk;
        else                   (*update)[loc_find(iunk,loc_inode,LOCAL_N)] = inode + Nnodes_per_proc*iunk;
     }

     /* Calculate max 'n min range of nodes for each dimension */

     node_to_ijk(inode, ijk);
     for (idim=0; idim<Ndim; idim++){
       if (ijk[idim] < Min_IJK[idim]) Min_IJK[idim] = ijk[idim];
       if (ijk[idim] > Max_IJK[idim]) Max_IJK[idim] = ijk[idim];
     }
  }

  /* print out bounding boxes */
  if (Iwrite == VERBOSE && Num_Proc>1) {
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Proc %2d: Position x: %g = %g  y: %g = %g\n",Proc,
            rcbbox.lo[0], rcbbox.hi[0], rcbbox.lo[1], rcbbox.hi[1]);

  /* get the Minimum and maximum ijk corresponding to the edges of the
     domain on this processor */

    printf("Proc %2d: Nnodes_per_proc = %d\n",Proc,Nnodes_per_proc);
    printf("Proc %2d: Position x: %g - %g  y: %g - %g\n",Proc,
            rcbbox.lo[0], rcbbox.hi[0], rcbbox.lo[1], rcbbox.hi[1]);
  }
/*********
     for (idim=0; idim<Ndim; idim++){
        Min_IJK[idim] =    (int)     ((rcbbox.lo[idim]+0.5*Size_x[idim])/Esize_x[idim]-0.25);
        Max_IJK[idim] = round_to_int ((rcbbox.hi[idim]+0.5*Size_x[idim])/Esize_x[idim]-0.25);
     }
     if (Ndim < 3){
         Max_IJK[2] = Min_IJK[2] = 0;
       if (Ndim < 2)
         Max_IJK[1] = Min_IJK[1] = 0;
     }
     printf("the bounds on the boxes for load balancing are:\n");
     for (idim=0; idim<Ndim; idim++){
        printf("Proc: %d  idim: %d :   Min_IJK: %d  Max_IJK: %d\n",
                Proc,idim,Min_IJK[idim],Max_IJK[idim]);
     }
**********/

  /* Finally, sort the new *update array in ascending order */

  sort_int_array(*N_update, *update - 1);

  /* Free the dots and treept arrays */

  safe_free((void *) &dots);
  safe_free((void *) &treept);

  MPI_Barrier(MPI_COMM_WORLD);
  if (Proc == 0 && Iwrite==VERBOSE) {
    printf("+++++++++++++++++++++++++++++++++++++");
    printf("+++++++++++++++++++++++++++++++++++++++\n");
  }
}

/******************************************************************************/
