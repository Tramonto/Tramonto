#define INT_GLOBAL_VEC    1000    /* see rf_comm.c            */
#define DOUBLE_GLOBAL_VEC 5000
#define INT_VEC_LENGTH    9000
#define FLAG_VEC_LENGTH    999
#define VACC_LOADS        1200
#define FLAG_DONE_READ    3500

/* System Include files */
#include <assert.h>
#include "dft_globals_const.h"
#include "rf_allo.h"

/* Here's a place to store communications routines for parallel processing */

/***************************************************************************/
void gather_global_vec(double *loc_vec, int *loc_index,int N_loc,
                                              double *global_vec)

/*
 * This routine was created by modifying an MPSalsa routine by the same name.
 *
 * NOTE: This routine assumes that there is enough memory on proc zero to
 * malloc a vector for an int and a double of length Num_Nodes.
 *
 *  Parameter List:
 *
 *       loc_vec[]    ==    vector containing the solution vector
 *                          index by the local node number on the
 *                          current processor.
 *
 *       local_index[]    ==    vector containing global node numbers
 *                              of each of the local processor nodes.
 * 
 *       N_loc        == Length of loc_vec and loc_index arrays.
 *
 *       global_vec[] ==    Output: global vector which has been 
 *                          gathered from all procs.
 *
 */

{
  int    max_proc_num_nodes;
  int    itype, rtype, ltype, flagtype;
  int    zero, flag=0;
  int    N_jproc,proc0;
  int type2,flag_done_read;

  int *loc_index_jproc;
  double *vec_jproc;
  int    from_p, i, cflag, read_flag;
  MPI_AZRequest st_wait;
  int st_write;

  int result;
  
  int max_length,msg_length,istart,n_msgs,imsg;

  char   *yo_error = "print_global_vec: ERROR on node ";

  /*
   * Find the maximum value of N_loc so that Node 0 knows how long the 
   * loc_index_jproc and vec_jproc arrays should be.
   */

  printf("calling AZ_gmax_int\n");
  max_proc_num_nodes = AZ_gmax_int(N_loc, Aztec.proc_config);
  

  if (Proc==0) {

    loc_index_jproc = (int *) array_alloc(1, max_proc_num_nodes, sizeof(int));
    vec_jproc = (double *) array_alloc(1, max_proc_num_nodes, sizeof(double));

    for (i=0; i<N_loc; i++) global_vec[loc_index[i]] = loc_vec[i];
  }

  /*
   * Loop through all procs sending ....
   */


  max_length = 100000*sizeof(double);

  flagtype = FLAG_VEC_LENGTH;
  flag_done_read = 0;
  type2 = FLAG_DONE_READ;
  if (Proc == 0) {

    for (from_p = Num_Proc-1; from_p > 0; from_p--) {

  printf("calling md_wrap_write: proc: %d from proc: %d\n",Proc,from_p);
      /* send a flag to from_p so that this processor knows to send its message */
      if (md_wrap_write((void *) &flag, sizeof(int), from_p, flagtype, &st_write) != 0) {
        (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
        (void) fprintf(stderr, "write failed: from_p\n");
      }
  printf("after md_wrap_write: proc: %d from proc: %d\n",Proc,from_p);

  printf("calling md_wrap_wait: proc: %d from proc: %d\n",Proc,from_p);
      itype = INT_VEC_LENGTH + from_p;
      read_flag = md_wrap_wait((void *) &N_jproc, sizeof(int), &from_p, &itype, 
                                        &cflag, &st_wait);
  printf("after md_wrap_wait: proc: %d from proc: %d\n",Proc,from_p);
 
      if (N_jproc*sizeof(double) > max_length){
           n_msgs = N_jproc*sizeof(double)/max_length;
           if (n_msgs*max_length/sizeof(double) != N_jproc) n_msgs++;
      }
      else n_msgs = 1;

      for (imsg = 0; imsg<n_msgs; imsg++){

         itype = INT_VEC_LENGTH + from_p;
         read_flag = md_wrap_wait((void *) &msg_length, sizeof(int), &from_p, &itype, 
                                           &cflag, &st_wait);
         istart = imsg*max_length/sizeof(double);

         itype = INT_GLOBAL_VEC + from_p;
         read_flag = md_wrap_wait((void *) &(loc_index_jproc[istart]), 
                                           msg_length*sizeof(int), 
                                      &from_p, &itype, &cflag, &st_wait);
      
         itype = DOUBLE_GLOBAL_VEC + from_p;
         read_flag = md_wrap_wait((void *) &(vec_jproc[istart]), 
                                  msg_length*sizeof(double), 
                                 &from_p, &itype, &cflag, &st_wait);

         result = md_write((char *) &flag_done_read, 
                  sizeof(int), from_p, type2, &cflag); 
	 assert( result == 0);
      }
     
      for (i=0; i<N_jproc; i++) 
         if (fabs(vec_jproc[i]) > 1.0e-10) global_vec[loc_index_jproc[i]] = vec_jproc[i];

    } /* processor loop */
  } /* proc == 0 */

  else {
    rtype = DOUBLE_GLOBAL_VEC;
    itype = INT_GLOBAL_VEC;
    ltype = INT_VEC_LENGTH;

    /* read a flag from proc zero */
    zero = 0;
    proc0 = 0;
  printf("before md_wrap_wait: proc: %d \n",Proc);
    read_flag = md_wrap_wait((void *) &flag, sizeof(int), &zero, &flagtype, 
                                                              &cflag, &st_wait);
  printf("after md_wrap_wait: proc: %d \n",Proc);

    if (md_wrap_write((void *) &N_loc, sizeof(int), 0, ltype+Proc, &st_write) != 0) {
      (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
      (void) fprintf(stderr, "write failed: N_loc\n");
    }

    if (N_loc*sizeof(double) > max_length){
         n_msgs = N_loc*sizeof(double)/max_length;
         if (n_msgs*max_length/sizeof(double) != N_loc) n_msgs++;
    }
    else n_msgs = 1;

    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(double);

       if(n_msgs == 1)        msg_length = N_loc;
       else if(imsg<n_msgs-1)msg_length = max_length/sizeof(double);
       else                   msg_length = N_loc - (n_msgs-1)*max_length/sizeof(double);

       if (md_wrap_write((void *) &msg_length, sizeof(int), 0, ltype+Proc, &st_write) != 0) {
         (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
         (void) fprintf(stderr, "write failed: N_loc\n");
       }

       if (md_wrap_write((void *) &(loc_index[istart]), 
                              msg_length*sizeof(int), 0, 
                              itype+Proc, &st_write) != 0) {
         (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
         (void) fprintf(stderr, "write failed: loc_index\n");
       }

       if (md_wrap_write((void *) &(loc_vec[istart]), 
                          msg_length*sizeof(double), 0, 
                            rtype+Proc, &st_write) != 0) {
         (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
         (void) fprintf(stderr, "write failed: loc_vec\n");
       }
       result = md_read((char *) &flag_done_read, 
                sizeof(int), &proc0, &type2, &cflag);
       assert( (result-sizeof(int)) == 0 );
    }

  } /* all procs != 0 */

  if (Proc==0){
    safe_free((void *) &loc_index_jproc);
    safe_free((void *) &vec_jproc);
  }
} 
/***************************************************************************/
void gather_global_vec_int(int *loc_vec, int *loc_index,int N_loc,
                                              int *global_vec)

/*
 * This routine was created by modifying an MPSalsa routine by the same name.
 *
 * NOTE: This routine assumes that there is enough memory on proc zero to
 * malloc a vector for an int and a int of length Num_Nodes.
 *
 *  Parameter List:
 *
 *       loc_vec[]    ==    vector containing the solution vector
 *                          index by the local node number on the
 *                          current processor.
 *
 *       local_index[]    ==    vector containing global node numbers
 *                              of each of the local processor nodes.
 * 
 *       N_loc        == Length of loc_vec and loc_index arrays.
 *
 *       global_vec[] ==    Output: global vector which has been 
 *                          gathered from all procs.
 *
 */

{
  int    max_proc_num_nodes;
  int    itype, rtype, ltype, flagtype;
  int    zero, flag=0;
  int    N_jproc,proc0;
  int type2,flag_done_read;

  int *loc_index_jproc;
  int *vec_jproc;
  int    from_p, i, cflag, read_flag;
  MPI_AZRequest st_wait;
  int st_write;

  int result;

  int max_length,msg_length,istart,n_msgs,imsg;

  char   *yo_error = "print_global_vec: ERROR on node ";

  /*
   * Find the maximum value of N_loc so that Node 0 knows how long the 
   * loc_index_jproc and vec_jproc arrays should be.
   */

  max_proc_num_nodes = AZ_gmax_int(N_loc, Aztec.proc_config);
  

  if (Proc==0) {

    loc_index_jproc = (int *) array_alloc(1, max_proc_num_nodes, sizeof(int));
    vec_jproc = (int *) array_alloc(1, max_proc_num_nodes, sizeof(int));

    for (i=0; i<N_loc; i++) global_vec[loc_index[i]] = loc_vec[i];
  }

  /*
   * Loop through all procs sending ....
   */


  max_length = 100000*sizeof(int);

  flagtype = FLAG_VEC_LENGTH;
  flag_done_read = 0;
  type2 = FLAG_DONE_READ;
  if (Proc == 0) {

    for (from_p = Num_Proc-1; from_p > 0; from_p--) {

      /* send a flag to from_p so that this processor knows to send its message */
      if (md_wrap_write((void *) &flag, sizeof(int), from_p, flagtype, &st_write) != 0) {
        (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
        (void) fprintf(stderr, "write failed: from_p\n");
      }

      itype = INT_VEC_LENGTH + from_p;
      read_flag = md_wrap_wait((void *) &N_jproc, sizeof(int), &from_p, &itype,
			       &cflag, &st_wait);
 
      if (N_jproc*sizeof(int) > max_length){
           n_msgs = N_jproc*sizeof(int)/max_length;
           if (n_msgs*max_length/sizeof(int) != N_jproc) n_msgs++;
      }
      else n_msgs = 1;

      for (imsg = 0; imsg<n_msgs; imsg++){

         itype = INT_VEC_LENGTH + from_p;
         read_flag = md_wrap_wait((void *) &msg_length, sizeof(int), &from_p, &itype, 
                                           &cflag, &st_wait);
         istart = imsg*max_length/sizeof(int);

         itype = INT_GLOBAL_VEC + from_p;
         read_flag = md_wrap_wait((void *) &(loc_index_jproc[istart]), 
                                           msg_length*sizeof(int), 
                                      &from_p, &itype, &cflag, &st_wait);
      
         itype = DOUBLE_GLOBAL_VEC + from_p;
         read_flag = md_wrap_wait((void *) &(vec_jproc[istart]), 
                                  msg_length*sizeof(int), 
                                 &from_p, &itype, &cflag, &st_wait);

         result = md_write((char *) &flag_done_read, 
                  sizeof(int), from_p, type2, &cflag); 
	 assert( result == 0);
      }
     
      for (i=0; i<N_jproc; i++) 
         if (abs(vec_jproc[i]) > 1.0e-10) global_vec[loc_index_jproc[i]] = vec_jproc[i];

    } /* processor loop */
  } /* proc == 0 */

  else {
    rtype = DOUBLE_GLOBAL_VEC;
    itype = INT_GLOBAL_VEC;
    ltype = INT_VEC_LENGTH;

    /* read a flag from proc zero */
    zero = 0;
    proc0 = 0;
    read_flag = md_wrap_wait((void *) &flag, sizeof(int), &zero, &flagtype, 
                                                              &cflag, &st_wait);

    if (md_wrap_write((void *) &N_loc, sizeof(int), 0, ltype+Proc, &st_write) != 0) {
      (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
      (void) fprintf(stderr, "write failed: N_loc\n");
    }

    if (N_loc*sizeof(int) > max_length){
         n_msgs = N_loc*sizeof(int)/max_length;
         if (n_msgs*max_length/sizeof(int) != N_loc) n_msgs++;
    }
    else n_msgs = 1;

    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(int);

       if(n_msgs == 1)        msg_length = N_loc;
       else if(imsg<n_msgs-1)msg_length = max_length/sizeof(int);
       else                   msg_length = N_loc - (n_msgs-1)*max_length/sizeof(int);

       if (md_wrap_write((void *) &msg_length, sizeof(int), 0, ltype+Proc, &st_write) != 0) {
         (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
         (void) fprintf(stderr, "write failed: N_loc\n");
       }

       if (md_wrap_write((void *) &(loc_index[istart]), 
                              msg_length*sizeof(int), 0, 
                              itype+Proc, &st_write) != 0) {
         (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
         (void) fprintf(stderr, "write failed: loc_index\n");
       }

       if (md_wrap_write((void *) &(loc_vec[istart]), 
                          msg_length*sizeof(int), 0, 
                            rtype+Proc, &st_write) != 0) {
         (void) fprintf(stderr, "%s%d\n", yo_error, Proc);
         (void) fprintf(stderr, "write failed: loc_vec\n");
       }
       result = md_read((char *) &flag_done_read, 
                sizeof(int), &proc0, &type2, &cflag);
       assert( (result-sizeof(int)) == 0 );
    }

  } /* all procs != 0 */

  if (Proc==0){
    safe_free((void *) &loc_index_jproc);
    safe_free((void *) &vec_jproc);
  }
} 
/****************************************************************************/
/* share_global_int_vec:  This routine was provided by Karen Devine.  
                    Given an integer vector that has been partially 
                    assembled on each processor and needs to be global, 
                    this routine assembles the global array in an 
                    efficient manner. A flag of -1 indicates the vector has
                    not yet been filled. */

void share_global_int_vec(int vec_length,int *int_vec)
{
int   j;
int   nprocs;       /* total number of processors */
int   size;
int   type;         /* type of next message */
int   type2;        /* type of next message */
int   partner;      /* processor I exchange with */
int   mask;         /* bit pattern identifying partner */
int   hbit;         /* largest nonzero bit in nprocs */
int   nprocs_small; /* largest power of 2 <= nprocs */
int   cflag;        /* dummy argument for compatability */
int   *tmp;         /* a temporary array for new info */

int max_length,msg_length,istart,n_msgs,imsg,flag_done_read;

int result;
 
  /*********************** first executable statment *****************/
 
  nprocs = Aztec.proc_config[AZ_N_procs];
  tmp = (int *) array_alloc (1, vec_length, sizeof(int));

  type = VACC_LOADS;
  type2 = FLAG_DONE_READ;
  size = vec_length * sizeof(int);
  flag_done_read = 0;
  type2 = VACC_LOADS;

  /*  Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  max_length = 50000*sizeof(int);
  if (size > max_length){
       n_msgs = size/max_length;
       if (n_msgs*max_length != size) n_msgs++;
  }
  else n_msgs = 1;

  partner = Proc ^ nprocs_small;
  if (Proc & nprocs_small) {

    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(int);

       if(n_msgs == 1)        msg_length = size;
       else if(imsg<n_msgs-1) msg_length = max_length;
       else                   msg_length = size - (n_msgs-1)*max_length;

       if (   md_write((char *) &(int_vec[istart]), 
                       msg_length, partner, type, &cflag)   ) {
         printf("1 WRITE ERROR ... see share_global_int_vec.\n");
         exit(-1);
       }
       result = md_read((char *) &flag_done_read, 
              sizeof(int), &partner, &type2, &cflag);
       assert( (result - sizeof(int)) == 0 );
       
    }
  }  
  else if ((Proc+nprocs_small) < nprocs) {
    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(int);

       if(n_msgs == 1)        msg_length = size;
       else if(imsg<n_msgs-1) msg_length = max_length;
       else                   msg_length = size - (n_msgs-1)*max_length;

       if (   md_read((char *) &(tmp[istart]), 
              msg_length, &partner, &type, &cflag) != msg_length) {
         printf("1 READ ERROR ... see share_global_int_vec.\n");
         exit(-1);
       }
       result = md_write((char *) &flag_done_read, 
                sizeof(int), partner, type2, &cflag); 
       assert( result == 0);

    }
    for (j = 0; j < vec_length; j++) {
      if (int_vec[j] == -1 && tmp[j] != -1) int_vec[j] = tmp[j];
    }
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(Proc & nprocs_small)) {
    for (mask=nprocs_small>>1; mask; mask>>=1) {
      type++;
      partner = Proc ^ mask;

      if (Proc < partner) {
      for (imsg=0; imsg<n_msgs; imsg++) {

         istart = imsg*max_length/sizeof(int);

         if(n_msgs == 1)        msg_length = size;
         else if(imsg<n_msgs-1) msg_length = max_length;
         else                   msg_length = size - (n_msgs-1)*max_length;

        if (md_write((char *) &(int_vec[istart]), 
             msg_length, partner, type, &cflag)) {
           printf("2 WRITE ERROR ... see share_global_int_vec.\n");
           exit(-1);
        }  
        result = md_read((char *) &flag_done_read, 
                sizeof(int), &partner, &type2, &cflag);
	 assert( (result - sizeof(int)) == 0 );
      }

      for (imsg=0; imsg<n_msgs; imsg++) {

         istart = imsg*max_length/sizeof(int);

         if(n_msgs == 1)        msg_length = size;
         else if(imsg<n_msgs-1) msg_length = max_length;
         else                   msg_length = size - (n_msgs-1)*max_length;

         if (md_read((char *) &(tmp[istart]), 
             msg_length, &partner, &type, &cflag) != msg_length) {
            printf("2 READ ERROR ... see share_global_int_vec.\n");
            exit(-1);
         }  
         result = md_write((char *) &flag_done_read, 
                sizeof(int), partner, type2, &cflag); 
	 assert( result == 0);
      }

      for (j = 0; j < vec_length; j++) {
        if (int_vec[j] == -1 && tmp[j] != -1) int_vec[j] = tmp[j];
      }
      } /* end of if Proc<partner */
      else {
      for (imsg=0; imsg<n_msgs; imsg++) {

         istart = imsg*max_length/sizeof(int);

         if(n_msgs == 1)        msg_length = size;
         else if(imsg<n_msgs-1) msg_length = max_length;
         else                   msg_length = size - (n_msgs-1)*max_length;

         if (md_read((char *) &(tmp[istart]), 
             msg_length, &partner, &type, &cflag) != msg_length) {
            printf("2 READ ERROR ... see share_global_int_vec.\n");
            exit(-1);
         }  
         result = md_write((char *) &flag_done_read, 
                sizeof(int), partner, type2, &cflag); 
	 assert( result == 0);
      }

      for (imsg=0; imsg<n_msgs; imsg++) {

         istart = imsg*max_length/sizeof(int);

         if(n_msgs == 1)        msg_length = size;
         else if(imsg<n_msgs-1) msg_length = max_length;
         else                   msg_length = size - (n_msgs-1)*max_length;

        if (md_write((char *) &(int_vec[istart]), 
             msg_length, partner, type, &cflag)) {
           printf("2 WRITE ERROR ... see share_global_int_vec.\n");
           exit(-1);
        }  
        result = md_read((char *) &flag_done_read, 
                sizeof(int), &partner, &type2, &cflag);
	 assert( (result - sizeof(int)) == 0);
      }

      for (j = 0; j < vec_length; j++) {
        if (int_vec[j] == -1 && tmp[j] != -1) int_vec[j] = tmp[j];
      }
      } /* end of if Proc<partner */
    }
  }  
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = Proc ^ nprocs_small;
  if (Proc & nprocs_small) {
    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(int);

       if(n_msgs == 1)        msg_length = size;
       else if(imsg<n_msgs-1) msg_length = max_length;
       else                   msg_length = size - (n_msgs-1)*max_length;

       if (md_read((char *) &(int_vec[istart]), msg_length, &partner, &type, &cflag) != msg_length) {
         printf("3 READ ERROR ... see share_global_int_vec.\n");
         exit(-1);
       }
       result = md_write((char *) &flag_done_read, 
                sizeof(int), partner, type2, &cflag); 
       assert( result == 0);
    }
  }  
  else if ((Proc+nprocs_small) < nprocs) {
    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(int);

       if(n_msgs == 1)        msg_length = size;
       else if(imsg<n_msgs-1) msg_length = max_length;
       else                   msg_length = size - (n_msgs-1)*max_length;

       if (md_write((char *) &(int_vec[istart]), msg_length, partner, type, &cflag)) {
         printf("3 WRITE ERROR ... see share_global_int_vec.\n");
         exit(-1);
       }
       result = md_read((char *) &flag_done_read, 
                sizeof(int), &partner, &type2, &cflag);
       assert( (result - sizeof(int)) == 0);
    }
  }  

  safe_free((void *) &tmp);
}
/****************************************************************************/
/* share_global_double_vec:  This routine was provided by Karen Devine.  
                    Given an integer vector that has been partially 
                    assembled on each processor and needs to be global, 
                    this routine assembles the global array in an 
                    efficient manner. A flag of -1 indicates the vector has
                    not yet been filled. */

void share_global_double_vec(int vec_length,double *double_vec)
{
int   j;
int   nprocs;       /* total number of processors */
int   size;
int   type;         /* type of next message */
int   type2;         /* type of next message */
int   flag_done_read;         /* type of next message */
int   partner;      /* processor I exchange with */
int   mask;         /* bit pattern identifying partner */
int   hbit;         /* largest nonzero bit in nprocs */
int   nprocs_small; /* largest power of 2 <= nprocs */
int   cflag;        /* dummy argument for compatability */
double   *tmp;         /* a temporary array for new info */
int max_length,msg_length,istart,n_msgs,imsg;

int result;
 
  /*********************** first executable statment *****************/
 
  nprocs = Aztec.proc_config[AZ_N_procs];
  tmp = (double *) array_alloc (1, vec_length, sizeof(double));

  type = VACC_LOADS;
  type2 = FLAG_DONE_READ;
  size = vec_length * sizeof(double);

  /*  Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  max_length = 100000*sizeof(double);
  if (size > max_length){
       n_msgs = size/max_length;
       if (n_msgs*max_length != size) n_msgs++;
  }

  else n_msgs = 1;
  partner = Proc ^ nprocs_small;
  if (Proc & nprocs_small) {
    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(double);

       if(n_msgs == 1)        msg_length = size;
       else if(imsg<n_msgs-1) msg_length = max_length;
       else                   msg_length = size - (n_msgs-1)*max_length;

       if (md_write((char *) &(double_vec[istart]), 
            msg_length, partner, type, &cflag)) {
         printf("1 WRITE ERROR ... see share_global_double_vec.\n");
         exit(-1);
       }
       result = md_read((char *) &flag_done_read, 
              sizeof(int), &partner, &type2, &cflag);
       assert( (result - sizeof(int)) == 0);
    }
  }  
  else if ((Proc+nprocs_small) < nprocs) {
    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(double);

       if(n_msgs == 1)        msg_length = size;
       else if(imsg<n_msgs-1) msg_length = max_length;
       else                   msg_length = size - (n_msgs-1)*max_length;

       if (md_read((char *) &(tmp[istart]), 
           msg_length, &partner, &type, &cflag) != msg_length) {
         printf("1 READ ERROR ... see share_global_double_vec.\n");
         exit(-1);
       }
       result = md_write((char *) &flag_done_read, 
                sizeof(int), partner, type2, &cflag); 
       assert( result == 0);
    }
    for (j = 0; j < vec_length; j++) {
      if (double_vec[j] == UNINIT_VEC && tmp[j] != UNINIT_VEC) double_vec[j] = tmp[j];
    }
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(Proc & nprocs_small)) {
    for (mask=nprocs_small>>1; mask; mask>>=1) {
      type++;
      partner = Proc ^ mask;

      if (Proc < partner) {
      for (imsg=0; imsg<n_msgs; imsg++) {

         istart = imsg*max_length/sizeof(double);

         if(n_msgs == 1)        msg_length = size;
         else if(imsg<n_msgs-1) msg_length = max_length;
         else                   msg_length = size - (n_msgs-1)*max_length;

        if (md_write((char *) &(double_vec[istart]), 
            msg_length, partner, type, &cflag)) {
           printf("2 WRITE ERROR ... see share_global_double_vec.\n");
           exit(-1);
        }  
        result = md_read((char *) &flag_done_read, 
              sizeof(int), &partner, &type2, &cflag);
	 assert( (result - sizeof(int)) == 0);
      }

      for (imsg=0; imsg<n_msgs; imsg++) {

         istart = imsg*max_length/sizeof(double);

         if(n_msgs == 1)        msg_length = size;
         else if(imsg<n_msgs-1) msg_length = max_length;
         else                   msg_length = size - (n_msgs-1)*max_length;

         if (md_read((char *) &(tmp[istart]), 
             msg_length, &partner, &type, &cflag) != msg_length) {
            printf("2 READ ERROR ... see share_global_double_vec.\n");
            exit(-1);
         }  
         result = md_write((char *) &flag_done_read, 
                sizeof(int), partner, type2, &cflag); 
	 assert( result == 0);
      }

      for (j = 0; j < vec_length; j++) {
        if (double_vec[j] == UNINIT_VEC && tmp[j] != UNINIT_VEC) double_vec[j] = tmp[j];
      }
      } /* end of Proc < partner condition */
      else {
      for (imsg=0; imsg<n_msgs; imsg++) {

         istart = imsg*max_length/sizeof(double);

         if(n_msgs == 1)        msg_length = size;
         else if(imsg<n_msgs-1) msg_length = max_length;
         else                   msg_length = size - (n_msgs-1)*max_length;

         if (md_read((char *) &(tmp[istart]), 
             msg_length, &partner, &type, &cflag) != msg_length) {
            printf("2 READ ERROR ... see share_global_double_vec.\n");
            exit(-1);
         }  
         result = md_write((char *) &flag_done_read, 
                sizeof(int), partner, type2, &cflag); 
	 assert( result == 0);
      }
      for (imsg=0; imsg<n_msgs; imsg++) {

         istart = imsg*max_length/sizeof(double);

         if(n_msgs == 1)        msg_length = size;
         else if(imsg<n_msgs-1) msg_length = max_length;
         else                   msg_length = size - (n_msgs-1)*max_length;

        if (md_write((char *) &(double_vec[istart]), 
            msg_length, partner, type, &cflag)) {
           printf("2 WRITE ERROR ... see share_global_double_vec.\n");
           exit(-1);
        }  
        result = md_read((char *) &flag_done_read, 
              sizeof(int), &partner, &type2, &cflag);
	 assert( (result - sizeof(int)) == 0);
      }


      for (j = 0; j < vec_length; j++) {
        if (double_vec[j] == UNINIT_VEC && tmp[j] != UNINIT_VEC) double_vec[j] = tmp[j];
      }
      }
    }
  }  
  else
    type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = Proc ^ nprocs_small;
  if (Proc & nprocs_small) {
    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(double);

       if(n_msgs == 1)        msg_length = size;
       else if(imsg<n_msgs-1) msg_length = max_length;
       else                   msg_length = size - (n_msgs-1)*max_length;

       if (md_read((char *) &(double_vec[istart]), 
         msg_length, &partner, &type, &cflag) != msg_length) {
         printf("3 READ ERROR ... see share_global_double_vec.\n");
         exit(-1);
       }
       result = md_write((char *) &flag_done_read, 
                sizeof(int), partner, type2, &cflag); 
       assert( result == 0 );
    }
  }  
  else if ((Proc+nprocs_small) < nprocs) {
    for (imsg=0; imsg<n_msgs; imsg++) {

       istart = imsg*max_length/sizeof(double);

       if(n_msgs == 1)        msg_length = size;
       else if(imsg<n_msgs-1) msg_length = max_length;
       else                   msg_length = size - (n_msgs-1)*max_length;

       if (md_write((char *) &(double_vec[istart]), 
         msg_length, partner, type, &cflag)) {
         printf("3 WRITE ERROR ... see share_global_double_vec.\n");
         exit(-1);
       }
       result = md_read((char *) &flag_done_read, 
              sizeof(int), &partner, &type2, &cflag);
       assert( (result - sizeof(int)) == 0);
    }
  }  

  safe_free((void *) &tmp);
}
