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
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_wrapmpi_id =
  "$Id$";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int gl_rbuf = 3;
int gl_sbuf = 3;
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
int the_proc_name = -1;

void get_parallel_info(int *proc, int *nprocs, int *dim)

{

  /* local variables */

  int i;

  MPI_Comm_size(MPI_COMM_WORLD, nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, proc);
  *dim = 0;
the_proc_name = *proc;

} /* get_parallel_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_read(char *buf, int bytes, int *source, int *type, int *flag)

{

  int        err, buffer = 1;
  MPI_Status status;

  if (*type   == -1) *type   = MPI_ANY_TAG;
  if (*source == -1) *source = MPI_ANY_SOURCE;

  if (bytes == 0) {
    err = MPI_Recv(&gl_rbuf, 1, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                   &status);
  }
  else {
    err = MPI_Recv(buf, bytes, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                   &status);
  }

  if (err != 0) (void) fprintf(stderr, "MPI_Recv error = %d\n", err);
  MPI_Get_count(&status,MPI_BYTE,&buffer);
  *source = status.MPI_SOURCE;
  *type   = status.MPI_TAG;
  if (bytes != 0) bytes = buffer;

  return bytes;

} /* md_read */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_write(char *buf, int bytes, int dest, int type, int *flag)

{

  int err;

  if (bytes == 0) {
    err = MPI_Send(&gl_sbuf, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }
  else {
    err = MPI_Send(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }

  if (err != 0) (void) fprintf(stderr, "MPI_Send error = %d\n", err);

  return 0;

} /* md_write */



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_iread(void *buf, int bytes, int *source, int *type,
                  MPI_Request *request)


/*******************************************************************************

  Machine dependent wrapped message-reading communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  source:          Source processor number.

  type:            Message type

*******************************************************************************/

{

  int err = 0;

  if (*type   == -1) *type   = MPI_ANY_TAG;
  if (*source == -1) *source = MPI_ANY_SOURCE;

  if (bytes == 0) {
    err = MPI_Irecv(&gl_rbuf, 1, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                    request);
  }
  else {
    err = MPI_Irecv(buf, bytes, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                    request);
  }

  return err;

} /* md_wrap_iread */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_write(void *buf, int bytes, int dest, int type, int *flag)

/*******************************************************************************

  Machine dependent wrapped message-sending communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  int err = 0;

  if (bytes == 0) {
    err = MPI_Send(&gl_sbuf, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }
  else {
    err = MPI_Send(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }

  return err;

} /* md_wrap_write */



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_wait(void *buf, int bytes, int *source, int *type, int *flag,
                 MPI_Request *request)

/*******************************************************************************

  Machine dependent wrapped message-wait communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.
  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  int        err, count;
  MPI_Status status;

  if ( MPI_Wait(request, &status) ) {
    (void) fprintf(stderr, "MPI_Wait error\n");
    exit(-1);
  }

  MPI_Get_count(&status, MPI_BYTE, &count);
  *source = status.MPI_SOURCE;
  *type   = status.MPI_TAG;

  /* return the count, which is in bytes */

  return count;

} /* md_wrap_wait */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_iwrite(void *buf, int bytes, int dest, int type, int *flag,
                  MPI_Request *request)

/*******************************************************************************

  Machine dependent wrapped message-sending (nonblocking) communication 
  routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  int err = 0;

  if (bytes == 0) {
    err = MPI_Isend(&gl_sbuf, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD,
                  request);
  }
  else {
    err = MPI_Isend(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD,
                  request);
  }

  return err;

} /* md_wrap_write */


/********************************************************************/
/*     NEW WRAPPERS to handle MPI Communicators                     */
/********************************************************************/

void parallel_info(int *proc,int *nprocs,int *dim, MPI_Comm comm)
{

  /* local variables */

  int i;

  MPI_Comm_size(comm, nprocs);
  MPI_Comm_rank(comm, proc);
  *dim = 0;
the_proc_name = *proc;

} /* get_parallel_info */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_mpi_iread(void *buf, int bytes, int *source, int *type,
                  MPI_Request *request, int *icomm)


/*******************************************************************************

  Machine dependent wrapped message-reading communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  source:          Source processor number.

  type:            Message type

  icomm:           MPI Communicator
*******************************************************************************/

{

  int err = 0;
  MPI_Comm *comm;

  comm = (MPI_Comm *) icomm;

  if (*type   == -1) *type   = MPI_ANY_TAG;
  if (*source == -1) *source = MPI_ANY_SOURCE;

  if (bytes == 0) {
    err = MPI_Irecv(&gl_rbuf, 1, MPI_BYTE, *source, *type, *comm,
                    request);
  }
  else {
    err = MPI_Irecv(buf, bytes, MPI_BYTE, *source, *type, *comm,
                    request);
  }

  return err;

} /* md_mpi_iread */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_mpi_write(void *buf, int bytes, int dest, int type, int *flag,
                  int *icomm)

/*******************************************************************************

  Machine dependent wrapped message-sending communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  dest:            Destination processor number.

  type:            Message type

  flag:

  icomm:           MPI Communicator

*******************************************************************************/

{

  int err = 0;
  MPI_Comm *comm;

  comm = (MPI_Comm *) icomm;

  if (bytes == 0) {
    err = MPI_Send(&gl_sbuf, 1, MPI_BYTE, dest, type, *comm);
  }
  else {
    err = MPI_Send(buf, bytes, MPI_BYTE, dest, type, *comm);
  }

  return err;

} /* md_wrap_write */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_mpi_wait(void *buf, int bytes, int *source, int *type, int *flag,
                 MPI_Request *request, int *icomm)

/*******************************************************************************

  Machine dependent wrapped message-wait communication routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.
  dest:            Destination processor number.

  type:            Message type

  flag:

  icomm:           MPI Communicator

*******************************************************************************/

{

  int        err, count;
  MPI_Status status;

  if ( MPI_Wait(request, &status) ) {
    (void) fprintf(stderr, "MPI_Wait error\n");
    exit(-1);
  }

  MPI_Get_count(&status, MPI_BYTE, &count);
  *source = status.MPI_SOURCE;
  *type   = status.MPI_TAG;

  /* return the count, which is in bytes */

  return count;

} /* md_mpi_wait */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_mpi_iwrite(void *buf, int bytes, int dest, int type, int *flag,
                  MPI_Request *request, int *icomm)

/*******************************************************************************

  Machine dependent wrapped message-sending (nonblocking) communication
  routine for MPI.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  dest:            Destination processor number.

  type:            Message type

  flag:

  icomm:           MPI Communicator

*******************************************************************************/
{

  int err = 0;
  MPI_Comm *comm;

  comm = (MPI_Comm *) icomm ;
  if (bytes == 0)
    err = MPI_Isend(&gl_sbuf, 1, MPI_BYTE, dest, type, *comm, request);
  else
    err = MPI_Isend(buf, bytes, MPI_BYTE, dest, type, *comm, request);

  return err;

} /* md_mpi_write */
