/* This file was automatically generated.  Do not edit! */
int md_mpi_iwrite(void *buf,int bytes,int dest,int type,int *flag,MPI_Request *request,int *icomm);
int md_mpi_wait(void *buf,int bytes,int *source,int *type,int *flag,MPI_Request *request,int *icomm);
int md_mpi_write(void *buf,int bytes,int dest,int type,int *flag,int *icomm);
int md_mpi_iread(void *buf,int bytes,int *source,int *type,MPI_Request *request,int *icomm);
void parallel_info(int *proc,int *nprocs,int *dim,MPI_Comm comm);
int md_wrap_iwrite(void *buf,int bytes,int dest,int type,int *flag,MPI_Request *request);
int md_wrap_wait(void *buf,int bytes,int *source,int *type,int *flag,MPI_Request *request);
int md_wrap_write(void *buf,int bytes,int dest,int type,int *flag);
int md_wrap_iread(void *buf,int bytes,int *source,int *type,MPI_Request *request);
int md_write(char *buf,int bytes,int dest,int type,int *flag);
int md_read(char *buf,int bytes,int *source,int *type,int *flag);
void get_parallel_info(int *proc,int *nprocs,int *dim);
extern int the_proc_name;
extern int gl_sbuf;
extern int gl_rbuf;
