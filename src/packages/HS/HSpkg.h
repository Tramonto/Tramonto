/* function definitions - only used in fill of the hard sphere (fundamental measures theory)  equations  */

external
struct RB_Struct *dphi_drb=NULL;

external
static struct RB_Struct d2phi_drb2_delta_rb_FMT1(int, int, double **,double,
                                            int *,double *,double,double,
                                            double);
external
static struct RB_Struct d2phi_drb2_theta_rb_FMT1(int, int, double **,double,int *);

external
static struct RB_Struct d2phi_drb2_delta_rb_FMT2(int, int, double **,double,
                                      int *,double *,double, double,double);

external
static struct RB_Struct d2phi_drb2_theta_rb_FMT2(int, int, double **,double,int *);

external
static struct RB_Struct d2phi_drb2_delta_rb_FMT3(int, int, double **,double,
                                      int *,double *,double, double,double);
external
static struct RB_Struct d2phi_drb2_theta_rb_FMT3(int, int, double **,double,int *);

extern void FMT1_1stderiv(double *,double,double,double *,double *);
extern void FMT2_1stderiv(double *,double,double,double *,double *);
extern void FMT3_1stderiv(double *,double,double,double *,double *);
extern void calc_FMT_derivatives(void(*fp_FMTderiv)(double *,double,double,double *,double *),
                                 int,double **,struct RB_Struct *);

