/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_quadrature.c
 *
 *  This file contains the quadrature schemes for the non-local
 *  interactions.  
 *
 */

#include "dft_quadrature.h"

/****************************************************************************/

void set_gauss_quad(int ngp, double *gp, double *gw)
{
  if (ngp == 40) {
    gp[ 0] = ( 0.038772417506050821933 + 1.0)/2.0;  gw[ 0] = 0.077505947978424811264/2.0;
    gp[ 1] = ( 0.116084070675255208483 + 1.0)/2.0;  gw[ 1] = 0.077039818164247965588/2.0;
    gp[ 2] = ( 0.192697580701371099716 + 1.0)/2.0;  gw[ 2] = 0.076110361900626242372/2.0;
    gp[ 3] = ( 0.268152185007253681141 + 1.0)/2.0;  gw[ 3] = 0.074723169057968264200/2.0;
    gp[ 4] = ( 0.341994090825758473007 + 1.0)/2.0;  gw[ 4] = 0.072886582395804059061/2.0;
    gp[ 5] = ( 0.413779204371605001525 + 1.0)/2.0;  gw[ 5] = 0.070611647391286779695/2.0;
    gp[ 6] = ( 0.483075801686178712909 + 1.0)/2.0;  gw[ 6] = 0.067912045815233903826/2.0;
    gp[ 7] = ( 0.549467125095128202076 + 1.0)/2.0;  gw[ 7] = 0.064804013456601038075/2.0;
    gp[ 8] = ( 0.612553889667980237953 + 1.0)/2.0;  gw[ 8] = 0.061306242492928939167/2.0;
    gp[ 9] = ( 0.671956684614179548379 + 1.0)/2.0;  gw[ 9] = 0.057439769099391551367/2.0;
    gp[10] = ( 0.727318255189927103281 + 1.0)/2.0;  gw[10] = 0.053227846983936824355/2.0;
    gp[11] = ( 0.778305651426519387695 + 1.0)/2.0;  gw[11] = 0.048695807635072232061/2.0;
    gp[12] = ( 0.824612230833311663196 + 1.0)/2.0;  gw[12] = 0.043870908185673271992/2.0;
    gp[13] = ( 0.865959503212259503821 + 1.0)/2.0;  gw[13] = 0.038782167974472017640/2.0;
    gp[14] = ( 0.902098806968874296728 + 1.0)/2.0;  gw[14] = 0.033460195282547847393/2.0;
    gp[15] = ( 0.932812808278676533361 + 1.0)/2.0;  gw[15] = 0.027937006980023401098/2.0;
    gp[16] = ( 0.957916819213791655805 + 1.0)/2.0;  gw[16] = 0.022245849194166957262/2.0;
    gp[17] = ( 0.977259949983774262663 + 1.0)/2.0;  gw[17] = 0.016421058381907888713/2.0;
    gp[18] = ( 0.990726238699457006453 + 1.0)/2.0;  gw[18] = 0.010498284531152813615/2.0;
    gp[19] = ( 0.998237709710559200350 + 1.0)/2.0;  gw[19] = 0.004521277098533191258/2.0;
    gp[20] = (-0.038772417506050821933 + 1.0)/2.0;  gw[20] = gw[0];
    gp[21] = (-0.116084070675255208483 + 1.0)/2.0;  gw[21] = gw[1];
    gp[22] = (-0.192697580701371099716 + 1.0)/2.0;  gw[22] = gw[2];
    gp[23] = (-0.268152185007253681141 + 1.0)/2.0;  gw[23] = gw[3];
    gp[24] = (-0.341994090825758473007 + 1.0)/2.0;  gw[24] = gw[4];
    gp[25] = (-0.413779204371605001525 + 1.0)/2.0;  gw[25] = gw[5];
    gp[26] = (-0.483075801686178712909 + 1.0)/2.0;  gw[26] = gw[6];
    gp[27] = (-0.549467125095128202076 + 1.0)/2.0;  gw[27] = gw[7];
    gp[28] = (-0.612553889667980237953 + 1.0)/2.0;  gw[28] = gw[8];
    gp[29] = (-0.671956684614179548379 + 1.0)/2.0;  gw[29] = gw[9];
    gp[30] = (-0.727318255189927103281 + 1.0)/2.0;  gw[30] = gw[10];
    gp[31] = (-0.778305651426519387695 + 1.0)/2.0;  gw[31] = gw[11];
    gp[32] = (-0.824612230833311663196 + 1.0)/2.0;  gw[32] = gw[12];
    gp[33] = (-0.865959503212259503821 + 1.0)/2.0;  gw[33] = gw[13];
    gp[34] = (-0.902098806968874296728 + 1.0)/2.0;  gw[34] = gw[14];
    gp[35] = (-0.932812808278676533361 + 1.0)/2.0;  gw[35] = gw[15];
    gp[36] = (-0.957916819213791655805 + 1.0)/2.0;  gw[36] = gw[16];
    gp[37] = (-0.977259949983774262663 + 1.0)/2.0;  gw[37] = gw[17];
    gp[38] = (-0.990726238699457006453 + 1.0)/2.0;  gw[38] = gw[18];
    gp[39] = (-0.998237709710559200350 + 1.0)/2.0;  gw[39] = gw[19];
  }
  else if (ngp == 20) {
    gp[ 0] = ( 0.076526521133497333755 +1.0 )/2.0;  gw[ 0] = 0.152753387130725850698/2.0;
    gp[ 1] = ( 0.227785851141645078080 +1.0 )/2.0;  gw[ 1] = 0.149172986472603746788/2.0;
    gp[ 2] = ( 0.373706088715419560673 + 1.0)/2.0;  gw[ 2] = 0.142096109318382051329/2.0;
    gp[ 3] = ( 0.510867001950827098004 + 1.0)/2.0;  gw[ 3] = 0.131688638449176626898/2.0;
    gp[ 4] = ( 0.636053680726515025453 + 1.0)/2.0;  gw[ 4] = 0.118194531961518417312/2.0;
    gp[ 5] = ( 0.746331906460150792614 + 1.0)/2.0;  gw[ 5] = 0.101930119817240435037/2.0;
    gp[ 6] = (0.839116971822218823395 + 1.0)/2.0;  gw[ 6] = 0.083276741576704748725/2.0;
    gp[ 7] = (0.912234428251325905868 + 1.0)/2.0;  gw[ 7] = 0.062672048334109063570/2.0;
    gp[ 8] = (0.963971927277913791268 + 1.0)/2.0;  gw[ 8] = 0.040601429800386941331/2.0;
    gp[ 9] = (0.993128599185094924786 + 1.0)/2.0;  gw[ 9] = 0.017614007139152118312/2.0;
    gp[10] = (-0.076526521133497333755 +1.0 )/2.0;  gw[10] = gw[0];
    gp[11] = (-0.227785851141645078080 +1.0 )/2.0;  gw[11] = gw[1];
    gp[12] = (-0.373706088715419560673 + 1.0)/2.0;  gw[12] = gw[2];
    gp[13] = (-0.510867001950827098004 + 1.0)/2.0;  gw[13] = gw[3];
    gp[14] = (-0.636053680726515025453 + 1.0)/2.0;  gw[14] = gw[4];
    gp[15] = (-0.746331906460150792614 + 1.0)/2.0;  gw[15] = gw[5];
    gp[16] = (-0.839116971822218823395 + 1.0)/2.0;  gw[16] = gw[6];
    gp[17] = (-0.912234428251325905868 + 1.0)/2.0;  gw[17] = gw[7];
    gp[18] = (-0.963971927277913791268 + 1.0)/2.0;  gw[18] = gw[8];
    gp[19] = (-0.993128599185094924786 + 1.0)/2.0;  gw[19] = gw[9];
  }

  else if (ngp == 12) {
    gp[ 0] = ( 0.125233408511469 + 1.0)/2.0;  gw[ 0] = 0.249147045813403/2.0;
    gp[ 1] = ( 0.367831498998180 + 1.0)/2.0;  gw[ 1] = 0.233492536538355/2.0;
    gp[ 2] = ( 0.587317954286617 + 1.0)/2.0;  gw[ 2] = 0.203167426723066/2.0;
    gp[ 3] = ( 0.769902674194305 + 1.0)/2.0;  gw[ 3] = 0.160078328543346/2.0;
    gp[ 4] = ( 0.904117256370475 + 1.0)/2.0;  gw[ 4] = 0.106939325995318/2.0;
    gp[ 5] = ( 0.981560634246719 + 1.0)/2.0;  gw[ 5] = 0.047175336386512/2.0;
    gp[ 6] = (-0.125233408511469 + 1.0)/2.0;  gw[ 6] = 0.249147045813403/2.0;
    gp[ 7] = (-0.367831498998180 + 1.0)/2.0;  gw[ 7] = 0.233492536538355/2.0;
    gp[ 8] = (-0.587317954286617 + 1.0)/2.0;  gw[ 8] = 0.203167426723066/2.0;
    gp[ 9] = (-0.769902674194305 + 1.0)/2.0;  gw[ 9] = 0.160078328543346/2.0;
    gp[10] = (-0.904117256370475 + 1.0)/2.0;  gw[10] = 0.106939325995318/2.0;
    gp[11] = (-0.981560634246719 + 1.0)/2.0;  gw[11] = 0.047175336386512/2.0;
  }
  else if (ngp == 3) {
    gp[ 0] = ( 0.7745966692      + 1.0)/2.0;  gw[ 0] = 0.555555555555556/2.0;
    gp[ 1] = ( 0.0000000         + 1.0)/2.0;  gw[ 1] = 0.888888888888888/2.0;
    gp[ 2] = (-0.7745966692      + 1.0)/2.0;  gw[ 2] = 0.555555555555556/2.0;
  }
  else if (ngp == 1) {
    gp[ 0] = ( 0.0000000         + 1.0)/2.0;  gw[ 0] = 2.000000000000000/2.0;
  }
  else if (ngp == 5) {
    gp[ 0] = ( 0.906179845936884 + 1.0)/2.0;  gw[ 0] = 0.236926885056189/2.0;
    gp[ 1] = ( 0.538469310105683 + 1.0)/2.0;  gw[ 1] = 0.478628670499366/2.0;
    gp[ 2] = ( 0.00000           + 1.0)/2.0;  gw[ 2] = 0.568888888888890/2.0;
    gp[ 3] = (-0.538469310105683 + 1.0)/2.0;  gw[ 3] = 0.478628670499366/2.0;
    gp[ 4] = (-0.906179845936884 + 1.0)/2.0;  gw[ 4] = 0.236926885056189/2.0;
  }
  else if (ngp == 6) {
    gp[ 0] = ( 0.238619186083197 + 1.0)/2.0;  gw[ 0] = 0.467913934572691/2.0;
    gp[ 1] = ( 0.661209386466265 + 1.0)/2.0;  gw[ 1] = 0.360761573048139/2.0;
    gp[ 2] = ( 0.932469514203152 + 1.0)/2.0;  gw[ 2] = 0.171324492379170/2.0;
    gp[ 3] = (-0.238619186083197 + 1.0)/2.0;  gw[ 3] = 0.467913934572691/2.0;
    gp[ 4] = (-0.661209386466265 + 1.0)/2.0;  gw[ 4] = 0.360761573048139/2.0;
    gp[ 5] = (-0.932469514203152 + 1.0)/2.0;  gw[ 5] = 0.171324492379170/2.0;
  }
  else if (ngp == 9) {
    gp[ 0] = ( 0.324253423403809 + 1.0)/2.0;  gw[ 0] = 0.312347077040004/2.0;
    gp[ 1] = ( 0.613371432700590 + 1.0)/2.0;  gw[ 1] = 0.260610696402935/2.0;
    gp[ 2] = ( 0.836031107326636 + 1.0)/2.0;  gw[ 2] = 0.180648160694857/2.0;
    gp[ 3] = ( 0.968160239507626 + 1.0)/2.0;  gw[ 3] = 0.081274388361574/2.0;
    gp[ 4] = ( 0.000000000000000 + 1.0)/2.0;  gw[ 4] = 0.330239355001260/2.0;
    gp[ 5] = (-0.324253423403809 + 1.0)/2.0;  gw[ 5] = 0.312347077040004/2.0;
    gp[ 6] = (-0.613371432700590 + 1.0)/2.0;  gw[ 6] = 0.260610696402935/2.0;
    gp[ 7] = (-0.836031107326636 + 1.0)/2.0;  gw[ 7] = 0.180648160694857/2.0;
    gp[ 8] = (-0.968160239507626 + 1.0)/2.0;  gw[ 8] = 0.081274388361574/2.0;
  }
  else {
    if (Iwrite_screen != SCREEN_NONE) printf("Requested Number of Gauss points not allowed (%d)\n",ngp);
    exit(-1);
  }
  return;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
int  get_integration_pts(int isten, int izone,
                 double ***point_ptr, double **wt_ptr)
/*
 * This routine returns the number of integration points in the
 * or quadrature scheme, the points normalized to a Radius of 1.0,
 * and weights which sum to one.  Note that these qudratures are NOT 
 * used in the code as of January 2007.
 */
{
  char *yo = "get_num_integration_pts";
  int num_int_pts, num_surf_pts, i, j, k;
  double gauss_pt[7], gauss_wt[7], **point_tmp, *wt_tmp;

  /********************** BEGIN EXECUTION ************************************/

  if (isten==DELTA_FN_R || isten==DELTA_FN_BOND) {

    switch (Sten_Choice_S[isten][izone]) {

      case TETRAHEDRON:
         num_int_pts = 4;
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         delta_tetrahedron(*point_ptr, *wt_ptr);
         return (num_int_pts);

      case OCTAHEDRON:
         num_int_pts = 6;
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         delta_octahedron(*point_ptr, *wt_ptr);
         return (num_int_pts);

      case CUBE:
         num_int_pts = 8;
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         delta_cube(*point_ptr, *wt_ptr);
         return (num_int_pts);

      case ICOSAHEDRON:
         num_int_pts = 12;
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         delta_icosahedron(*point_ptr, *wt_ptr);
         return (num_int_pts);

      case DODECAHEDRON:
         num_int_pts = 20;
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         delta_dodecahedron(*point_ptr, *wt_ptr);
         return (num_int_pts);

      case SEVENTY_TWO:
         num_int_pts = 72;
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         delta_seventy_two(*point_ptr, *wt_ptr);
         return (num_int_pts);

      case ONE_D_TWELVE:
         if (Ndim!=1) {
            if (Iwrite_screen != SCREEN_NONE) printf("ERROR: THIS STENCIL ONLY FOR 1D problems\n"); 
            exit(-1);
         }
         num_int_pts = 12;
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         delta_one_D_twelve(*point_ptr, *wt_ptr);
         return (num_int_pts);

      case MIDPOINT_RULE:
         if (Ndim==1) 
            num_int_pts = MPsten_Npts_R[izone];
         else if (Ndim==2) 
            num_int_pts = MPsten_Npts_R[izone] * MPsten_Npts_arc[izone];
         else {
            num_int_pts = MPsten_Npts_R[izone] * MPsten_Npts_arc[izone]
                          *MPsten_Npts_phi[izone];
/*            if (Iwrite_screen != SCREEN_NONE) printf("ERROR: THIS STENCIL ONLY FOR 1D or 2D problems\n");
            exit(-1);*/
         }
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         delta_midpoint(*point_ptr, *wt_ptr, izone);
         return (num_int_pts);

      default:
         if (Iwrite_screen != SCREEN_NONE) {
            printf("%s: Unknown or uninitialized choice of integration\n",yo);
            printf("\t\tscheme for isten: %d, izone: %d, Sten_Choice_S: %d\n",
                                isten, izone, Sten_Choice_S[isten][izone]);
         }
         exit(-1); 
         return(-1);break;
    }
  }

  else if (isten==THETA_FN_R || isten==THETA_PAIRPOT_RCUT || 
           isten==THETA_CR_RPM_MSA || isten==THETA_CR_GENERAL_MSA) {

    switch (Sten_Choice_S[isten][izone]) {

      case TETRAHEDRON:

         /* 
          * Theta_Fn stencil is a product of surface and radial stencils
          * Sten_Choice_R should hold the number of radial Gauss point (1..7)
             (Surface*Radial) 
          */

         num_surf_pts = 4;
         num_int_pts = num_surf_pts * Sten_Choice_R[isten][izone];

        /* Allocate arrays called point and wt in dft_stencil.c */

        *point_ptr =(double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    =(double *)  array_alloc(1, num_int_pts,    sizeof(double));

        /* Get weights and radial positions for radial quadrature */

         get_radial_quadrature(gauss_pt, gauss_wt, Sten_Choice_R[isten][izone]);

        /*
         * For each radial Gauss point, get the surface quadrature point,
         * and then scale them by the radius and weight of the radial scheme
         */

         for (i=0; i<Sten_Choice_R[isten][izone]; i++) {

           point_tmp = &((*point_ptr)[i*num_surf_pts]);
           wt_tmp    = &((*wt_ptr   )[i*num_surf_pts]);

           delta_tetrahedron(point_tmp, wt_tmp);

           for (j=0; j<num_surf_pts; j++) {
             wt_tmp[j] *= gauss_wt[i];
             for (k=0; k<Ndim; k++) point_tmp[j][k] *= gauss_pt[i];
           }
         }
         return (num_int_pts);

      case OCTAHEDRON: /* See comments in tetrahedon section */
         num_surf_pts = 6;
         num_int_pts = num_surf_pts * Sten_Choice_R[isten][izone]; /* Surface*Radial */
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         get_radial_quadrature(gauss_pt, gauss_wt, Sten_Choice_R[isten][izone]);
         for (i=0; i<Sten_Choice_R[isten][izone]; i++) {
           point_tmp = &((*point_ptr)[i*num_surf_pts]);
           wt_tmp    = &((*wt_ptr   )[i*num_surf_pts]);
           delta_octahedron(point_tmp, wt_tmp);
           for (j=0; j<num_surf_pts; j++) {
             wt_tmp[j] *= gauss_wt[i];
             for (k=0; k<Ndim; k++) point_tmp[j][k] *= gauss_pt[i];
           }
         }
         return (num_int_pts);

      case CUBE: /* See comments in tetrahedon section */
         num_surf_pts = 8;
         num_int_pts = num_surf_pts * Sten_Choice_R[isten][izone]; /* Surface*Radial */
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         get_radial_quadrature(gauss_pt, gauss_wt, Sten_Choice_R[isten][izone]);
         for (i=0; i<Sten_Choice_R[isten][izone]; i++) {
           point_tmp = &((*point_ptr)[i*num_surf_pts]);
           wt_tmp    = &((*wt_ptr   )[i*num_surf_pts]);
           delta_cube(point_tmp, wt_tmp);
           for (j=0; j<num_surf_pts; j++) {
             wt_tmp[j] *= gauss_wt[i];
             for (k=0; k<Ndim; k++) point_tmp[j][k] *= gauss_pt[i];
           }
         }
         return (num_int_pts);

      case ICOSAHEDRON: /* See comments in tetrahedon section */
         num_surf_pts = 12;
         num_int_pts = num_surf_pts * Sten_Choice_R[isten][izone]; /* Surface*Radial */
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         get_radial_quadrature(gauss_pt, gauss_wt, Sten_Choice_R[isten][izone]);
         for (i=0; i<Sten_Choice_R[isten][izone]; i++) {
           point_tmp = &((*point_ptr)[i*num_surf_pts]);
           wt_tmp    = &((*wt_ptr   )[i*num_surf_pts]);
           delta_icosahedron(point_tmp, wt_tmp);
           for (j=0; j<num_surf_pts; j++) {
             wt_tmp[j] *= gauss_wt[i];
             for (k=0; k<Ndim; k++) point_tmp[j][k] *= gauss_pt[i];
           }
         }
         return (num_int_pts);

      case DODECAHEDRON: /* See comments in tetrahedon section */
         num_surf_pts = 20;
         num_int_pts = num_surf_pts * Sten_Choice_R[isten][izone]; /* Surface*Radial */
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         get_radial_quadrature(gauss_pt, gauss_wt, Sten_Choice_R[isten][izone]);
         for (i=0; i<Sten_Choice_R[isten][izone]; i++) {
           point_tmp = &((*point_ptr)[i*num_surf_pts]);
           wt_tmp    = &((*wt_ptr   )[i*num_surf_pts]);
           delta_dodecahedron(point_tmp, wt_tmp);
           for (j=0; j<num_surf_pts; j++) {
             wt_tmp[j] *= gauss_wt[i];
             for (k=0; k<Ndim; k++) point_tmp[j][k] *= gauss_pt[i];
           }
         }
         return (num_int_pts);

      case SEVENTY_TWO: /* See comments in tetrahedon section */
         num_surf_pts = 72;
         num_int_pts = num_surf_pts * Sten_Choice_R[isten][izone]; /* Surface*Radial */
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         get_radial_quadrature(gauss_pt, gauss_wt, Sten_Choice_R[isten][izone]);
         for (i=0; i<Sten_Choice_R[isten][izone]; i++) {
           point_tmp = &((*point_ptr)[i*num_surf_pts]);
           wt_tmp    = &((*wt_ptr   )[i*num_surf_pts]);
           delta_seventy_two(point_tmp, wt_tmp);
           for (j=0; j<num_surf_pts; j++) {
             wt_tmp[j] *= gauss_wt[i];
             for (k=0; k<Ndim; k++) point_tmp[j][k] *= gauss_pt[i];
           }
         }
         return (num_int_pts);

      case ONE_D_TWELVE: /* See comments in tetrahedon section */
         if (Ndim!=1) {
            if (Iwrite_screen != SCREEN_NONE) printf("ERROR: STENCIL ONE_D_TWELVE ONLY FOR 1D problems\n"); 
            exit(-1);
         }
         num_surf_pts = 12;
         num_int_pts = num_surf_pts * Sten_Choice_R[isten][izone]; /* Surface*Radial */
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         get_radial_quadrature(gauss_pt, gauss_wt, Sten_Choice_R[isten][izone]);
         for (i=0; i<Sten_Choice_R[isten][izone]; i++) {
           point_tmp = &((*point_ptr)[i*num_surf_pts]);
           wt_tmp    = &((*wt_ptr   )[i*num_surf_pts]);
           delta_one_D_twelve(point_tmp, wt_tmp);
           for (j=0; j<num_surf_pts; j++) {
             wt_tmp[j] *= gauss_wt[i];
             for (k=0; k<Ndim; k++) point_tmp[j][k] *= gauss_pt[i];
           }
         }
         return (num_int_pts);

      case MIDPOINT_RULE:

         if (isten==THETA_FN_R) { 
           if (Ndim==1) 
              num_int_pts = MPsten_Npts_R[izone];
           else if (Ndim==2) 
              num_int_pts = MPsten_Npts_R[izone] * MPsten_Npts_arc[izone];
           else  
              num_int_pts = MPsten_Npts_R[izone] * MPsten_Npts_arc[izone] 
                                                 * MPsten_Npts_phi[izone];
         }
         else if (isten==THETA_PAIRPOT_RCUT || isten==THETA_CR_RPM_MSA || isten==THETA_CR_GENERAL_MSA)
           num_int_pts = MPsten_Npts_R[izone] * MPsten_Npts_arc[izone] 
                                              * MPsten_Npts_phi[izone];
        *point_ptr = (double **) array_alloc(2, num_int_pts, 3, sizeof(double));
        *wt_ptr    = (double *)  array_alloc(1, num_int_pts,    sizeof(double));
         if (isten==THETA_FN_R) 
           theta_midpoint(*point_ptr, *wt_ptr, izone, Ndim);
         else if (isten==THETA_PAIRPOT_RCUT || isten==THETA_CR_RPM_MSA || isten==THETA_CR_GENERAL_MSA )
           theta_midpoint(*point_ptr, *wt_ptr, izone, 3);
         return (num_int_pts);

      default:
         if (Iwrite_screen != SCREEN_NONE) printf("%s: Unknown or uninitialized choice of integration\n",yo);
         if (Iwrite_screen != SCREEN_NONE) printf("\t\tscheme for isten: %d, izone: %d, Sten_Choice_S: %d\n",
                                isten, izone, Sten_Choice_S[isten][izone]);
         exit(-1); 
         return(-1);break;
    }
  }

}

/****************************************************************************/

void delta_tetrahedron(double **point, double *wt)
{
/* one over the square root of three */
#define SQRT1_3  5.7735026918962584e-01

  /********************** BEGIN EXECUTION ************************************/

  point[0][0]= SQRT1_3; point[0][1]= SQRT1_3; point[0][2]= SQRT1_3; wt[0]=0.25;
  point[1][0]= SQRT1_3; point[1][1]=-SQRT1_3; point[1][2]=-SQRT1_3; wt[1]=0.25;
  point[2][0]=-SQRT1_3; point[2][1]= SQRT1_3; point[2][2]=-SQRT1_3; wt[2]=0.25;
  point[3][0]=-SQRT1_3; point[3][1]=-SQRT1_3; point[3][2]= SQRT1_3; wt[3]=0.25;
}
/****************************************************************************/

void delta_octahedron(double **point, double *wt)
{
#define SIXTH 1.6666666666666666e-01
  /********************** BEGIN EXECUTION ************************************/

  point[0][0]= 1.0;  point[0][1]= 0.0;  point[0][2]= 0.0; wt[0]=SIXTH;
  point[1][0]=-1.0;  point[1][1]= 0.0;  point[1][2]= 0.0; wt[1]=SIXTH;
  point[2][0]= 0.0;  point[2][1]= 1.0;  point[2][2]= 0.0; wt[2]=SIXTH;
  point[3][0]= 0.0;  point[3][1]=-1.0;  point[3][2]= 0.0; wt[3]=SIXTH;
  point[4][0]= 0.0;  point[4][1]= 0.0;  point[4][2]= 1.0; wt[4]=SIXTH;
  point[5][0]= 0.0;  point[5][1]= 0.0;  point[5][2]=-1.0; wt[5]=SIXTH;
}
/****************************************************************************/

void delta_cube(double **point, double *wt)
{
  /********************** BEGIN EXECUTION ************************************/

  point[0][0]= SQRT1_3; point[0][1]= SQRT1_3; point[0][2]= SQRT1_3; wt[0]=0.125;
  point[1][0]= SQRT1_3; point[1][1]= SQRT1_3; point[1][2]=-SQRT1_3; wt[1]=0.125;
  point[2][0]= SQRT1_3; point[2][1]=-SQRT1_3; point[2][2]= SQRT1_3; wt[2]=0.125;
  point[3][0]= SQRT1_3; point[3][1]=-SQRT1_3; point[3][2]=-SQRT1_3; wt[3]=0.125;
  point[4][0]=-SQRT1_3; point[4][1]= SQRT1_3; point[4][2]= SQRT1_3; wt[4]=0.125;
  point[5][0]=-SQRT1_3; point[5][1]= SQRT1_3; point[5][2]=-SQRT1_3; wt[5]=0.125;
  point[6][0]=-SQRT1_3; point[6][1]=-SQRT1_3; point[6][2]= SQRT1_3; wt[6]=0.125;
  point[7][0]=-SQRT1_3; point[7][1]=-SQRT1_3; point[7][2]=-SQRT1_3; wt[7]=0.125;
}
/****************************************************************************/

void delta_icosahedron(double **point, double *wt)
{
#define TWELFTH 8.3333333333333329e-02
/*
 * The vertices are (0, 1, GM) where GM= (1+sqrt(5)/2. 
 * Normalized, these are 0 and:
 */
#define NRM_ONE 5.2573111211913359e-01
#define NRM_GM  8.5065080835203999e-01

  /********************** BEGIN EXECUTION ************************************/

  point[0][0]= NRM_GM; point[0][1]= NRM_ONE;point[0][2]= 0.0;     wt[0]=TWELFTH;
  point[1][0]=-NRM_GM; point[1][1]= NRM_ONE;point[1][2]= 0.0;     wt[1]=TWELFTH;
  point[2][0]= NRM_GM; point[2][1]=-NRM_ONE;point[2][2]= 0.0;     wt[2]=TWELFTH;
  point[3][0]=-NRM_GM; point[3][1]=-NRM_ONE;point[3][2]= 0.0;     wt[3]=TWELFTH;
  point[4][0]= NRM_ONE;point[4][1]= 0.0;    point[4][2]= NRM_GM;  wt[4]=TWELFTH;
  point[5][0]= NRM_ONE;point[5][1]= 0.0;    point[5][2]=-NRM_GM;  wt[5]=TWELFTH;
  point[6][0]=-NRM_ONE;point[6][1]= 0.0;    point[6][2]= NRM_GM;  wt[6]=TWELFTH;
  point[7][0]=-NRM_ONE;point[7][1]= 0.0;    point[7][2]=-NRM_GM;  wt[7]=TWELFTH;
  point[8][0]= 0.0;    point[8][1]= NRM_GM; point[8][2]= NRM_ONE; wt[8]=TWELFTH;
  point[9][0]= 0.0;    point[9][1]=-NRM_GM; point[9][2]= NRM_ONE; wt[9]=TWELFTH;
  point[10][0]= 0.0;  point[10][1]= NRM_GM;point[10][2]=-NRM_ONE;wt[10]=TWELFTH;
  point[11][0]= 0.0;  point[11][1]=-NRM_GM;point[11][2]=-NRM_ONE;wt[11]=TWELFTH;
}
/****************************************************************************/

void delta_dodecahedron(double **point, double *wt)
{
/*
 * Some vertices are (1, 1, 1), which normalized gives SQRT1_3 defined above
 * The vertices are (0, GM-1, GM) where GM= (1+sqrt(5)/2. 
 * Normalized, these are 0 and:
 */
#define NRM_GM1 3.5682208977308999e-01
#define NRM_GMB 9.3417235896271578e-01
  /********************** BEGIN EXECUTION ************************************/

  point[ 0][0]= SQRT1_3;point[ 0][1]= SQRT1_3;point[ 0][2]= SQRT1_3;wt[ 0]=0.05;
  point[ 1][0]= SQRT1_3;point[ 1][1]= SQRT1_3;point[ 1][2]=-SQRT1_3;wt[ 1]=0.05;
  point[ 2][0]= SQRT1_3;point[ 2][1]=-SQRT1_3;point[ 2][2]= SQRT1_3;wt[ 2]=0.05;
  point[ 3][0]= SQRT1_3;point[ 3][1]=-SQRT1_3;point[ 3][2]=-SQRT1_3;wt[ 3]=0.05;
  point[ 4][0]=-SQRT1_3;point[ 4][1]= SQRT1_3;point[ 4][2]= SQRT1_3;wt[ 4]=0.05;
  point[ 5][0]=-SQRT1_3;point[ 5][1]= SQRT1_3;point[ 5][2]=-SQRT1_3;wt[ 5]=0.05;
  point[ 6][0]=-SQRT1_3;point[ 6][1]=-SQRT1_3;point[ 6][2]= SQRT1_3;wt[ 6]=0.05;
  point[ 7][0]=-SQRT1_3;point[ 7][1]=-SQRT1_3;point[ 7][2]=-SQRT1_3;wt[ 7]=0.05;
  point[ 8][0]= NRM_GM1;point[ 8][1]= NRM_GMB;point[ 8][2]= 0.0;    wt[ 8]=0.05;
  point[ 9][0]=-NRM_GM1;point[ 9][1]= NRM_GMB;point[ 9][2]= 0.0;    wt[ 9]=0.05;
  point[10][0]= NRM_GM1;point[10][1]=-NRM_GMB;point[10][2]= 0.0;    wt[10]=0.05;
  point[11][0]=-NRM_GM1;point[11][1]=-NRM_GMB;point[11][2]= 0.0;    wt[11]=0.05;
  point[12][0]= NRM_GMB;point[12][1]= 0.0;    point[12][2]= NRM_GM1;wt[12]=0.05;
  point[13][0]= NRM_GMB;point[13][1]= 0.0;    point[13][2]=-NRM_GM1;wt[13]=0.05;
  point[14][0]=-NRM_GMB;point[14][1]= 0.0;    point[14][2]= NRM_GM1;wt[14]=0.05;
  point[15][0]=-NRM_GMB;point[15][1]= 0.0;    point[15][2]=-NRM_GM1;wt[15]=0.05;
  point[16][0]= 0.0;    point[16][1]= NRM_GM1;point[16][2]= NRM_GMB;wt[16]=0.05;
  point[17][0]= 0.0;    point[17][1]=-NRM_GM1;point[17][2]= NRM_GMB;wt[17]=0.05;
  point[18][0]= 0.0;    point[18][1]= NRM_GM1;point[18][2]=-NRM_GMB;wt[18]=0.05;
  point[19][0]= 0.0;    point[19][1]=-NRM_GM1;point[19][2]=-NRM_GMB;wt[19]=0.05;
}
/****************************************************************************/

void delta_one_D_twelve(double **point, double *wt)
{
/*
 * Twelve point gauss quadrature -- only for 1D
 */
  /********************** BEGIN EXECUTION ************************************/

 point[ 0][0]= 0.125233408511469; wt[ 0]=0.249147045813403/2.0;
 point[ 1][0]= 0.367831498998180; wt[ 1]=0.233492536538355/2.0;
 point[ 2][0]= 0.587317954286617; wt[ 2]=0.203167426723066/2.0;
 point[ 3][0]= 0.769902674194305; wt[ 3]=0.160078328543346/2.0;
 point[ 4][0]= 0.904117256370475; wt[ 4]=0.106939325995318/2.0;
 point[ 5][0]= 0.981560634246719; wt[ 5]=0.047175336386512/2.0;
 point[ 6][0]=-0.125233408511469; wt[ 6]=0.249147045813403/2.0;
 point[ 7][0]=-0.367831498998180; wt[ 7]=0.233492536538355/2.0;
 point[ 8][0]=-0.587317954286617; wt[ 8]=0.203167426723066/2.0;
 point[ 9][0]=-0.769902674194305; wt[ 9]=0.160078328543346/2.0;
 point[10][0]=-0.904117256370475; wt[10]=0.106939325995318/2.0;
 point[11][0]=-0.981560634246719; wt[11]=0.047175336386512/2.0;
}
/****************************************************************************/

void delta_midpoint(double **point, double *wt, int izone)
/*
 * Integration of a line (1D) or circle (2D) with equal spacing of
 * quadrature points. Weights are a function of r, and normalized 
 * to sum to 1. The weights are set so that the SURFACE of a sphere
 * is being integrated.
 *
 * 1D: integrate from -1 to 1 with weight  1
 * 2D: integrate from 0 to 1  and 0 to 2PI with weight  2*x/sqrt(1-x*x)
 */
{
  int i,j,k;
  double r, arc,  delr, delarc, wtsum;
  
  if (Ndim==1) {

     delr = 2.0 / ((double) MPsten_Npts_R[izone]);
     r = -1.0 + delr/2.0;
     wtsum = 0.0;
     for (i=0; i< MPsten_Npts_R[izone]; i++) {
        point[i][0] = r;
        wt[i] = 1.0;
        r += delr;
        wtsum += wt[i];
     }
     for (k=0; k<MPsten_Npts_R[izone]; k++) wt[k] /= wtsum;
  }
  else {

     delr = 1.0 / ((double) MPsten_Npts_R[izone]);
     r = delr/2.0;
 
     delarc = 2.0*PI/((double) MPsten_Npts_arc[izone]);

     wtsum = 0;

     for (i=0; i< MPsten_Npts_R[izone]; i++) {
        arc = delarc/2.0;
        for (j=0; j< MPsten_Npts_arc[izone]; j++) {
           k = j + i * MPsten_Npts_arc[izone];
           point[k][0] = r * cos(arc);
           point[k][1] = r * sin(arc);
           wt[k] = 2.0 * r / sqrt(1.0 - r*r);
           wtsum += wt[k];
           arc += delarc;
        }
        r += delr;
     }
     for (k=0; k<MPsten_Npts_R[izone]
                 *MPsten_Npts_arc[izone]; k++) wt[k] /= wtsum;
  }
}
/****************************************************************************/

void theta_midpoint(double **point, double *wt, int izone, int num_dim)
/*
 * Integration of a line (1D) or circle (2D) with equal spacing of
 * quadrature points. Weights are a function of r, and normalized 
 * to sum to 1. The weights are set so that the VOLUME of a sphere
 * is being integrated.
 *
 * 1D: integrate from -1 to 1 with weight  1-x*x
 * 2D: integrate from 0 to 1  and 0 to 2PI with weight  2*x*sqrt(1-x*x)
 */
{
  int i,j,k,l;
  double r, arc, phi, delr, delarc, delphi, wtsum;
  
  if (num_dim==1) {

     delr = 2.0 / ((double) MPsten_Npts_R[izone]);
     r = -1.0 + delr/2.0;
     wtsum = 0.0;
     for (i=0; i< MPsten_Npts_R[izone]; i++) {
        point[i][0] = r;
        wt[i] = 1.0 - r*r;
        r += delr;
        wtsum += wt[i];
     }
     for (k=0; k<MPsten_Npts_R[izone]; k++) wt[k] /= wtsum;
  }
  else if (num_dim==2) {

     delr = 1.0 / ((double) MPsten_Npts_R[izone]);
     r = delr/2.0;
 
     delarc = 2.0*PI/((double) MPsten_Npts_arc[izone]);

     wtsum = 0;

     for (i=0; i< MPsten_Npts_R[izone]; i++) {
        arc = delarc/2.0;
        for (j=0; j< MPsten_Npts_arc[izone]; j++) {
           k = j + i * MPsten_Npts_arc[izone];
           point[k][0] = r * cos(arc);
           point[k][1] = r * sin(arc);
           wt[k] = 2.0 * r * sqrt(1.0 - r*r);
           wtsum += wt[k];
           arc += delarc;
        }
        r += delr;
     }
     for (k=0; k<MPsten_Npts_R[izone]*MPsten_Npts_arc[izone]; k++) wt[k] /= wtsum;
  }
  else {

     delr = 1.0 / ((double) MPsten_Npts_R[izone]);
     delarc = 2.0*PI/((double) MPsten_Npts_arc[izone]);
     delphi =     PI/((double) MPsten_Npts_phi[izone]);

     wtsum = 0;
     r = delr/2.0;
 
     for (i=0; i< MPsten_Npts_R[izone]; i++) {
        arc = delarc/2.0;
        for (j=0; j< MPsten_Npts_arc[izone]; j++) {
           phi = delphi/2.0;
           for (l=0; l< MPsten_Npts_phi[izone]; l++) {
              k = l + (j + i * MPsten_Npts_arc[izone])*MPsten_Npts_phi[izone];
              point[k][0] = r * cos(arc) * sin(phi);
              point[k][1] = r * sin(arc) * sin(phi);
              point[k][2] = r * cos(phi);
              wt[k] = r*r*sin(phi);
              wtsum += wt[k];
              phi += delphi;
           }
           arc += delarc;
        }
        r += delr;
     }
     for (k=0; k< MPsten_Npts_R[izone]*MPsten_Npts_arc[izone]
                 *MPsten_Npts_phi[izone]; k++) wt[k] /= wtsum;
  }
}
/****************************************************************************/

void delta_seventy_two(double **point, double *wt)
{
  if (Iwrite_screen != SCREEN_NONE) printf("NO QUADRATURE SCHEME FOR SEVENTY_TWO PTS YET\n");
}

/****************************************************************************/

void get_radial_quadrature(double gauss_pt[], double gauss_wt[], int num_gp)
/*
 * This routine gives Gauss points and weights for integration of the
 * the radiall component of a sphere. I derived these so that the
 * Integral of {f(r) r^2 dr} is exact for f(r) a polynomial of
 * degree 2*num_gp.
 * Coupled with a quadrature schem for the surface of a sphere, this will
 * lead to a quadrature scheme for the volume of a sphere.
 */
{
  char *yo ="get_radial_quadrature";

  switch (num_gp) {

   case 1:
     gauss_wt[0]=1.0; gauss_pt[0]=  0.75;
     break;

   case 2:
     gauss_wt[0]=3.0235764623947340e-01; gauss_pt[0]= 4.5584815598877187e-01;
     gauss_wt[1]=6.9764235376052660e-01; gauss_pt[0]= 8.7748517734455789e-01;
     break;

   case 3:
     gauss_wt[0]=8.9852109025777560e-02; gauss_pt[0]= 2.9499779011155186e-01;
     gauss_wt[1]=4.3873880777961416e-01; gauss_pt[1]= 6.5299623396168860e-01;
     gauss_wt[2]=4.7140908319460828e-01; gauss_pt[2]= 9.2700597592686007e-01;
     break;

   case 4:
     gauss_wt[0]=3.1056722249766130e-02; gauss_pt[0]= 2.0414858210322620e-01;
     gauss_wt[1]=2.0590166151894631e-01; gauss_pt[1]= 4.8295270489575676e-01;
     gauss_wt[2]=4.3037636939768730e-01; gauss_pt[2]= 7.6139926244827938e-01;
     gauss_wt[3]=3.3266524683360026e-01; gauss_pt[3]= 9.5149945055304430e-01;
     break;

   case 5:
     gauss_wt[0]=1.2341475610014505e-02; gauss_pt[0]= 1.4894578705567910e-01;
     gauss_wt[1]=9.6166802174040925e-02; gauss_pt[1]= 3.6566652737615607e-01;
     gauss_wt[2]=2.6760048367080157e-01; gauss_pt[2]= 6.1011361294286182e-01;
     gauss_wt[3]=3.7859688569565880e-01; gauss_pt[3]= 8.2651967923358127e-01;
     gauss_wt[4]=2.4529435284948420e-01; gauss_pt[4]= 9.6542106008302164e-01;
     break;

   case 6:
     gauss_wt[0]=5.4932276161256360e-03; gauss_pt[0]= 1.1319438529429510e-01;
     gauss_wt[1]=4.7160892501962144e-02; gauss_pt[1]= 2.8431887515809884e-01;
     gauss_wt[2]=1.5386871464866680e-01; gauss_pt[2]= 4.9096358933247874e-01;
     gauss_wt[3]=2.8373156036950175e-01; gauss_pt[3]= 6.9756308376186738e-01;
     gauss_wt[4]=3.2212949788121320e-01; gauss_pt[4]= 8.6843605919817024e-01;
     gauss_wt[5]=1.8761610698253047e-01; gauss_pt[5]= 9.7409544508282242e-01;
     break;

   case 7:
     gauss_wt[0]=2.6780648140958850e-03; gauss_pt[0]= 8.8816841944067915e-02;
     gauss_wt[1]=2.4488780875018760e-02; gauss_pt[1]= 2.2648276869235773e-01;
     gauss_wt[2]=8.8266640653712615e-02; gauss_pt[2]= 3.9997850347348030e-01;
     gauss_wt[3]=1.8943914009760058e-01; gauss_pt[3]= 5.8599786894734329e-01;
     gauss_wt[4]=2.7520140687435119e-01; gauss_pt[4]= 7.5944588235711541e-01;
     gauss_wt[5]=2.7209646644871578e-01; gauss_pt[5]= 8.9691097453725066e-01;
     gauss_wt[6]=1.4782950023650519e-01; gauss_pt[6]= 9.7986722695055362e-01;
     break;

   default:
     if (Iwrite_screen != SCREEN_NONE) 
        printf("%s: Error! Number of radial Gauss Pts for Theta_fn must be 1-7...it is %d\n", yo, num_gp);
     if (num_gp>7) printf
          ("\t\tWhat, a 14th order polynomial isn't good enough for you???\n");
     exit(-1);
  }
}
/****************************************************************************/
