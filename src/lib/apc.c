/* FST - a Fast Shapelet Transformer
 *
   Copyright (C) 2006-2024 Sarod Yatawatta <sarod@users.sf.net>  
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 $Id$
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "shapelet.h"

/* various possibilities of the projection matrix */
#define PROJ_MAT_ZERO 2 /* when Ai full rank, I-A'(A A')^-1 A = 0 */
#define PROJ_MAT_I 1 /* when Ai=0 or ||Ai|| small, proj=I */
#define PROJ_MAT_NOR 0 /* other cases, where we do need to use the projection matrix */

/* lgrid,mgrid: l,m coordinate grid, each npix x 1 , l (deg),m (deg)
 * l0,m0 : coords of center (deg)
 * b: data vector, npix x 1 
 * P: projection matrix, modes x modes, output
 * x: solution vector, modes x 1, output
 * modes: number of modes
 * beta: shapelet scale 
 * n0: order, modes = n0*n0
 * Nt: number of threads
 *
 *
 * Calculate projection matrix and initial solution for this subproblem
 */
static int
calculate_projection_matrix_and_solution(double *lgrid, double *mgrid, double l0, double m0, int npix, float *b, float *P, float *x, int *sflag, int modes, double beta, int n0, int Nt) {

  if (modes != n0*n0) {
    fprintf(stderr,"%s: %d: number of modes should agree with model order\n",__FILE__,__LINE__);
    exit(1);
  }

  double *l,*m,*Av;
  /* allocate memory for l,m grid points (radians) and basis functions */
  if ((l=(double*)calloc((size_t)npix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((m=(double*)calloc((size_t)npix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  /* copy l,m coords */
  for (int ci=0; ci<npix; ci++) {
    l[ci]=-(lgrid[ci]-l0)*M_PI/180.0;
    m[ci]=(mgrid[ci]-m0)*M_PI/180.0;
  }

  /* Av (npix x modes) storage will be allocated within the routine */
#ifndef HAVE_CUDA
  calculate_mode_vectors_thread(l,m,npix,beta,n0,&Av,Nt);
#else
  calculate_mode_vectors_cuda(l,m,npix,beta,n0,&Av);
#endif /* HAVE_CUDA */

  *sflag=PROJ_MAT_NOR;
  /* check norm of Av, if too small skip calculation and set
   * projection to I and original solution to zero */
  double a_norm=my_dnrm2(modes*npix,Av);
  if (a_norm<1e-9) {
    *sflag=PROJ_MAT_I;
    for (int ci=0; ci<modes; ci++) {
      x[ci]=0.0f;
    }
    free(l);
    free(m);
    free(Av);

    printf("skipping finding projection/initial solution because norm is too low\n");
    return 0;
  }
  /* if Av is full rank, projection matrix = 0 */
  if (npix > modes) {
    *sflag=PROJ_MAT_ZERO;
  }

  if (*sflag==PROJ_MAT_NOR) {
  /* projection P = eye(modes) - A' * (A * A')^{-1} A */
  double *AAt;
  if ((AAt=(double*)calloc((size_t)npix*npix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }

  /* find product AAt = A * A' */
  my_dgemm('N','T',npix,npix,modes,1.0,Av,npix,Av,npix,0.0,AAt,npix);

  double w[1],*WORK,*U,*S,*VT;
  int status,lwork=1;
  /* find inv(A * A') using SVD */
  /* allocate memory for SVD */
  if ((U=(double*)calloc((size_t)npix*npix,sizeof(double)))==0) {
    printf("%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((VT=(double*)calloc((size_t)npix*npix,sizeof(double)))==0) {
    printf("%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((S=(double*)calloc((size_t)npix,sizeof(double)))==0) {
    printf("%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }


  /* memory for SVD: use first location of Bi */
  status=my_dgesvd('A','A',npix,npix,AAt,npix,S,U,npix,VT,npix,w,-1);
  if (!status) {
    lwork=(int)w[0];
  } else {
    printf("%s: %d: LAPACK error %d\n",__FILE__,__LINE__,status);
    exit(1);
  }
  if ((WORK=(double*)calloc((size_t)lwork,sizeof(double)))==0) {
    printf("%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  status=my_dgesvd('A','A',npix,npix,AAt,npix,S,U,npix,VT,npix,WORK,lwork);
  if (status) {
    printf("%s: %d: LAPACK error %d\n",__FILE__,__LINE__,status);
    exit(1);
  }
  /* find 1/singular values, and multiply columns of U with new singular values */
  for (int ci=0; ci<npix; ci++) {
    if (S[ci]>1e-9) {
     S[ci]=1.0/S[ci];
     my_dscal(npix,S[ci],&U[ci*npix]);//only do this for nonzero S[]
    } else {
     S[ci]=0.0;
     memset(&U[ci*npix],0,npix*sizeof(double));
    }
  }

  /* find product U 1/S V^T */
  my_dgemm('N','N',npix,npix,npix,1.0,U,npix,VT,npix,0.0,AAt,npix);
  free(WORK);
  free(S);
  free(U);
  free(VT);

  /* now AAt=(A A')^{-1}, size npix x npix */
  /* find A' * (A A')^{-1}, size modes x npix */
  double *A_AAt_inv;
  if ((A_AAt_inv=(double*)calloc((size_t)modes*npix,sizeof(double)))==0) {
    printf("%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  my_dgemm('T','N',modes,npix,npix,1.0,Av,npix,AAt,npix,0.0,A_AAt_inv,modes);
  /* find (A' *(A A')^{-1}) * A, size modes x modes */
  double *A_AAt_inv_A;
  if ((A_AAt_inv_A=(double*)calloc((size_t)modes*modes,sizeof(double)))==0) {
    printf("%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  my_dgemm('N','N',modes,modes,npix,1.0,A_AAt_inv,modes,Av,npix,0.0,A_AAt_inv_A,modes);


  free(AAt);
  free(A_AAt_inv);

  /* Projection = I - A' * (A A')^{-1} * A */
  my_dscal(modes*modes,-1.0,A_AAt_inv_A);
  for (int ci=0; ci<modes; ci++) {
    A_AAt_inv_A[ci*modes+ci]+=1.0;
  }
  for (int ci=0; ci<modes*modes; ci++) {
    P[ci]=(float)A_AAt_inv_A[ci];
  }
  free(A_AAt_inv_A);

#ifdef DEBUG
  printf("P=[\n");
  for (int ci=0; ci<modes; ci++) {
  for (int cj=0; cj<modes; cj++) {
    printf("%f ",P[ci*modes+cj]);
  }
  printf("\n");
  }
  printf("];\n");
#endif
  } /* PROJ_MAT_NOR */

  /* x = pinv(A)*b = pinv(A' * A) * A' * b */
  double *bd,*xd;
  if ((bd=(double*)calloc((size_t)npix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((xd=(double*)calloc((size_t)modes,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }

  for (int ci=0; ci<npix; ci++) {
    bd[ci]=(double)b[ci];
  }

  /* find least squares estimate for x */
  //lsq_lapack(Av,bd,xd,npix,modes);
  // normalize b ~ ||b||=1 and scale back solution x
  elasticnet_fista(Av,bd,xd,npix,modes,1e-3,1e-9,300);

  printf("||A||=%lf ||x||=%lf\n",a_norm,my_dnrm2(modes,xd));
  for (int ci=0; ci<modes; ci++) {
    x[ci]=(float)xd[ci];
  }

  free(bd);
  free(xd);

  free(l);
  free(m);
  free(Av);
  return 0;
}



/* calculate model (coeffients given by z) 
 * over a subimage b, the right offset should be given to the subimage
 * imgrid: coordinate grid, 4xnpix x 1 , each 4 values: l (deg),m (deg),Freq,Stokes (last two not used)
 * l0,m0: coords of image center (deg)
 * b: npix x 1 pixel values (output)
 * z: modes x 1, model
 * beta: shapelet scale
 * n0: modes = n0*n0
 * Nt: number of threads
 */
static int
evaluate_model_over_subimage(double *imgrid, double l0, double m0, int npix, double *b, double *z, int modes, double beta, int n0, int Nt) 
{

  if (modes != n0*n0) {
    fprintf(stderr,"%s: %d: number of modes should agree with model order\n",__FILE__,__LINE__);
    exit(1);
  }

  double *l,*m,*Av;
  /* allocate memory for l,m grid points (radians) and basis functions */
  if ((l=(double*)calloc((size_t)npix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((m=(double*)calloc((size_t)npix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  /* copy l,m coords */
  for (int ci=0; ci<npix; ci++) {
    l[ci]=-(imgrid[4*ci]-l0)*M_PI/180.0;
    m[ci]=(imgrid[4*ci+1]-m0)*M_PI/180.0;
  }

  /* Av storage will be allocated within the routine */
#ifndef HAVE_CUDA
  calculate_mode_vectors_thread(l,m,npix,beta,n0,&Av,Nt);
#else
  calculate_mode_vectors_cuda(l,m,npix,beta,n0,&Av);
#endif /* HAVE_CUDA */

  for (int ci=0; ci<modes; ci++) {
    my_daxpy(npix,&Av[ci*npix],z[ci],b);
  }
  
  /*
  for (int ci=0; ci<npix; ci++) {
    b[ci]=Av[2*npix+ci];
  }
  */
  

  free(l);
  free(m);
  free(Av);
  return 0;
}

/* divide d into J equally sized amounts, low and high values are stored 
 * in arrays lowp, highp
 * lowp,highp: Jx1 arrays 
 */
static long int
divide_into_subsets(int J,long int d, long int *lowp, long int *highp) {
  if (J>d) {
    fprintf(stderr,"%s: %d: number of subtasks cannot be higher than the total\n",__FILE__,__LINE__);
    exit(1);
  }

  long int p =d / J;
  long int r =d % J;
  long int counter=1;
  for (int ci=0; ci<r; ci++) {
    lowp[ci]=counter;
    highp[ci]=counter+p;
    counter+=p+1;
  }
  for (int ci=r; ci<J; ci++) {
    lowp[ci]=counter;
    highp[ci]=counter+p-1;
    counter+=p;
  }

  return p;
}

/* solve large linear system using accelerated projection based consensus
 * (APC)
 */
int
apc_decompose_fits_file(char* filename, double cutoff, int *Nx, int *Ny, double *beta, int *M, int *n0, double **img, double **av, double **z, position *cen, char* outfile, int J, int Nt) {

  /* open the file once, get all the metadata to make a plan to divide the pixels */
  io_buff fitsref;
  int status=0;
  long int increment[4]={1,1,1,1};
  int null_flag=0;
  float nullval=0.0;

  /* WCSlib */
  char *header;
  int ncard,nreject,nwcs;
  int stat[NWCSFIX];

  /* shapelet mode parameters (must be pre-defined) */
  if (*M==-1) {
    if (*n0 > 0) {
      /* depending on M, update n0 ~ sqrt(M), and M =n0*n0 */
      *M=(*n0)*(*n0);
    } else {
      fprintf(stderr,"Warning: overriding default value for number of modes %d\n",*M);
      *M=100;
      *n0=10;
    }
  } else {
    /* update n0 ~=sqrt(M) and M=n0*n0 */
    *n0=(int)(sqrt((double)*M)+1);
    *M=(*n0)*(*n0);
  }
  if (*beta==-1.0) {
    fprintf(stderr,"Warning: overriding default value for scale (beta) %lf\n",*beta);
    *beta=0.01;
  }
  printf("M=%d n0=%d beta=%lf\n",*M,*n0,*beta);
  int modes=*M;
  

  /* image data, each Npix x 1 */
  float **b;
  if ((b=(float**)calloc((size_t)J,sizeof(float*)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  /* projection matrices, each modes x modes */
  float **P;
  if ((P=(float**)calloc((size_t)J,sizeof(float*)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  /* local solution (per subtask), each modes x 1 */
  float **xb;
  if ((xb=(float**)calloc((size_t)J,sizeof(float*)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  /* flag to indicate the usability of each subtask, default is PROJ_MAT_NOR */
  int *sflag;
  if ((sflag=(int*)calloc((size_t)J,sizeof(int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }

  /* for division into (almost equal) subimages  - reconstruction,
   * keep this lower than the number of columns */
  int Jimg=J;
  long int *lowp,*highp;
  if ((lowp=(long int*)calloc((size_t)Jimg,sizeof(long int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((highp=(long int*)calloc((size_t)Jimg,sizeof(long int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }


  fits_open_file(&fitsref.fptr,filename,READONLY,&status);
  if (status) {
     fits_report_error(stderr,status);
     exit(1);
  }
  /* get WCS */
  if ((status = fits_hdr2str(fitsref.fptr, 1, NULL, 0, &header, &ncard, &status))) {
     fits_report_error(stderr, status);
     exit(1);
  }
  if ((status = wcspih(header, ncard, WCSHDR_all, 2, &nreject, &nwcs, &fitsref.wcs))) {
      fprintf(stderr, "wcspih ERROR %d, ignoring WCS\n", status);
      exit(1);
  }
  free(header);
  status=0;
  /* Fix non-standard WCS keyvalues. */
  if ((status = wcsfix(7, 0, fitsref.wcs, stat))) {
     fprintf(stderr,"wcsfix ERROR, status returns: (");
     for (int ii = 0; ii < NWCSFIX; ii++) {
          fprintf(stderr,ii ? ", %d" : "%d", stat[ii]);
     }
     printf(")\n\n");
  }
  if ((status = wcsset(fitsref.wcs))) {
     fprintf(stderr, "wcsset ERROR %d:\n", status);
     exit(1);
  }
  /* set scale */
  double bscale,bzero;
  fits_read_key(fitsref.fptr,TDOUBLE,"BSCALE",&bscale,0,&status);
  if (status) { status=0; bscale=1.0; }
  fits_read_key(fitsref.fptr,TDOUBLE,"BZERO",&bzero,0,&status);
  if (status) { status=0; bzero=0.0; }
  fits_set_bscale(fitsref.fptr, bscale, bzero, &status);
  /* get dimensions */
  int naxis,bitpix;
  fits_get_img_dim(fitsref.fptr, &naxis, &status);
  /* fix zero length axes */
  if (naxis<4) naxis=4;
  /* allocate mem for axes info */
  if ((fitsref.arr_dims.d=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((fitsref.arr_dims.lpix=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((fitsref.arr_dims.hpix=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  /* get axis sizes */
  fits_get_img_size(fitsref.fptr, naxis, fitsref.arr_dims.d, &status);
  if (status) fits_report_error(stderr,status);
  /* reset zero length axes to 1 */
  if (fitsref.arr_dims.d[2]==0) { fitsref.arr_dims.d[2]=1; }
  if (fitsref.arr_dims.d[3]==0) { fitsref.arr_dims.d[3]=1; }
  fitsref.arr_dims.naxis=naxis;
  fitsref.arr_dims.tol=cutoff;
  /* get data type */
  fits_get_img_type(fitsref.fptr, &bitpix, &status);
  if (status) fits_report_error(stderr,status);
  /* use the whole x,y-axes */
  /* correct the coordinates for 1 indexing */
  fitsref.arr_dims.lpix[0]=1;
  fitsref.arr_dims.hpix[0]=fitsref.arr_dims.d[0];
  fitsref.arr_dims.lpix[1]=1;
  fitsref.arr_dims.hpix[1]=fitsref.arr_dims.d[1];
  /* only work with stokes I */
  fitsref.arr_dims.hpix[2]=fitsref.arr_dims.hpix[3]=fitsref.arr_dims.lpix[2]=fitsref.arr_dims.lpix[3]=1;
  /* since data is row major, divide the data into rows (axis 1 or y axis)
  * to distribute the work */
  if (Jimg>fitsref.arr_dims.d[1]) { Jimg=fitsref.arr_dims.d[1]; }
  divide_into_subsets(Jimg,fitsref.arr_dims.d[1],lowp,highp);

  /* find l,m of image center because 
   * we need to shit the l,m grid to have this as origin (0,0) */
  double cpixelc[4], cimgc[4], cworldc[4], cphic[1], cthetac[1];
  int cstatc[1];
  cpixelc[0]=(0.5*(double)fitsref.arr_dims.d[0]);
  cpixelc[1]=(0.5*(double)fitsref.arr_dims.d[1]);
  cpixelc[2]=cpixelc[3]=1.0;
  if ((status = wcsp2s(fitsref.wcs, 1, naxis, cpixelc, cimgc, cphic, cthetac,
       cworldc, cstatc))) {
       fprintf(stderr,"wcsp2s ERROR %2d\n", status);
       /* Handle Invalid pixel coordinates. */
       if (status == 8) status = 0;
  }
  double l0,m0;
  l0=cimgc[0];
  m0=cimgc[1];
  *Nx=fitsref.arr_dims.d[0];
  *Ny=fitsref.arr_dims.d[1];

  /*****************************************************************/
  /* read the whole image and create coordinate for the whole image
   * here */
  long int totalpix=(fitsref.arr_dims.hpix[0]-fitsref.arr_dims.lpix[0]+1)
      *(fitsref.arr_dims.hpix[1]-fitsref.arr_dims.lpix[1]+1);
  float *image;
  if ((image=(float*)calloc((size_t)totalpix,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  /* read whole-image */
  fits_read_subset(fitsref.fptr, TFLOAT, fitsref.arr_dims.lpix, fitsref.arr_dims.hpix, increment,
        &nullval, image, &null_flag, &status);
  if (status) { fits_report_error(stderr,status); exit(1); }
  /* create grid for whole-image */
  double *pixelc,*imgc,*worldc,*phic,*thetac;
  int *statc;
  if ((pixelc=(double*)calloc((size_t)totalpix*4,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((imgc=(double*)calloc((size_t)totalpix*4,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((worldc=(double*)calloc((size_t)totalpix*4,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((phic=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((thetac=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((statc=(int*)calloc((size_t)totalpix,sizeof(int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }

  int kk=0;
  for (int jj=fitsref.arr_dims.lpix[1];jj<=fitsref.arr_dims.hpix[1];jj++) {
  for (int ii=fitsref.arr_dims.lpix[0];ii<=fitsref.arr_dims.hpix[0];ii++) {
             pixelc[kk+0]=(double)ii;
             pixelc[kk+1]=(double)jj;
             pixelc[kk+2]=pixelc[kk+3]=1.0;
             kk+=4;
   }
  }
  if ((status = wcsp2s(fitsref.wcs, totalpix, naxis, pixelc, imgc, phic, thetac,
       worldc, statc))) {
       fprintf(stderr,"wcsp2s ERROR %2d\n", status);
       /* Handle Invalid pixel coordinates. */
       if (status == 8) status = 0;
  }

  /* copy l,m coords (size 4 pixelsx1) into smaller arrays */
  double *lgrid,*mgrid;
  if ((lgrid=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  } 
  if ((mgrid=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  } 
  for(int ci=0; ci<totalpix; ci++) {
     lgrid[ci]=imgc[4*ci]; // l (deg)
     mgrid[ci]=imgc[4*ci+1]; // m (deg)
  }
  free(pixelc);
  free(imgc);
  free(worldc);
  free(phic);
  free(thetac);
  free(statc);
  /*****************************************************************/
  printf("original pixels %ld\n",totalpix);
  totalpix=(fitsref.arr_dims.hpix[0]-fitsref.arr_dims.lpix[0]+1)
      *(fitsref.arr_dims.hpix[1]-fitsref.arr_dims.lpix[1]+1)/J; // tail will be truncated
  printf("subtask pixels %ld\n",totalpix);

  for (int ci=0; ci<J; ci++) {
    /*****************************************************************/
    /* select the rows by incremental access of increment J */
    if ((b[ci]=(float*)calloc((size_t)totalpix,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    double *lcoord,*mcoord;
    if ((lcoord=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    if ((mcoord=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }

    /* copy subimage from full image */
    my_scopy(totalpix,&image[ci],J,b[ci],1);
    my_dcopy(totalpix,&lgrid[ci],J,lcoord,1);
    my_dcopy(totalpix,&mgrid[ci],J,mcoord,1);
    /* check the tail of b[] to see if any are zero (cannot be used for fitting) */
    int tail=0;
    for (int nt=0; nt<J; nt++) {
      if (fabsf(b[ci][totalpix-nt-1])<=1e-5f){
        tail++;
      }
    }

    if ((P[ci]=(float*)calloc((size_t)modes*modes,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    if ((xb[ci]=(float*)calloc((size_t)modes,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    /* use the coordinate values to calculate the basis functions, and
     * the projection matrix, and the initial solution */
    calculate_projection_matrix_and_solution(lcoord,mcoord,l0,m0,totalpix-tail,b[ci],P[ci],xb[ci], &sflag[ci], modes,*beta, *n0, Nt);

    printf("norm %d %f\n",ci,my_snrm2(modes*modes,P[ci]));
    /* if projection matrix is not used, free up memory */
    if (sflag[ci] != PROJ_MAT_NOR) {
      free(P[ci]);
    }

    free(lcoord);
    free(mcoord);
    /*****************************************************************/
  }
  free(image);
  free(lgrid);
  free(mgrid);

  /*****************************************************************/
  float gamma=0.1f;
  float eta=0.1f;
  int Nadmm=100;
  /* solution */
  float *x;
  if ((x=(float*)calloc((size_t)modes,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  float *xdiff;
  if ((xdiff=(float*)calloc((size_t)modes,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  float *xold;
  if ((xold=(float*)calloc((size_t)modes,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }

  /* count subtasks where valid projection matrix (!PROJ_MAT_I) */
  int Jvalid=0;
  for (int ci=0; ci<J; ci++) {
    if (sflag[ci]!=PROJ_MAT_I) {
      Jvalid++;
    }
  }
  printf("valid subtasks %d\n",Jvalid);
  /* ADMM iterations */
  for (int admm=0; admm<Nadmm; admm++) {
    /* workers update their estimate */
    for (int ci=0; ci<J; ci++) {
       /* x_i <= x_i + gamma Proj_i (x - x_i) */
       if (sflag[ci]==PROJ_MAT_NOR) {
        my_scopy(modes,x,1,xdiff,1);
        my_saxpy(modes,xb[ci],-1.0f,xdiff);
        my_sgemv('N',modes,modes,gamma,P[ci],modes,xdiff,1,1.0f,xb[ci],1);
       } else if (sflag[ci]==PROJ_MAT_I) {
        my_scopy(modes,x,1,xdiff,1);
        my_saxpy(modes,xb[ci],-1.0f,xdiff);
        my_saxpy(modes,xdiff,gamma,xb[ci]);
       } //(sflag[ci]==PROJ_MAT_ZERO) no change
    }
    /* find mean x_i */
    memset(xdiff,0,modes*sizeof(float));
    for (int ci=0; ci<J; ci++) {
      if (sflag[ci]!=PROJ_MAT_I) {
       my_saxpy(modes,xb[ci],1.0f,xdiff);
      }
    }
    if(Jvalid>0) {
     my_sscal(modes,1.0f/(float)J,xdiff);
    }
    /* update current estimate (with momentum) */
    my_scopy(modes,x,1,xold,1); /* backup old for bookkeeping */
    /* x <= eta xnew + (1-eta) x */
    my_sscal(modes,1.0f-eta,x);
    my_saxpy(modes,xdiff,eta,x);

    /* find ||xold-x|| */
    my_saxpy(modes,x,-1.0f,xold);
    printf("%d %f\n",admm,my_snrm2(modes,xold));
  }
  free(xdiff);
  free(xold);
  /*****************************************************************/

  /* copy the solution */
  if ((*av=(double*)calloc((size_t)modes,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  for (int ci=0; ci<modes; ci++) {
    (*av)[ci]=(double)x[ci];
  }
  for (int ci=0; ci<J; ci++) {
    free(b[ci]);
    if (sflag[ci]==PROJ_MAT_NOR) {
     free(P[ci]);
    }
    free(xb[ci]);
  }
  free(b);
  free(P);
  free(xb);
  free(x);
  free(sflag);

  /* recreate pixel values based on the model, also read image */
  if ((*z=(double*)calloc((size_t) fitsref.arr_dims.d[0]*fitsref.arr_dims.d[1],sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  if ((*img=(double*)calloc((size_t) fitsref.arr_dims.d[0]*fitsref.arr_dims.d[1],sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
  }
  long int offset=0;
  for (int ci=0; ci<Jimg; ci++) {
    /*****************************************************************/
    /* select the rows (1 indexing) */
    fitsref.arr_dims.lpix[1]=lowp[ci];
    fitsref.arr_dims.hpix[1]=highp[ci];
    long int totalpix=(fitsref.arr_dims.hpix[0]-fitsref.arr_dims.lpix[0]+1)
      *(fitsref.arr_dims.hpix[1]-fitsref.arr_dims.lpix[1]+1);
    /* create grid for this sub-image */
    double *pixelc,*imgc,*worldc,*phic,*thetac;
    int *statc;
    if ((pixelc=(double*)calloc((size_t)totalpix*4,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    if ((imgc=(double*)calloc((size_t)totalpix*4,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    if ((worldc=(double*)calloc((size_t)totalpix*4,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    if ((phic=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    if ((thetac=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }
    if ((statc=(int*)calloc((size_t)totalpix,sizeof(int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }

    int kk=0;
    for (int jj=fitsref.arr_dims.lpix[1];jj<=fitsref.arr_dims.hpix[1];jj++) { 
    for (int ii=fitsref.arr_dims.lpix[0];ii<=fitsref.arr_dims.hpix[0];ii++) {
             pixelc[kk+0]=(double)ii;
             pixelc[kk+1]=(double)jj;
             pixelc[kk+2]=pixelc[kk+3]=1.0;
             kk+=4;
    }
    }
    if ((status = wcsp2s(fitsref.wcs, totalpix, naxis, pixelc, imgc, phic, thetac,
       worldc, statc))) {
       fprintf(stderr,"wcsp2s ERROR %2d\n", status);
       /* Handle Invalid pixel coordinates. */
       if (status == 8) status = 0;
    }

    fits_read_subset(fitsref.fptr, TDOUBLE, fitsref.arr_dims.lpix, fitsref.arr_dims.hpix, increment,
        &nullval, &((*img)[offset]), &null_flag, &status);


    /* use the coordinate values to calculate the basis functions, and
     * the model image */
    evaluate_model_over_subimage(imgc, l0, m0, totalpix, &((*z)[offset]), *av, modes, *beta, *n0, Nt);
    offset+=totalpix;
    free(pixelc);
    free(imgc);
    free(worldc);
    free(phic);
    free(thetac);
    free(statc);
    /*****************************************************************/
  }

  free(lowp);
  free(highp);

  /* update image center information */
  cworldc[0]*=24.0/360.0;
  cen->ra_h=(int)((int)cworldc[0]%24);
  cen->ra_m=(int)((cworldc[0]-cen->ra_h)*60.0);
  cen->ra_s=(cworldc[0]-cen->ra_h-cen->ra_m/60.0)*3600.0;
  cen->dec_d=(int)((int)cworldc[1]%180);
  cen->dec_m=(int)((cworldc[1]-cen->dec_d)*60.0);
  cen->dec_s=(cworldc[1]-cen->dec_d-cen->dec_m/60.0)*3600.0;

  /* find the scale difference between input image and model image */
  double img_norm=my_dnrm2(fitsref.arr_dims.d[0]*fitsref.arr_dims.d[1],*img);
  double model_norm=my_dnrm2(fitsref.arr_dims.d[0]*fitsref.arr_dims.d[1],*z);
  if (model_norm>0.0) {
   /* rescale model to match image */
   my_dscal(fitsref.arr_dims.d[0]*fitsref.arr_dims.d[1],img_norm/model_norm,*z);
   /* also rescale model coefficients */
   my_dscal(modes,img_norm/model_norm,*av);
  }


  /* if output file is given, write model to output */
  if (outfile) {
    fitsref.arr_dims.lpix[0]=1;
    fitsref.arr_dims.hpix[0]=fitsref.arr_dims.d[0];
    fitsref.arr_dims.lpix[1]=1;
    fitsref.arr_dims.hpix[1]=fitsref.arr_dims.d[1];
    write_fits_file(filename, outfile,*z,fitsref);
  }

  fits_close_file(fitsref.fptr,&status);
  if (status) fits_report_error(stderr,status);
  wcsfree(fitsref.wcs);
  free(fitsref.wcs);
  free(fitsref.arr_dims.d);
  free(fitsref.arr_dims.lpix);
  free(fitsref.arr_dims.hpix);

  return 0;
}
