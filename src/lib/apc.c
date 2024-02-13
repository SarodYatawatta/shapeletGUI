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


/* imgrid: coordinate grid, 4xnpix x 1 , each 4 values: l (deg),m (deg),Freq,Stokes (last two not used)
 * b: data vector, npix x 1 
 * P: projection matrix, modes x modes, output
 * x: solution vector, modes x 1, output
 * modes: number of modes
 * beta: shapelet scale 
 * n0: order, modes = n0*n0
 *
 *
 * Calculate projection matrix and initial solution for this subproblem
 */
static int
calculate_projection_matrix_and_solution(double *imgrid, int npix, float *b, float *P, float *x, int modes, double beta, int n0) 
{

  if (modes != n0*n0) {
    fprintf(stderr,"%s: %d: number of modes should agree with model order\n",__FILE__,__LINE__);
    exit(1);
  }

  double *l,*m,*Av;
  /* allocate memory for l,m grid points (radians) and basis functions */
  if ((l=(double*)calloc((size_t)npix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      return 1;
  }
  if ((m=(double*)calloc((size_t)npix,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      return 1;
  }
  /* copy l,m coords */
  for (int ci=0; ci<npix; ci++) {
    l[ci]=-imgrid[4*ci]*M_PI/180.0;
    m[ci]=imgrid[4*ci+1]*M_PI/180.0;
  }

  /* Av storage will be allocated within the routine */
  calculate_mode_vectors_bi(l,m,npix,beta,n0,&Av);

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
  lsq_lapack(Av,bd,xd,npix,modes);

  for (int ci=0; ci<modes; ci++) {
    x[ci]=(float)xd[ci];
    printf("ci=%d x[ci]=%f\n",ci,x[ci]);
  }

  free(bd);
  free(xd);

  free(l);
  free(m);
  free(Av);
  return 0;
}



/* solve large linear system using accelerated projection based consensus
 * (APC)
 */
int
apc_decompose_fits_file(char* filename, double cutoff, double *beta, int *M, int *n0, double **av, double **z) {

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
  

  /* number of subtasks */
  int J=10;
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
  /* use the whole y-axis */
  /* correct the coordinates for 1 indexing */
  fitsref.arr_dims.lpix[1]=1;
  fitsref.arr_dims.hpix[1]=fitsref.arr_dims.d[1];
  /* only work with stokes I */
  fitsref.arr_dims.hpix[2]=fitsref.arr_dims.hpix[3]=fitsref.arr_dims.lpix[2]=fitsref.arr_dims.lpix[3]=1;
  /* since data is column major, divide the data into columns (axis 0 or x axis)
   * to distribute the work */
  long int Ncol=(fitsref.arr_dims.d[0]+J-1)/J;
  printf("divide columns %ld into %ld\n",fitsref.arr_dims.d[0],Ncol);

  for (int ci=0; ci<J; ci++) {
    /*****************************************************************/
    /* select the rows (1 indexing) */
    fitsref.arr_dims.lpix[0]=ci*Ncol+1;
    fitsref.arr_dims.hpix[0]=((ci+1)*Ncol<fitsref.arr_dims.d[0]?(ci+1)*Ncol:fitsref.arr_dims.d[0]);
    long int totalpix=(fitsref.arr_dims.hpix[0]-fitsref.arr_dims.lpix[0]+1)
      *(fitsref.arr_dims.hpix[1]-fitsref.arr_dims.lpix[1]+1);
    printf("%ld %ld %ld\n",fitsref.arr_dims.lpix[0],fitsref.arr_dims.hpix[0],totalpix);
    if ((b[ci]=(float*)calloc((size_t)totalpix,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      return 1;
    }
    /* read sub-image */
    fits_read_subset(fitsref.fptr, TFLOAT, fitsref.arr_dims.lpix, fitsref.arr_dims.hpix, increment,
        &nullval, b[ci], &null_flag, &status);
    if (status) fits_report_error(stderr,status);
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
    for (int ii=fitsref.arr_dims.lpix[0];ii<=fitsref.arr_dims.hpix[0];ii++)
     for (int jj=fitsref.arr_dims.lpix[1];jj<=fitsref.arr_dims.hpix[1];jj++) {
             pixelc[kk+0]=(double)ii;
             pixelc[kk+1]=(double)jj;
             pixelc[kk+2]=(double)1.0;
             pixelc[kk+3]=(double)1.0;
             kk+=4;
    }
    if ((status = wcsp2s(fitsref.wcs, totalpix, naxis, pixelc, imgc, phic, thetac,
       worldc, statc))) {
       fprintf(stderr,"wcsp2s ERROR %2d\n", status);
       /* Handle Invalid pixel coordinates. */
       if (status == 8) status = 0;
    }
    if ((P[ci]=(float*)calloc((size_t)modes*modes,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      return 1;
    }
    if ((xb[ci]=(float*)calloc((size_t)modes,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      return 1;
    }
    /* use the coordinate values to calculate the basis functions, and
     * the projection matrix, and the initial solution */
    calculate_projection_matrix_and_solution(imgc,totalpix,b[ci],P[ci],xb[ci],modes,*beta, *n0);


    free(pixelc);
    free(imgc);
    free(worldc);
    free(phic);
    free(thetac);
    free(statc);
    /*****************************************************************/
  }


  fits_close_file(fitsref.fptr,&status);
  if (status) fits_report_error(stderr,status);
  wcsfree(fitsref.wcs);
  free(fitsref.wcs);
  free(fitsref.arr_dims.d);
  free(fitsref.arr_dims.lpix);
  free(fitsref.arr_dims.hpix);

  for (int ci=0; ci<J; ci++) {
    free(b[ci]);
    free(P[ci]);
    free(xb[ci]);
  }
  free(b);
  free(P);
  free(xb);
  exit(1);
  return 0;
}
