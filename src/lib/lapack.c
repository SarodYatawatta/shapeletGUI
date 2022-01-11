/* FST - a Fast Shapelet Transformer
 *
   Copyright (C) 2006-2011 Sarod Yatawatta <sarod@users.sf.net>  
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

#include "shapelet.h"

/* solve linear least squares problem using divide and conquor SVD
 */
int 
dgelsd(int M, int N, int NRHS, double *A, int LDA, double *B, int LDB, double *S,
	double RCOND, int *RANK, double *WORK, int LWORK, int *IWORK) {
  int info;
	extern void dgelsd_(int *M,int *N,int *NRHS, double *A, int *LDA, double *B,
		int *LDB, double *S, double *RCOND, int *RANK, double *WORK, int *LWORK, int *IWORK, int *INFO);
	dgelsd_(&M,&N,&NRHS, A, &LDA, B,
		&LDB, S, &RCOND, RANK, WORK, &LWORK, IWORK, &info);
	return info;
}

/* get machine precision for noise eigenvalues */
static double
dlamch (char CMACH)
{
				  extern  double  dlamch_ (char *CMACHp);
					  return  dlamch_ (&CMACH);
}


/* choose problem parameters */
int
ilaenv(int ISPEC, char *NAME, char *OPTS, int N1, int N2, int N3, int N4) {
	/* for dgelsdm ISPEC=9 */
	extern int ilaenv_(int *ISPEC, char *NAME, char *OPTS, int *N1, int *N2, int *N3, int *N4);
	return ilaenv_(&ISPEC,NAME,OPTS,&N1,&N2,&N3,&N4);
}

/* y = a.x + y */
void
daxpy(int N, double *x, double a, double *y) {
		int i=1; /* strides */
		extern void daxpy_(int *N, double *alpha, double *x, int *incx, double *y, int *incy);
		daxpy_(&N,&a,x,&i,y,&i);
}

/* y = x */
void
dcopy(int N, double *x, double *y) {
  extern void dcopy_(int *N, double *x, int *incx, double *y, int *incy);
  int i=1;
  dcopy_(&N,x,&i,y,&i);
}

/* scale */
void
dscal(int N, double a, double *x) {
  extern void dscal_(int *N, double *alpha, double *x, int *incx);
  int i=1;
  dscal_(&N,&a,x,&i);
}

/* norm */
double
dnrm2(int N, double *x) {
  extern double  dnrm2_(int *N, double *x, int *incx);
  int i=1;
  return(dnrm2_(&N,x,&i));
}

/* max id, start from 1... */
int
idamax(int N, double *x, int incx) {
    extern int idamax_(int *N, double *x, int *incx);
    return idamax_(&N,x,&incx);
}

/* min id, start from 1... */
int
idamin(int N, double *x, int incx) {
    extern int idamin_(int *N, double *x, int *incx);
    return idamin_(&N,x,&incx);
}

/* solve the linear least squares problem using LAPACK */
/* min_x |Ax-b|_2 norm */
/* A: N by M, N> M 
 * b: N by 1 vector 
 * x: M by 1 vector */
int 
lsq_lapack(double *Av,double *b,double *x, int N, int M) {
	int smlsiz,nlvl,vec_len;
	double  *sigma, *WORK, *A, *bb;
	int *IWORK,LWORK,info,rank;
	double w[1];
#ifdef DEBUG
	int i;
#endif

	/* create a copy of the matrix as it is destroyed */
  if ((A=(double*)calloc((size_t)(M*N),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	memcpy((void*)A,(void*)Av,(size_t)(M*N)*sizeof(double));


	smlsiz=ilaenv(9,"DGELSD","",N,M,1,-1);
	printf("smlsiz=%d\n",smlsiz);
	nlvl=(int)(log((double)(M)/(smlsiz+1.0))/log(2.0))+1;
	if (nlvl<0) nlvl=1.0;
	/* allocate memory for LAPACK */
  if ((bb=(double*)calloc((size_t)(N),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

	memcpy((void*)bb,(void*)b,(size_t)(N)*sizeof(double));

  if ((sigma=(double*)calloc((size_t)(N),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	vec_len=3*N*nlvl+11*(M);
	printf("nlvl=%d vec_len=%d\n",nlvl,vec_len);
  if ((IWORK=(int*)calloc((size_t)(vec_len),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}


	info=dgelsd(N,M,1,A,N,bb,N,sigma,2*dlamch('S'),&rank,w,-1,IWORK);
	if ( info ) {
	 fprintf(stderr,"something is wrong. flag=%d\n",info);
	}


	/*get actual work size */
	LWORK=(int)w[0];
	printf("work size =%d\n",LWORK);
  if ((WORK=(double*)calloc((size_t)(LWORK*2),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

	info=dgelsd(N,M,1,A,N,bb,N,sigma,2*dlamch('S'),&rank,WORK,LWORK,IWORK);
	if ( info ) {
	 fprintf(stderr,"something is wrong. flag=%d\n",info);
	}
	printf("optimal work size =%lf\n",WORK[0]);
	printf("rank=%d\n",rank);
#ifdef DEBUG
  printf("singular values\n");
  for(i=0;i<M; i++) {
		printf("%d: %lf\n",i,sigma[i]);
	}
  printf("solution\n");
  for(i=0;i<N; i++) {
		printf("%d: %lf\n",i,bb[i]);
	}
#endif
  /* copy solution back to x */
	memcpy((void*)x,(void*)bb,(size_t)(M)*sizeof(double));
	free(A);
	free(bb);
	free(sigma);
	free(WORK);
	free(IWORK);

	return info;
}
