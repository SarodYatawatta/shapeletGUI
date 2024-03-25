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
#include <fitsio.h>
#include <math.h>

#include "shapelet.h"

/* struct for sorting coordinates */
typedef struct coordval_{
 double val;
 int idx;
} coordval;

/* comparison */
static int
compare_coordinates(const void *a, const void *b) {
 const coordval *da=(const coordval *)a;
 const coordval *db=(const coordval *)b;

 return(da->val>=db->val?1:-1);
}

/** calculate the UV mode vectors (block version, needs at least 2 grid points)
 * in: u,v: arrays of the grid points in UV domain
 *      M: number of modes
 *      beta: scale factor
 *      n0: number of modes in each dimension
 * out:
 *      Av: array of mode vectors size Nu.Nv times n0.n0, in column major order
 *      cplx: array of integers, size n0*n0, if 1 this mode is imaginary, else real
 *
 */
int
calculate_uv_mode_vectors(double *u, int Nu, double *v, int Nv, double beta, int n0, double **Av, int **cplx) {

	double *grid;
	int *xindex,*yindex;
	int xci,yci,zci,Ntot;
	int *neg_grid;

	double **shpvl, *fact;
	int n1,n2,start;
	//double pi_4r=1/sqrt(sqrt(M_PI));
	double xval;
	int signval;

  /* image size: Nu by Nv pixels
   */
	/* allocate memory to store grid points */
  if ((grid=(double*)calloc((size_t)(Nu+Nv),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

  if ((xindex=(int*)calloc((size_t)(Nu),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  if ((yindex=(int*)calloc((size_t)(Nv),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

	/* merge axes coordinates */
	xci=yci=zci=0;
	while(xci<Nu && yci<Nv ) {
	 if (u[xci]==v[yci]){
		/* common index */
		grid[zci]=u[xci];
		xindex[xci]=zci;
		yindex[yci]=zci;
	  zci++;
	  xci++;	 
	  yci++;	 
	 } else if (u[xci]<v[yci]){
		 grid[zci]=u[xci];
		 xindex[xci]=zci;
	   zci++;
	   xci++;	 
	 } else {
		 grid[zci]=v[yci];
		 yindex[yci]=zci;
	   zci++;
	   yci++;	 
	 }	 
	}
	/* copy the tail */
	if(xci<Nu && yci==Nv) {
		/* tail from x */
		while(xci<Nu) {
		 grid[zci]=u[xci];
		 xindex[xci]=zci;
	   zci++;
	   xci++;	 
		}
	} else if (xci==Nu && yci<Nv) {
		/* tail from y */
		while(yci<Nv) {
		 grid[zci]=v[yci];
		 yindex[yci]=zci;
	   zci++;
	   yci++;	 
		}
	}
	Ntot=zci;

	if (Ntot<2) {
		fprintf(stderr,"Error: Need at least 2 grid points\n");
		exit(1);
	}
#ifdef DEBUG
	printf("\n");
	for (xci=0; xci<Nu; xci++) {
	 printf("[%d]=%lf ",xci,u[xci]);
	}
	printf("\n");
	for (xci=0; xci<Nv; xci++) {
	 printf("[%d]=%lf ",xci,v[xci]);
	}
	printf("\n");
	for (xci=0; xci<Ntot; xci++) {
	 printf("[%d]=%lf ",xci,grid[xci]);
	}
	printf("\n");
	for (xci=0; xci<Nu; xci++) {
	 printf("[%d]=%d ",xci,xindex[xci]);
	}
	printf("\n");
	for (xci=0; xci<Nv; xci++) {
	 printf("[%d]=%d ",xci,yindex[xci]);
	}
	printf("\n");
#endif
	/* wrap up negative values to positive ones if possible */
  if ((neg_grid=(int*)calloc((size_t)(Ntot),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	for (zci=0; zci<Ntot; zci++) {
			neg_grid[zci]=-1;
	}
	zci=Ntot-1;
	xci=0;
	/* find positive values to all negative values if possible */
	while(xci<Ntot && grid[xci]<0) {
	 /* try to find a positive value for this is possible */
	 while(zci>=0 && grid[zci]>0 && -grid[xci]<grid[zci]) {
				zci--;
	 }
	 /* now we might have found a correct value */
	 if (zci>=0 && grid[zci]>0 && -grid[xci]==grid[zci]) {
			neg_grid[xci]=zci;
	 }
	 xci++;
	}

#ifdef DEBUG
	for (xci=0; xci<Ntot; xci++) {
	 printf("[%d]=%d ",xci,neg_grid[xci]);
	}
	printf("\n");
#endif


	/* set up factorial array */
  if ((fact=(double*)calloc((size_t)(n0),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  fact[0]=1;
	for (xci=1; xci<(n0); xci++) {
		fact[xci]=((double)xci)*fact[xci-1];
	}

#ifdef DEBUG
	printf("Fact\n");
	for (xci=0; xci<(n0); xci++) {
	 printf("[%d]=%lf ",xci,fact[xci]);
	}
	printf("\n");
#endif

	/* setup array to store calculated shapelet value */
	/* need max storage Ntot x n0 */
  if ((shpvl=(double**)calloc((size_t)(Ntot),sizeof(double*)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	for (xci=0; xci<Ntot; xci++) {
   if ((shpvl[xci]=(double*)calloc((size_t)(n0),sizeof(double)))==0) {
	   fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	   exit(1);
	 }
	}

	/* start filling in the array from the positive values */
	for (zci=Ntot-1; zci>=0; zci--) {
    /* check to see if there are any positive values */
     if (neg_grid[zci] !=-1) {
			/* copy in the values from positive one, with appropriate sign change */
	     for (xci=0; xci<(n0); xci++) {
				 shpvl[zci][xci]=(xci%2==1?-shpvl[neg_grid[zci]][xci]:shpvl[neg_grid[zci]][xci]);
			 }
		 } else {
	     for (xci=0; xci<(n0); xci++) {
				/*take into account the scaling - in Fourier domain
				*/
				 xval=grid[zci]*beta;
				 //shpvl[zci][xci]=sqr_beta*pi_4r*H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt((2<<xci)*fact[xci]);
				 shpvl[zci][xci]=H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt(pow(2.0,(double)xci+1)*fact[xci]);
		   }
		 }
	}


#ifdef DEBUG
	for (zci=0; zci<Ntot; zci++) {
		printf("%lf= ",grid[zci]);
	  for (xci=0; xci<(n0); xci++) {
		  printf("%lf, ",shpvl[zci][xci]);
		}
		printf("\n");
	}
#endif

	/* now calculate the mode vectors */
	/* each vector is Nu x Nv length and there are n0*n0 of them */
  if ((*Av=(double*)calloc((size_t)(Nu*Nv*(n0)*(n0)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  if ((*cplx=(int*)calloc((size_t)((n0)*(n0)),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

  for (n2=0; n2<(n0); n2++) {
	 for (n1=0; n1<(n0); n1++) {
	  (*cplx)[n2*n0+n1]=((n1+n2)%2==0?0:1) /* even (real) or odd (imaginary)*/;
		/* sign */
		if ((*cplx)[n2*n0+n1]==0) {
			signval=(((n1+n2)/2)%2==0?1:-1);
		} else {
			signval=(((n1+n2-1)/2)%2==0?1:-1);
		}

    /* fill in Nu*Nv*(zci) to Nu*Nv*(zci+1)-1 */
		start=Nu*Nv*(n2*(n0)+n1);
		if (signval==-1) {
	  for (yci=0; yci<Nv; yci++) {
	    for (xci=0; xci<Nu; xci++) {
        //??FIXME(*Av)[start+Nv*yci+xci]=-shpvl[xindex[xci]][n1]*shpvl[yindex[yci]][n2];
        (*Av)[start+Nu*yci+xci]=-shpvl[xindex[xci]][n1]*shpvl[yindex[yci]][n2];
			}
		}
		} else {
	  for (yci=0; yci<Nv; yci++) {
	    for (xci=0; xci<Nu; xci++) {
        //??FIXME(*Av)[start+Nv*yci+xci]=shpvl[xindex[xci]][n1]*shpvl[yindex[yci]][n2];
        (*Av)[start+Nu*yci+xci]=shpvl[xindex[xci]][n1]*shpvl[yindex[yci]][n2];
			}
		}
		}
	 }
	}

#ifdef DEBUG
	printf("Matrix dimension=%d by %d\n",Nu*Nv,(n0)*(n0));
#endif
#ifdef DEBUG
	for (n1=0; n1<(n0); n1++) {
	 for (n2=0; n2<(n0); n2++) {
		start=Nu*Nv*(n1*(n0)+n2);
	  for (xci=0; xci<Nu; xci++) {
	    for (yci=0; yci<Nv; yci++) {
        printf("%lf ",(*Av)[start+Nv*xci+yci]);
			}
		}
		printf("\n");
	 }
	}
#endif
	free(grid);
	free(xindex);
	free(yindex);
	free(neg_grid);
	free(fact);
	for (xci=0; xci<Ntot; xci++) {
	 free(shpvl[xci]);
	}
	free(shpvl);

  return 0;
}

/** calculate the UV mode vectors (scalar version, only 1 point)
 * in: u,v: arrays of the grid points in UV domain
 *      M: number of modes
 *      beta: scale factor
 *      n0: number of modes in each dimension
 * out:
 *      Av: array of mode vectors size Nu.Nv times n0.n0, in column major order
 *      cplx: array of integers, size n0*n0, if 1 this mode is imaginary, else real
 *
 */
int
calculate_uv_mode_vectors_scalar(double *u, int Nu, double *v, int Nv, double beta, int n0, double **Av, int **cplx) {

	int xci,zci,Ntot;

	double **shpvl, *fact;
	int n1,n2,start;
	//double pi_4r=1/sqrt(sqrt(M_PI));
	double xval;
	int signval;

  /* image size: Nu by Nv pixels
   */
	if (Nu!=1 || Nv != 1) { 
	  fprintf(stderr,"%s: %d: U,V arrays must be length 1\n",__FILE__,__LINE__);
	  exit(1);
	}
  Ntot=2; /* u,v seperately */
	/* set up factorial array */
  if ((fact=(double*)calloc((size_t)(n0),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  fact[0]=1;
	for (xci=1; xci<(n0); xci++) {
		fact[xci]=((double)xci)*fact[xci-1];
	}

#ifdef DEBUG
	printf("Fact\n");
	for (xci=0; xci<(n0); xci++) {
	 printf("[%d]=%lf ",xci,fact[xci]);
	}
	printf("\n");
#endif

	/* setup array to store calculated shapelet value */
	/* need max storage Ntot x n0 */
  if ((shpvl=(double**)calloc((size_t)(Ntot),sizeof(double*)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	for (xci=0; xci<Ntot; xci++) {
   if ((shpvl[xci]=(double*)calloc((size_t)(n0),sizeof(double)))==0) {
	   fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	   exit(1);
	 }
	}


	/* start filling in the array from the positive values */
	zci=0;
	for (xci=0; xci<n0; xci++) {
		xval=u[0]*beta;
		//shpvl[zci][xci]=sqr_beta*pi_4r*H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt((2<<xci)*fact[xci]);
		shpvl[zci][xci]=H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt(pow(2.0,(double)xci+1)*fact[xci]);
	}
	zci=1;
	for (xci=0; xci<n0; xci++) {
		xval=v[0]*beta;
		//shpvl[zci][xci]=sqr_beta*pi_4r*H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt((2<<xci)*fact[xci]);
		shpvl[zci][xci]=H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt(pow(2.0,(double)xci+1)*fact[xci]);
	}


#ifdef DEBUG
	for (zci=0; zci<Ntot; zci++) {
	  for (xci=0; xci<(n0); xci++) {
		  printf("%lf, ",shpvl[zci][xci]);
		}
		printf("\n");
	}
#endif

	/* now calculate the mode vectors */
	/* each vector is Nu x Nv length and there are n0*n0 of them */
  if ((*Av=(double*)calloc((size_t)(Nu*Nv*(n0)*(n0)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  if ((*cplx=(int*)calloc((size_t)((n0)*(n0)),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

  for (n2=0; n2<(n0); n2++) {
	 for (n1=0; n1<(n0); n1++) {
	  (*cplx)[n2*n0+n1]=((n1+n2)%2==0?0:1) /* even (real) or odd (imaginary)*/;
		/* sign */
		if ((*cplx)[n2*n0+n1]==0) {
			signval=(((n1+n2)/2)%2==0?1:-1);
		} else {
			signval=(((n1+n2-1)/2)%2==0?1:-1);
		}

    /* fill in Nu*Nv*(zci) to Nu*Nv*(zci+1)-1 */
		start=Nu*Nv*(n2*(n0)+n1);
		if (signval==-1) {
        (*Av)[start]=-shpvl[0][n1]*shpvl[1][n2];
		} else {
        (*Av)[start]=shpvl[0][n1]*shpvl[1][n2];
		}
	 }
	}

#ifdef DEBUG
	printf("Matrix dimension=%d by %d\n",Nu*Nv,(n0)*(n0));
#endif
#ifdef DEBUG
	for (n1=0; n1<(n0); n1++) {
	 for (n2=0; n2<(n0); n2++) {
		start=Nu*Nv*(n1*(n0)+n2);
        printf("%lf ",(*Av)[start]);
		printf("\n");
	 }
	}
#endif
	free(fact);
	for (xci=0; xci<Ntot; xci++) {
	 free(shpvl[xci]);
	}
	free(shpvl);

  return 0;
}


/* calculate mode vectors for each (u,v) point given by the arrays u, v
 * of equal length.
 *
 * in: u,v: arrays of the grid points
 *      N: number of grid points
 *      beta: scale factor
 *      n0: number of modes in each dimension
 * out:        
 *      Av: array of mode vectors size N times n0.n0, in column major order
 *      cplx: array of integers, size n0*n0, if 1 this mode is imaginary, else real
 *
 */
int
calculate_uv_mode_vectors_bi(double *u, double *v, int N,  double beta, int n0, double **Av, int **cplx) {

	double *grid;
	int *xindex,*yindex;
	int xci,yci,zci,Ntot;
	int *neg_grid;

	double **shpvl, *fact;
	int n1,n2,start;
	double xval;

  /* for sorting */
  coordval *cx_val,*cy_val; 
	int signval;

  /* image size: N pixels
   */
	/* allocate memory to store grid points: max unique points are 2N */
  if ((grid=(double*)calloc((size_t)(2*N),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

  if ((xindex=(int*)calloc((size_t)(N),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  if ((yindex=(int*)calloc((size_t)(N),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

  /* sort coordinates */
  if ((cx_val=(coordval*)calloc((size_t)(N),sizeof(coordval)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  if ((cy_val=(coordval*)calloc((size_t)(N),sizeof(coordval)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

  for (xci=0; xci<N; xci++) {
   cx_val[xci].idx=xci;
   cx_val[xci].val=u[xci];
   cy_val[xci].idx=xci;
   cy_val[xci].val=v[xci];
  }

#ifdef DEBUG
  printf("Before sort id x, y\n");
  for (xci=0; xci<N; xci++) {
   printf("%d %lf %lf\n",xci,cx_val[xci].val,cy_val[xci].val);
  }
#endif
  qsort(cx_val,N,sizeof(coordval),compare_coordinates);
  qsort(cy_val,N,sizeof(coordval),compare_coordinates);
#ifdef DEBUG
  printf("After sort id x, y\n");
  for (xci=0; xci<N; xci++) {
   printf("%d %lf %lf\n",xci,cx_val[xci].val,cy_val[xci].val);
  }
#endif

	/* merge axes coordinates */
	xci=yci=zci=0;
	while(xci<N && yci<N ) {
	 if (cx_val[xci].val==cy_val[yci].val){
		/* common index */
		grid[zci]=cx_val[xci].val;
		/*xindex[xci]=zci;
		yindex[yci]=zci; */
    xindex[cx_val[xci].idx]=zci;
    yindex[cy_val[yci].idx]=zci;
	  zci++;
	  xci++;	 
	  yci++;	 
	 } else if (cx_val[xci].val<cy_val[yci].val){
		 grid[zci]=cx_val[xci].val;
		 //xindex[xci]=zci;
     xindex[cx_val[xci].idx]=zci;
	   zci++;
	   xci++;	 
	 } else {
		 grid[zci]=cy_val[yci].val;
		 //yindex[yci]=zci;
     yindex[cy_val[yci].idx]=zci;
	   zci++;
	   yci++;	 
	 }	 
	}
	/* copy the tail */
	if(xci<N && yci==N) {
		/* tail from x */
		while(xci<N) {
		 grid[zci]=cx_val[xci].val;
		 //xindex[xci]=zci;
     xindex[cx_val[xci].idx]=zci;
	   zci++;
	   xci++;	 
		}
	} else if (xci==N && yci<N) {
		/* tail from y */
		while(yci<N) {
		 grid[zci]=cy_val[yci].val;
		 //yindex[yci]=zci;
     yindex[cy_val[yci].idx]=zci;
	   zci++;
	   yci++;	 
		}
	}
	Ntot=zci;

	if (Ntot<2) {
		fprintf(stderr,"Error: Need at least 2 grid points\n");
		exit(1);
	}
#ifdef DEBUG
	printf("Input coord points\n");
	for (xci=0; xci<N; xci++) {
	 printf("[%d]=%lf %lf ",xci,u[xci],v[xci]);
	}
	printf("Grid\n");
	for (xci=0; xci<Ntot; xci++) {
	 printf("[%d]=%lf ",xci,grid[xci]);
	}
	printf("Index x,y\n");
	for (xci=0; xci<N; xci++) {
	 printf("[%d]=%d %d ",xci,xindex[xci],yindex[xci]);
	}
	printf("\n");
#endif
	/* wrap up negative values to positive ones if possible */
  if ((neg_grid=(int*)calloc((size_t)(Ntot),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	for (zci=0; zci<Ntot; zci++) {
			neg_grid[zci]=-1;
	}
	zci=Ntot-1;
	xci=0;
	/* find positive values to all negative values if possible */
	while(xci<Ntot && grid[xci]<0) {
	 /* try to find a positive value for this is possible */
	 while(zci>=0 && grid[zci]>0 && -grid[xci]<grid[zci]) {
				zci--;
	 }
	 /* now we might have found a correct value */
	 if (zci>=0 && grid[zci]>0 && -grid[xci]==grid[zci]) {
			neg_grid[xci]=zci;
	 }
	 xci++;
	}

#ifdef DEBUG
	printf("Neg grid\n");
	for (xci=0; xci<Ntot; xci++) {
	 printf("[%d]=%d ",xci,neg_grid[xci]);
	}
	printf("\n");
#endif


	/* set up factorial array */
  if ((fact=(double*)calloc((size_t)(n0),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  fact[0]=1;
	for (xci=1; xci<(n0); xci++) {
		fact[xci]=((double)xci)*fact[xci-1];
	}

#ifdef DEBUG
	printf("Factorials\n");
	for (xci=0; xci<(n0); xci++) {
	 printf("[%d]=%lf ",xci,fact[xci]);
	}
	printf("\n");
#endif

	/* setup array to store calculated shapelet value */
	/* need max storage Ntot x n0 */
  if ((shpvl=(double**)calloc((size_t)(Ntot),sizeof(double*)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	for (xci=0; xci<Ntot; xci++) {
   if ((shpvl[xci]=(double*)calloc((size_t)(n0),sizeof(double)))==0) {
	   fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	   exit(1);
	 }
	}

  if ((*cplx=(int*)calloc((size_t)((n0)*(n0)),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

	/* start filling in the array from the positive values */
	for (zci=Ntot-1; zci>=0; zci--) {
    /* check to see if there are any positive values */
     if (neg_grid[zci] !=-1) {
			/* copy in the values from positive one, with appropriate sign change */
	     for (xci=0; xci<(n0); xci++) {
				 shpvl[zci][xci]=(xci%2==1?-shpvl[neg_grid[zci]][xci]:shpvl[neg_grid[zci]][xci]);
			 }
		 } else {
	     for (xci=0; xci<(n0); xci++) {
				/*take into account the scaling
				*/
				 xval=grid[zci]*(beta);
				 //shpvl[zci][xci]=inv_beta*pi_4r*H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt((2<<xci)*fact[xci]);
				 shpvl[zci][xci]=H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt(pow(2.0,(double)xci+1)*fact[xci]);
		   }
		 }
	}


#ifdef DEBUG
  printf("x, shapelet val\n");
	for (zci=0; zci<Ntot; zci++) {
		printf("%lf= ",grid[zci]);
	  for (xci=0; xci<(n0); xci++) {
		  printf("%lf, ",shpvl[zci][xci]);
		}
		printf("\n");
	}
#endif

	/* now calculate the mode vectors */
	/* each vector is Nx x Ny length and there are n0*n0 of them */
  if ((*Av=(double*)calloc((size_t)(N*(n0)*(n0)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	for (n2=0; n2<(n0); n2++) {
	 for (n1=0; n1<(n0); n1++) {
	  (*cplx)[n2*n0+n1]=((n1+n2)%2==0?0:1) /* even (real) or odd (imaginary)*/;
		/* sign */
		if ((*cplx)[n2*n0+n1]==0) {
			signval=(((n1+n2)/2)%2==0?1:-1);
		} else {
			signval=(((n1+n2-1)/2)%2==0?1:-1);
		}
    /* fill in N*(zci) to N*(zci+1)-1 */
		start=N*(n2*(n0)+n1);
		if (signval==-1) {
	    for (xci=0; xci<N; xci++) {
        (*Av)[start+xci]=-shpvl[xindex[xci]][n1]*shpvl[yindex[xci]][n2];
			}
		} else {
	    for (xci=0; xci<N; xci++) {
        (*Av)[start+xci]=shpvl[xindex[xci]][n1]*shpvl[yindex[xci]][n2];
			}
		}

	 }
	}

#ifdef DEBUG
	printf("Matrix dimension=%d by %d\n",N,(n0)*(n0));
	for (n1=0; n1<(n0); n1++) {
	 for (n2=0; n2<(n0); n2++) {
    /* fill in N*N*(zci) to N*N*(zci+1)-1 */
		start=N*(n1*(n0)+n2);
	  for (xci=0; xci<N; xci++) {
        printf("%lf ",(*Av)[start+xci]);
		}
		printf("\n");
	 }
	}
#endif
	free(grid);
	free(xindex);
	free(yindex);
  free(cx_val);
  free(cy_val);
	free(neg_grid);
	free(fact);
	for (xci=0; xci<Ntot; xci++) {
	 free(shpvl[xci]);
	}
	free(shpvl);
  return 0;
}




/* calculate mode vectors for the regular grid given by the u,v arrays after 
 * performing the given linear transform. i.e.,
 * |X| =    |a  0| |cos(th)  sin(th)| { |u| - |p| }
 * |Y|      |0  b| |-sin(th) cos(th)| { |v|   |q| }
 * where the new coords are (X,Y)
 * a,b: scaling (keep a=1, and make b smaller: similar to major/minor axes)
 * theta: rotation
 * p,q: shift of origin, in pixels, not used here (grid is shifted)
 * Note that even when we have regular x,y grid, we get non regular X,Y
 * grid. so we have to evaluate each new X,Y point individually.
 *
 * in: 
 *     u,v: coordinate arrays
 *     n0: no. of modes in each dimension
 *     beta: scale
 *     a,b: scaling
 *     theta: rotation 
 *      (all in radians)
 *     p,q:shift of origin (pixels) - not used here
 * out:
 *     Av: size Nu*Nv time n0*n0 matrix in column major order
 *      cplx: array of integers, size n0*n0, if 1 this mode is imaginary, else real
 */
int
calculate_uv_mode_vectors_tf(double *u, int Nu, double *v, int Nv,
                   double a, double b, double theta, double p, double q,
                        int n0, double beta, double **Av, int **cplx) {


 double ct,st,tempx,tempy;
 double *xx,*yy;
 int kk,ci,cj;

 /* UV mode, theta+pi/2 */
 ct=cos(theta+M_PI/2.0);
 st=sin(theta+M_PI/2.0);

 printf("Linear tf %lf %lf %lf\n",theta,a,b);
 /* allocate memory to store new grid points size Nx*Ny */
 if ((xx=(double*)calloc((size_t)(Nu*Nv),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
 }
 if ((yy=(double*)calloc((size_t)(Nu*Nv),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
 }

 kk=0;
 for(cj=0; cj<Nv; cj++) {
   tempy=v[cj];
   for (ci=0; ci<Nu; ci++) { 
    tempx=u[ci];
    xx[kk]=a*(tempx*ct+tempy*st);
    yy[kk]=b*(-tempx*st+tempy*ct);
    kk++;
  }
 }
 calculate_uv_mode_vectors_bi(xx, yy, Nu*Nv,  beta, n0, Av, cplx);
 free(xx);
 free(yy);
 return 0;
}
