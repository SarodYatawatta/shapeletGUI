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
#include <pthread.h>
#include <locale.h>

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

/** calculate the mode vectors
 * in: x,y: arrays of the grid points
 * out:        
 *      M: number of modes
 *      beta: scale factor
 *      Av: array of mode vectors size Nx.Ny times n0.n0, in column major order
 *      n0: number of modes in each dimension
 *
 */
int 
calculate_mode_vectors(double *x, int Nx, double *y, int Ny, int *M, double 
								*beta, double **Av, int *n0) {

  /*a little quirk: the shapelet decom was done on a RA,Dec grid. In the l,m 
	grid the sign of l is reversed. So we need to change the coefficients to 
	the l,m grid. This means each coeff has to be multiplied bu (-1)^n1
	where n1 is the order of the basis functio in l axis
  FIXME: check this!
	*/
	double *grid;
	int *xindex,*yindex;
	int xci,yci,zci,Ntot;
	int *neg_grid;
	double th_min,th_max;

	double **shpvl, *fact;
	int n1,n2,start;
	//double pi_4r=1/sqrt(sqrt(M_PI));
	double xval;
	double tmpval;

  /* image size: Nx by Ny pixels
   */
	/* allocate memory to store grid points */
  if ((grid=(double*)calloc((size_t)(Nx+Ny),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

  if ((xindex=(int*)calloc((size_t)(Nx),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  if ((yindex=(int*)calloc((size_t)(Ny),sizeof(int)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

	/* merge axes coordinates */
	xci=yci=zci=0;
	while(xci<Nx && yci<Ny ) {
	 if (x[xci]==y[yci]){
		/* common index */
		grid[zci]=x[xci];
		xindex[xci]=zci;
		yindex[yci]=zci;
	  zci++;
	  xci++;	 
	  yci++;	 
	 } else if (x[xci]<y[yci]){
		 grid[zci]=x[xci];
		 xindex[xci]=zci;
	   zci++;
	   xci++;	 
	 } else {
		 grid[zci]=y[yci];
		 yindex[yci]=zci;
	   zci++;
	   yci++;	 
	 }	 
	}
	/* copy the tail */
	if(xci<Nx && yci==Ny) {
		/* tail from x */
		while(xci<Nx) {
		 grid[zci]=x[xci];
		 xindex[xci]=zci;
	   zci++;
	   xci++;	 
		}
	} else if (xci==Nx && yci<Ny) {
		/* tail from y */
		while(yci<Ny) {
		 grid[zci]=y[yci];
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
	for (xci=0; xci<Nx; xci++) {
	 printf("[%d]=%le ",xci,x[xci]);
	}
	printf("\n");
	for (xci=0; xci<Ny; xci++) {
	 printf("[%d]=%le ",xci,y[xci]);
	}
	printf("\n");
	for (xci=0; xci<Ntot; xci++) {
	 printf("[%d]=%lf ",xci,grid[xci]);
	}
	printf("\n");
	for (xci=0; xci<Nx; xci++) {
	 printf("[%d]=%d ",xci,xindex[xci]);
	}
	printf("\n");
	for (xci=0; xci<Ny; xci++) {
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
	while(xci<Ntot && grid[xci]<0.0) {
	 /* try to find a positive value for this is possible */
	 while(zci>=0 && grid[zci]>0.0 && -grid[xci]<grid[zci]) {
				zci--;
	 }
	 /* now we might have found a correct value */
	 if (zci>=0 && grid[zci]>0.0 && -grid[xci]==grid[zci]) {
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

	/* determine the scales and the number of modes */
	th_min=MAX_TOL;
	/* search the grid for smallest (not zero) scale */
  for (xci=1; xci<Ntot; xci++) {
	 tmpval=fabs(grid[xci]-grid[xci-1]);
	 if (tmpval>TOL && th_min>tmpval) th_min=tmpval;
	}

  /* beta ~ lambda/R, where lambda=wavelength, R= max UV distance in m */
	th_max=fabs(grid[Ntot-1]-grid[0]); /* largest scale */
	if (*beta==-1.0) {
	 *beta=sqrt(th_min*th_max);
	} /* else use given value */ 
	if (*M==-1) {
   if (*n0 >0) {
    *M=(*n0)*(*n0);
   } else {
  printf("Theta max,min %lf,%lf\n",th_max,th_min);
	  *M=(int)(th_max/th_min)+1; /* number of modes needed */
	  /* individual modes */
	  *n0=(int)(sqrt((double)*M)+1);
    *M=(*n0)*(*n0);
    }
	} else { /* else use given value */
	 *n0=(int)(sqrt((double)*M)+1);
   *M=(*n0)*(*n0);
  }

#ifdef DEBUG
	printf("min, max dimensions =%lf, %lf scale=%lf, modes=%d <(%d)^2\n",th_min,th_max,*beta,*M,*n0);
#endif
	/* set up factorial array */
  if ((fact=(double*)calloc((size_t)(*n0),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  fact[0]=1.0;
	for (xci=1; xci<(*n0); xci++) {
		fact[xci]=((double)xci)*fact[xci-1];
	}

#ifdef DEBUG
	printf("Fact\n");
	for (xci=0; xci<(*n0); xci++) {
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
   if ((shpvl[xci]=(double*)calloc((size_t)(*n0),sizeof(double)))==0) {
	   fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	   exit(1);
	 }
	}

	/* start filling in the array from the positive values */
	for (zci=Ntot-1; zci>=0; zci--) {
    /* check to see if there are any positive values */
     if (neg_grid[zci] !=-1) {
			/* copy in the values from positive one, with appropriate sign change */
	     for (xci=0; xci<(*n0); xci++) {
				 shpvl[zci][xci]=(xci%2==1?-shpvl[neg_grid[zci]][xci]:shpvl[neg_grid[zci]][xci]);
			 }
		 } else {
	     for (xci=0; xci<(*n0); xci++) {
				/*take into account the scaling
				*/
				 xval=grid[zci]/(*beta);
				 //shpvl[zci][xci]=inv_beta*pi_4r*H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt((2<<xci)*fact[xci]);
				 shpvl[zci][xci]=H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt((double)(2<<xci)*fact[xci]);
		   }
		 }
	}


#ifdef DEBUG
	for (zci=0; zci<Ntot; zci++) {
		printf("%lf= ",grid[zci]);
	  for (xci=0; xci<(*n0); xci++) {
		  printf("%lf, ",shpvl[zci][xci]);
		}
		printf("\n");
	}
#endif

	/* now calculate the mode vectors */
	/* each vector is Nx x Ny length and there are n0*n0 of them */
  if ((*Av=(double*)calloc((size_t)(Nx*Ny*(*n0)*(*n0)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
	for (n2=0; n2<(*n0); n2++) {
	 for (n1=0; n1<(*n0); n1++) {
    /* fill in Nx*Ny*(zci) to Nx*Ny*(zci+1)-1 */
		start=Nx*Ny*(n2*(*n0)+n1);
	  for (yci=0; yci<Ny; yci++) {
#pragma GCC ivdep
	   for (xci=0; xci<Nx; xci++) {
        //(*Av)[start+Ny*xci+yci]=-1;
        (*Av)[start+Nx*yci+xci]=shpvl[xindex[xci]][n1]*shpvl[yindex[yci]][n2];
			}
		}
	 }
	}

#ifdef DEBUG
	printf("Matrix dimension=%d by %d\n",Nx*Ny,(*n0)*(*n0));
#endif
//#define DEBUG
#ifdef DEBUG
	for (n1=0; n1<(*n0); n1++) {
	 for (n2=0; n2<(*n0); n2++) {
    /* fill in Nx*Ny*(zci) to Nx*Ny*(zci+1)-1 */
		start=Nx*Ny*(n1*(*n0)+n2);
	  for (xci=0; xci<Nx; xci++) {
	    for (yci=0; yci<Ny; yci++) {
        printf("%lf ",(*Av)[start+Ny*xci+yci]);
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

/* calculate mode vectors for each (x,y) point given by the arrays x, y
 * of equal length.
 *
 * in: x,y: arrays of the grid points
 *      N: number of grid points
 *      beta: scale factor
 *      n0: number of modes in each dimension
 * out:        
 *      Av: array of mode vectors size N times n0.n0, in column major order, 
 *      i.e., first N values for 0-th mode, second N values for 1-st mode ...
 *
 */
int
calculate_mode_vectors_bi(double *x, double *y, int N,  double beta, int n0, double **Av) {

	double *grid;
	int *xindex,*yindex;
	int xci,yci,zci,Ntot;
	int *neg_grid;

	double **shpvl, *fact;
	int n1,n2,start;

  /* for sorting */
  coordval *cx_val,*cy_val; 

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
   cx_val[xci].val=x[xci];
   cy_val[xci].idx=xci;
   cy_val[xci].val=y[xci];
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
    xindex[cx_val[xci].idx]=zci;
    yindex[cy_val[yci].idx]=zci;
	  zci++;
	  xci++;	 
	  yci++;	 
	 } else if (cx_val[xci].val<cy_val[yci].val){
		 grid[zci]=cx_val[xci].val;
     xindex[cx_val[xci].idx]=zci;
	   zci++;
	   xci++;	 
	 } else {
		 grid[zci]=cy_val[yci].val;
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
     xindex[cx_val[xci].idx]=zci;
	   zci++;
	   xci++;	 
		}
	} else if (xci==N && yci<N) {
		/* tail from y */
		while(yci<N) {
		 grid[zci]=cy_val[yci].val;
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
	 printf("[%d]=%lf %lf ",xci,x[xci],y[xci]);
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
				 double xvalt=grid[zci]/(beta);
				 //shpvl[zci][xci]=inv_beta*pi_4r*H_e(xval,xci)*exp(-0.5*xval*xval)/sqrt((2<<xci)*fact[xci]);
				 shpvl[zci][xci]=H_e(xvalt,xci)*exp(-0.5*xvalt*xvalt)/sqrt((2<<xci)*fact[xci]);
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
    /* fill in N*(zci) to N*(zci+1)-1 */
		start=N*(n2*(n0)+n1);
#pragma GCC ivdep
	  for (xci=0; xci<N; xci++) {
      (*Av)[start+xci]=shpvl[xindex[xci]][n1]*shpvl[yindex[xci]][n2];
		}
	 }
	}

#ifdef DEBUG
	printf("%%Matrix dimension=%d by %d\n",N,(n0)*(n0));
  printf("A=[\n");
	for (xci=0; xci<N; xci++) {
	for (n1=0; n1<(n0); n1++) {
	 for (n2=0; n2<(n0); n2++) {
    /* N*(zci) to N*(zci+1)-1 */
		start=N*(n1*(n0)+n2);
    printf("%lf ",(*Av)[start+xci]);
    }
	 }
	 printf("\n");
	}
  printf("];\n");
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




/* calculate mode vectors for the regular grid given by the x,y arrays after 
 * performing the given linear transform. i.e.,
 * |X| =    |a  0| |cos(th)  sin(th)| { |x| - |p| }
 * |Y|      |0  b| |-sin(th) cos(th)| { |y|   |q| }
 * where the new coords are (X,Y)
 * a,b: scaling (keep a=1, and make b smaller: similar to major/minor axes)
 * theta: rotation
 * p,q: shift of origin, in pixels, not used here (grid is shifted)
 * Note that even when we have regular x,y grid, we get non regular X,Y
 * grid. so we have to evaluate each new X,Y point individually.
 *
 * in: 
 *     x,y: coordinate arrays
 *     n0: no. of modes in each dimension
 *     beta: scale
 *     a,b: scaling
 *     theta: rotation 
 *      (all in radians)
 * out:
 *     Av: size Nx*Ny time n0*n0 matrix in column major order
 */
int
calculate_mode_vectors_tf(double *x, int Nx, double *y, int Ny,
                   double a, double b, double theta,
                        int *n0, double *beta, double **Av) {


 double ct,st,tempx,tempy;
 double *xx,*yy;
 int kk,ci,cj;

 ct=cos(theta);
 st=sin(theta);
 /* allocate memory to store new grid points size Nx*Ny */
 if ((xx=(double*)calloc((size_t)(Nx*Ny),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
 }
 if ((yy=(double*)calloc((size_t)(Nx*Ny),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
 }

 kk=0;
 for(cj=0; cj<Ny; cj++) {
   tempy=y[cj];
   for (ci=0; ci<Nx; ci++) { 
    tempx=x[ci];
    xx[kk]=a*(tempx*ct+tempy*st);
    yy[kk]=b*(tempy*ct-tempx*st);
    kk++;
  }
 }
 if (*n0<0 || *beta<0.0) {
   /* find default values */
   tempx=1e6; /* min */
   tempy=-1e6; /* max */
   for (kk=0; kk<Nx*Ny; kk++) {
    if (xx[kk]<tempx) tempx=xx[kk];
    if (xx[kk]>tempy) tempy=xx[kk];
    if (yy[kk]<tempx) tempx=yy[kk];
    if (yy[kk]>tempy) tempy=yy[kk];
   }
 tempx=fabs(tempy-tempx); /* largest scale */
 tempy=fabs(x[1]-x[0]); /* smallest ?? */

 if (*n0<0)  *n0=(int)sqrt(tempx/tempy)+1;
 if (*beta<0.0)  *beta=sqrt(tempx*tempy);
 }
#ifdef DEBUG
 printf("beta =%lf n0=%d\n",*beta,*n0);
#endif
 calculate_mode_vectors_bi(xx, yy, Nx*Ny,  *beta, *n0, Av);
 free(xx);
 free(yy);
 return 0;
}



int
calculate_mode_vectors_simple(double *x, double *y, int N,  double beta, int n0, double **Av) {

  /* set up factorial array */
  double *fact;
  if ((fact=(double*)calloc((size_t)(n0),sizeof(double)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  fact[0]=1.0;
  for (int ci=1; ci<(n0); ci++) {
    fact[ci]=((double)ci)*fact[ci-1];
  }

  if ((*Av=(double*)calloc((size_t)(N*(n0)*(n0)),sizeof(double)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }

  for (int ci=0; ci<N; ci++) {
    double xx=x[ci]/beta;
    double yy=y[ci]/beta;
    int cj=0;
    for (int n2=0; n2<n0; n2++) {
      for (int n1=0; n1<n0; n1++) {
        (*Av)[ci+cj*N]=H_e(xx,n1)*exp(-0.5*xx*xx)/sqrt((double)(2<<n1)*fact[n1])
          *H_e(yy,n2)*exp(-0.5*yy*yy)/(sqrt((double)(2<<n2)*fact[n2]));
        cj++;
      }
    }
  }

  free(fact);
  return 0;
}


typedef struct thread_data_pix_t_ {
  double *x;
  double *y;
  double beta;
  int start,end,N,n0;
  double *A;
  double *fact;
} thread_data_pix_t;

static void*
calculate_modes_th(void *data) {
  thread_data_pix_t *t=(thread_data_pix_t*)data;
  for (int ci=t->start; ci<=t->end; ci++) {
    double xx=t->x[ci]/t->beta;
    double yy=t->y[ci]/t->beta;
    int cj=0;
    for (int n2=0; n2<t->n0; n2++) {
      for (int n1=0; n1<t->n0; n1++) {
        t->A[ci+cj*t->N]=H_e(xx,n1)*exp(-0.5*xx*xx)/sqrt((double)(2<<n1)*t->fact[n1])
          *H_e(yy,n2)*exp(-0.5*yy*yy)/(sqrt((double)(2<<n2)*t->fact[n2]));
        cj++;
      }
    }
  }

  return NULL;
}

/* multi threaded version */
int
calculate_mode_vectors_thread(double *x, double *y, int N,  double beta, int n0, double **Av, int Nt) {

  /* set up factorial array */
  double *fact;
  if ((fact=(double*)calloc((size_t)(n0),sizeof(double)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  fact[0]=1.0;
  for (int ci=1; ci<(n0); ci++) {
    fact[ci]=((double)ci)*fact[ci-1];
  }

  if ((*Av=(double*)calloc((size_t)(N*(n0)*(n0)),sizeof(double)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }

  /* divide N pixels into Nt subsets */
  int Nthb0,Nthb;
  pthread_attr_t attr;
  pthread_t *th_array;
  thread_data_pix_t *threaddata;
  /* setup threads */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);
  if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((threaddata=(thread_data_pix_t*)malloc((size_t)Nt*sizeof(thread_data_pix_t)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }

  /* calculate min values a thread can handle */
  Nthb0=(N+Nt-1)/Nt;
  /* iterate over threads, allocating indices per thread */
  int ci=0;
  int nth;
  for (nth=0;  nth<Nt && ci<N; nth++) {
    if (ci+Nthb0<N) {
     Nthb=Nthb0;
    } else {
     Nthb=N-ci;
    }
    threaddata[nth].start=ci;
    threaddata[nth].end=ci+Nthb-1;
    threaddata[nth].A=*Av;
    threaddata[nth].fact=fact;
    threaddata[nth].beta=beta;
    threaddata[nth].n0=n0;
    threaddata[nth].x=x;
    threaddata[nth].y=y;
    threaddata[nth].N=N;
    pthread_create(&th_array[nth],&attr,calculate_modes_th,(void*)(&threaddata[nth]));
    /* next baseline set */
    ci=ci+Nthb;
  }
  /* now wait for threads to finish */
  for(int nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

  pthread_attr_destroy(&attr);
  free(th_array);
  free(threaddata);

  free(fact);
  return 0;
}



int
save_decomposition(const char* filename, double beta, int n0, double *modes, position cen) {
  setlocale(LC_NUMERIC,"POSIX"); /* dont print ',' instead of decimal '.' */
  FILE *fp;
  int i;
  fp=fopen(filename,"w");

  /* write to file */
  /* first line RA (h,m,s) Dec (d,m,s) */
  fprintf(fp,"%d %d %lf %d %d %lf\n",cen.ra_h,cen.ra_m,cen.ra_s, cen.dec_d,cen.dec_m,cen.dec_s);
  /* first line: n0 beta scaled by 2pi*/
  fprintf(fp,"%d %le\n",n0,beta*M_PI*2.0);
  /* now each indiviaul mode, not ?? scaled by beta^2 */
  for (i=0; i<n0*n0; i++)
   fprintf(fp,"%d %le\n",i,modes[i]);
  /* save LT parameters in parsable format */
  fprintf(fp,"L %lf %lf %lf\n",1.0,1.0,0.0);
  /* last lines additional info, save info on any linear transform used */
  fprintf(fp,"#\n#\n");
  fprintf(fp,"#a=%lf b=%lf theta=%lf p=%lf q=%lf\n",1.0,1.0,0.0,0.0,0.0);
  /* save filename, original beta */
  fprintf(fp,"#file=%s beta=%lf\n",filename,beta);
  /* as help save full line to be included in sky model */
  /*  name h m s d m s I Q U V spectral_index RM extent_X(rad) extent_Y(rad) pos_angle(rad) freq0 */
  fprintf(fp,"# LSM format:\n");
  fprintf(fp,"## %s %d %d %lf %d %d %lf 1 0 0 0 0 0 %lf %lf %lf 1000000.0\n",filename,cen.ra_h,cen.ra_m,cen.ra_s, cen.dec_d,cen.dec_m,cen.dec_s,1.0,1.0,0.0);
  fclose(fp);
  return 0;
}
