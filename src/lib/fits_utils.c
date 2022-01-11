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


#include "shapelet.h"

#define DEBUG
/* Execute the command using this shell program.  */
#define SHELL "/bin/sh"
#include <unistd.h>
#include <sys/wait.h>
/* for fstat */
#include <stdio.h>
int
my_system (const char *command)
{
       int status;
       pid_t pid;

       pid = fork ();
       if (pid == 0)
         {
           /* This is the child process.  Execute the shell command. */
           execl (SHELL, SHELL, "-c", command, NULL);
           _exit (EXIT_FAILURE);
         }
       else if (pid < 0)
         /* The fork failed.  Report failure.  */
         status = -1;
       else
         /* This is the parent process.  Wait for the child to complete.  */
         if (waitpid (pid, &status, 0) != pid)
           status = -1;
       return status;
}

int zero_image_float(long totalrows, long offset, long firstrow, long nrows,
   int ncols, iteratorCol *cols, void *user_struct) {
    int ii;

    /* declare counts static to preserve the values between calls */
		/* so it traverses the whole array */
		static char *charp;
		static short int *sintp;
		static long int *lintp;
		static float *fptr;
		static double *dptr;

		static double tmpval;
		static long int pt,d1,d2,d3,d4;
    static long int xmin;
		static long int ymin,xmax,ymax;

    nlims *arr_dims=(nlims*)user_struct;
		int datatype=arr_dims->datatype;
    /*for (ii=0;ii<arr_dims->naxis;ii++) {
			printf("%d %ld\n",ii,arr_dims->d[ii]);
		}*/


    /*--------------------------------------------------------*/
    /*  Initialization procedures: execute on the first call  */
    /*--------------------------------------------------------*/
    if (firstrow == 1)
    {
       if (ncols != 1)
           return(-1);  /* number of columns incorrect */
#ifdef DEBUG
    for (ii=0;ii<arr_dims->naxis;ii++) {
			printf("%d %ld\n",ii,arr_dims->d[ii]);
		}
#endif
       /* assign the input pointers to the appropriate arrays and null ptrs*/
    switch (datatype) {
			case TBYTE:
       charp= (char *)  fits_iter_get_array(&cols[0]);
			 break;
		  case TSHORT:
       sintp= (short int*)  fits_iter_get_array(&cols[0]);
       break;
			case TLONG:
       lintp= (long int*)  fits_iter_get_array(&cols[0]);
       break;
			case TFLOAT:
       fptr= (float *)  fits_iter_get_array(&cols[0]);
			 break;
			case TDOUBLE:
       dptr= (double *)  fits_iter_get_array(&cols[0]);
			 
		}
			 /* initialize the limits */
			 xmin=arr_dims->d[0];
			 xmax=-1;
			 ymin=arr_dims->d[1];
			 ymax=-1;
    }

    //printf("limit (%ld,%ld)---(%ld,%ld)\n",xmin,ymin,xmax,ymax); 
    /*--------------------------------------------*/
    /*  Main loop: process all the rows of data */
    /*--------------------------------------------*/

    /*  NOTE: 1st element of array is the null pixel value!  */
    /*  Loop from 1 to nrows, not 0 to nrows - 1.  */

    switch (datatype) {
			case TBYTE:
         for (ii = 1; ii <= nrows; ii++) {
			     //printf("arr =%f\n",counts[ii]);
           //counts[ii] = 1.;
			     tmpval=(double)charp[ii];
			     if (arr_dims->tol<=fabs(tmpval)) {
			       //printf("arr =%lf\n",tmpval);
				     /* calculate 4D coords */
				     pt=firstrow+ii-1;
				     //printf("coord point=%ld ",pt);
				     d4=pt/(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2])*d4;
				     d3=pt/(arr_dims->d[0]*arr_dims->d[1]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1])*d3;
				     d2=pt/(arr_dims->d[0]);
				     pt-=(arr_dims->d[0])*d2;
				     d1=pt;
				     //printf("coords =(%ld,%ld,%ld,%ld)\n",d1,d2,d3,d4);
				     /* find current limit */
				     if (xmin>d1) {
						    xmin=d1;
				     }
				     if(xmax<d1) {
					      xmax=d1;
				     }
				     if (ymin>d2) {
						    ymin=d2;
				     }
				     if(ymax<d2) {
					      ymax=d2;
				     }
			
			    }

       }
			 break;

			case TSHORT:
         for (ii = 1; ii <= nrows; ii++) {
           //counts[ii] = 1.;
			     tmpval=(double)sintp[ii];
			     //printf("arr =%lf\n",tmpval);
			     if (arr_dims->tol<=fabs(tmpval)) {
			       //printf("arr =%lf\n",tmpval);
				     /* calculate 4D coords */
				     pt=firstrow+ii-1;
				     //printf("coord point=%ld ",pt);
				     d4=pt/(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2])*d4;
				     d3=pt/(arr_dims->d[0]*arr_dims->d[1]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1])*d3;
				     d2=pt/(arr_dims->d[0]);
				     pt-=(arr_dims->d[0])*d2;
				     d1=pt;
				     //printf("coords =(%ld,%ld,%ld,%ld)\n",d1,d2,d3,d4);
				     /* find current limit */
				     if (xmin>d1) {
						    xmin=d1;
				     }
				     if(xmax<d1) {
					      xmax=d1;
				     }
				     if (ymin>d2) {
						    ymin=d2;
				     }
				     if(ymax<d2) {
					      ymax=d2;
				     }
			
			    }

       }
			 break;

			case TLONG:
         for (ii = 1; ii <= nrows; ii++) {
			     //printf("arr =%f\n",counts[ii]);
           //counts[ii] = 1.;
			     tmpval=(double)lintp[ii];
			     if (arr_dims->tol<=fabs(tmpval)) {
			       //printf("arr =%lf\n",tmpval);
				     /* calculate 4D coords */
				     pt=firstrow+ii-1;
				     //printf("coord point=%ld ",pt);
				     d4=pt/(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2])*d4;
				     d3=pt/(arr_dims->d[0]*arr_dims->d[1]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1])*d3;
				     d2=pt/(arr_dims->d[0]);
				     pt-=(arr_dims->d[0])*d2;
				     d1=pt;
				     //printf("coords =(%ld,%ld,%ld,%ld)\n",d1,d2,d3,d4);
				     /* find current limit */
				     if (xmin>d1) {
						    xmin=d1;
				     }
				     if(xmax<d1) {
					      xmax=d1;
				     }
				     if (ymin>d2) {
						    ymin=d2;
				     }
				     if(ymax<d2) {
					      ymax=d2;
				     }
			
			    }

       }
			 break;

			case TFLOAT:
         for (ii = 1; ii <= nrows; ii++) {
			     //printf("arr =%f\n",counts[ii]);
           //counts[ii] = 1.;
			     tmpval=(double)fptr[ii];
			     if (arr_dims->tol<=fabs(tmpval)) {
			       //printf("arr =%lf\n",tmpval);
				     /* calculate 4D coords */
				     pt=firstrow+ii-1;
				     //printf("coord point=%ld ",pt);
				     d4=pt/(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2])*d4;
				     d3=pt/(arr_dims->d[0]*arr_dims->d[1]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1])*d3;
				     d2=pt/(arr_dims->d[0]);
				     pt-=(arr_dims->d[0])*d2;
				     d1=pt;
				     //printf("coords =(%ld,%ld,%ld,%ld)\n",d1,d2,d3,d4);
				     /* find current limit */
				     if (xmin>d1) {
						    xmin=d1;
				     }
				     if(xmax<d1) {
					      xmax=d1;
				     }
				     if (ymin>d2) {
						    ymin=d2;
				     }
				     if(ymax<d2) {
					      ymax=d2;
				     }
			
			    }

       }
			 break;

			case TDOUBLE:
         for (ii = 1; ii <= nrows; ii++) {
			     //printf("arr =%f\n",counts[ii]);
           //counts[ii] = 1.;
			     tmpval=(double)dptr[ii];
			     if (arr_dims->tol<=fabs(tmpval)) {
			       //printf("arr =%lf\n",tmpval);
				     /* calculate 4D coords */
				     pt=firstrow+ii-1;
				     //printf("coord point=%ld ",pt);
				     d4=pt/(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1]*arr_dims->d[2])*d4;
				     d3=pt/(arr_dims->d[0]*arr_dims->d[1]);
				     pt-=(arr_dims->d[0]*arr_dims->d[1])*d3;
				     d2=pt/(arr_dims->d[0]);
				     pt-=(arr_dims->d[0])*d2;
				     d1=pt;
				     //printf("coords =(%ld,%ld,%ld,%ld)\n",d1,d2,d3,d4);
				     /* find current limit */
				     if (xmin>d1) {
						    xmin=d1;
				     }
				     if(xmax<d1) {
					      xmax=d1;
				     }
				     if (ymin>d2) {
						    ymin=d2;
				     }
				     if(ymax<d2) {
					      ymax=d2;
				     }
			
			    }

       }
			 break;


		}
    //printf("cols =%d, starting row=%ld, nrows = %ld\n", ncols, firstrow, nrows);
    //printf("limit (%ld,%ld)---(%ld,%ld)\n",xmin,ymin,xmax,ymax); 
		/* set the current limit */
    arr_dims->lpix[0]=xmin;
    arr_dims->lpix[1]=ymin;
    arr_dims->lpix[2]=0;
    arr_dims->lpix[3]=0;
    arr_dims->hpix[0]=xmax;
    arr_dims->hpix[1]=ymax;
    arr_dims->hpix[2]=arr_dims->d[2]-1;
    arr_dims->hpix[3]=arr_dims->d[3]-1;
 
    return(0);  /* return successful status */

}

/* function to read min-max values of fits file */
int get_min_max(long totalrows, long offset, long firstrow, long nrows,
   int ncols, iteratorCol *cols, void *user_struct) {

		static double min_val;
		static double max_val;
		static double tmpval;
		int ii;

    drange *xylims=(drange *)user_struct;
		static char *charp;
		static short int *sintp;
		static long int *lintp;
		static float *fptr;
		static double *dptr;
	  int		datatype=xylims->datatype;


    if (firstrow == 1)
    {
       if (ncols != 1)
           return(-1);  /* number of columns incorrect */
       /* assign the input pointers to the appropriate arrays and null ptrs*/
    switch (datatype) {
			case TBYTE:
       charp= (char *)  fits_iter_get_array(&cols[0]);
			 break;
		  case TSHORT:
       sintp= (short int*)  fits_iter_get_array(&cols[0]);
       break;
			case TLONG:
       lintp= (long int*)  fits_iter_get_array(&cols[0]);
       break;
			case TFLOAT:
       fptr= (float *)  fits_iter_get_array(&cols[0]);
			 break;
			case TDOUBLE:
       dptr= (double *)  fits_iter_get_array(&cols[0]);
			 
		}

		min_val=1e6;
		max_val=-1e6;
    }

   switch (datatype) {
			case TBYTE:
        for (ii = 1; ii <= nrows; ii++) {
			    tmpval=(double)charp[ii];
		      if (min_val>tmpval) min_val=tmpval;
		      if (max_val<tmpval) max_val=tmpval;
        }
				break;
		  case TSHORT:
        for (ii = 1; ii <= nrows; ii++) {
			    tmpval=(double)sintp[ii];
					//printf("%lf==%d\n",tmpval,sintp[ii]);
		      if (min_val>tmpval) min_val=tmpval;
		      if (max_val<tmpval) max_val=tmpval;
        }
				break;
			case TLONG:
        for (ii = 1; ii <= nrows; ii++) {
			    tmpval=(double)lintp[ii];
		      if (min_val>tmpval) min_val=tmpval;
		      if (max_val<tmpval) max_val=tmpval;
        }
				break;
			case TFLOAT:
        for (ii = 1; ii <= nrows; ii++) {
			    tmpval=(double)fptr[ii];
		      if (min_val>tmpval) min_val=tmpval;
		      if (max_val<tmpval) max_val=tmpval;
        }
				break;
			case TDOUBLE:
        for (ii = 1; ii <= nrows; ii++) {
			    tmpval=(double)dptr[ii];
		      if (min_val>tmpval) min_val=tmpval;
		      if (max_val<tmpval) max_val=tmpval;
        }
				break;
	 }

		xylims->lims[0]=min_val;
		xylims->lims[1]=max_val;
#ifdef DEBUG
		//printf("min_max: min=%lf max=%lf\n",min_val,max_val);
#endif
		return 0;
}
 


/* filename: file name
 * cutoff: cutoff to truncate the image, -1 to ignore this, 1.0 to read whole image
 * myarr: 4D data array of truncated image
 * new_naxis: array of dimensions of each axis
 * lgrid: grid points in l axis
 * mgrid: grid points in m axis
 * ignore_wcs: sometimes too small images screw up WCS, set 1 to use manual grid
   cen : position (RA (h,m,s) Dec (h,m,s) 
   xlow,xhigh,ylow,yhigh: if cutoff==-1, then will use these to read a sub image
 p,q: if not zero, will shift center to these pixels (p,q)
* clipmin,clipmax: if not zero, clip values out this range before decomposition
   use_mask: if 1, look for fitfile.MASK.fits and use it as a mask  NOTE: mask file is assumed to have same dimensions as the image fits file
   nmaskpix: total pixels in the mask
 */
int read_fits_file(const char *filename,double cutoff, double**myarr, long int *new_naxis, double **lgrid, double **mgrid, io_buff *fbuff, int ignore_wcs, position *cen, int xlow, int xhigh, int ylow, int yhigh,double p, double q, double clipmin, double clipmax, int use_mask, int *nmaskpix) {
    iteratorCol cols[3];  /* structure used by the iterator function */
    int n_cols;
    long rows_per_loop, offset;
		drange arr_limits;

    int status;
		int naxis;
		int bitpix;

		int ii,jj,kk;
		int datatype=0;
		long int totalpix;
		double bscale,bzero;
		long int increment[4]={1,1,1,1};
		int null_flag=0;
    double nullval=0.0;

    /* for reading the mask */
    char *maskname=0;
    fitsfile *maskptr=0;
    double *maskarr=0;
    FILE *maskd;

		/* stuctures from WCSLIB */
		char *header;
		int ncard,nreject,nwcs;
		int ncoord;
		double *pixelc, *imgc, *worldc, *phic, *thetac;
		int *statc;
		double l0,m0;
    double del_lm; 
    int create_grid=0;

		int stat[NWCSFIX];

    int ra_h,ra_m,dec_d,dec_m;
    double ra_s,dec_s;


    /*** for calculating the center pixel values **/
		double cpixelc[4], cimgc[4], cworldc[4], cphic[1], cthetac[1];
    int cstatc[1];
    
    
    *nmaskpix=0; /* no mask, 0 pixels */
		int use_wcs=1;
    if (ignore_wcs) use_wcs=0;
		
    status = 0; 

#ifdef DEBUG
    printf("File =%s\n",filename);
#endif
    fits_open_file(&fbuff->fptr, filename, READWRITE, &status); /* open file */

/* WCSLIB et al. */
		/* read FITS header */
		if ((status = fits_hdr2str(fbuff->fptr, 1, NULL, 0, &header, &ncard, &status))) {
		 fits_report_error(stderr, status);
		 return 1;
		}

#ifdef DEBUG
		printf("header %s\n",header); 
#endif
/* try to Parse the primary header of the FITS file. */
    if ((status = wcspih(header, ncard, WCSHDR_all, 2, &nreject, &nwcs, &fbuff->wcs))) {
	      fprintf(stderr, "wcspih ERROR %d, ignoring WCS\n", status);
				use_wcs=0;
				/* to avoid compiler warnings do this : */
				pixelc=imgc=worldc=phic=thetac=0;
				statc=0;
				l0=m0=0;
		}

		status=0;
		if (use_wcs) {
		/* Fix non-standard WCS keyvalues. */
		if ((status = wcsfix(7, 0, fbuff->wcs, stat))) {
		  printf("wcsfix ERROR, status returns: (");
			  for (ii = 0; ii < NWCSFIX; ii++) {
					printf(ii ? ", %d" : "%d", stat[ii]);
				}
				printf(")\n\n");
	  }

		if ((status = wcsset(fbuff->wcs))) {
		  fprintf(stderr, "wcsset ERROR %d:\n", status);
		  return 1;
		}

#ifdef DEBUG
	  /* Print the struct. */
	  if ((status = wcsprt(fbuff->wcs))) return status;
#endif
		}

		/* turn off scaling so that we copy the pixel values */
    /* find exact scaling parameters */
    fits_read_key(fbuff->fptr,TDOUBLE,"BSCALE",&bscale,0,&status);
    if (status) { status=0; bscale=1.0; }
    fits_read_key(fbuff->fptr,TDOUBLE,"BZERO",&bzero,0,&status);
    if (status) { status=0; bzero=0.0; }
#ifdef DEBUG
		printf("Settin scale=%lf zero=%lf\n",bscale,bzero);
#endif
		//bscale=1.0; bzero=0.0;
    fits_set_bscale(fbuff->fptr,  bscale, bzero, &status);


		fits_get_img_dim(fbuff->fptr, &naxis, &status);
#ifdef DEBUG
		printf("Axis=%d\n",naxis);
#endif
    /* fix zero length axes */
		if (naxis<4) naxis=4;

		if ((fbuff->arr_dims.d=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((fbuff->arr_dims.lpix=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((fbuff->arr_dims.hpix=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		/* get axis dimensions */
		fits_get_img_size(fbuff->fptr, naxis, fbuff->arr_dims.d, &status);
    /* reset zero length axes to 1 */
    if (fbuff->arr_dims.d[2]==0) { fbuff->arr_dims.d[2]=1; }
    if (fbuff->arr_dims.d[3]==0) { fbuff->arr_dims.d[3]=1; }
		fbuff->arr_dims.naxis=naxis;
		fbuff->arr_dims.tol=cutoff;
		/* get data type */
		fits_get_img_type(fbuff->fptr, &bitpix, &status);
		if(bitpix==BYTE_IMG) {
#ifdef DEBUG
			printf("Type Bytes\n");
#endif
			datatype=TBYTE;
		}else if(bitpix==SHORT_IMG) {
#ifdef DEBUG
			printf("Type Short Int\n");
#endif
			datatype=TSHORT;
		}else if(bitpix==LONG_IMG) {
#ifdef DEBUG
			printf("Type Long Int\n");
#endif
			datatype=TLONG;
		}else if(bitpix==FLOAT_IMG) {
#ifdef DEBUG
			printf("Type Float\n");
#endif
			datatype=TFLOAT;
		}else if(bitpix==DOUBLE_IMG) {
#ifdef DEBUG
			printf("Type Double\n");
#endif
			datatype=TDOUBLE;
		}


    if (cutoff>0 && cutoff<1.0) {
    /* use the cutoff to find the rigt image size */
    n_cols = 1;

    /* define input column structure members for the iterator function */
    fits_iter_set_file(&cols[0], fbuff->fptr);
    fits_iter_set_iotype(&cols[0], InputOutputCol);
    fits_iter_set_datatype(&cols[0], 0);

    rows_per_loop = 0;  /* use default optimum number of rows */
    offset = 0;         /* process all the rows */


		/* determine limits of image data */
		arr_limits.datatype=datatype;
    fits_iterate_data(n_cols, cols, offset, rows_per_loop,
                      get_min_max, (void*)&arr_limits, &status);

#ifdef DEBUG
		printf("Limits Min %lf, Max %lf\n",arr_limits.lims[0],arr_limits.lims[1]);
#endif
		fbuff->arr_dims.tol=(1-cutoff)*(arr_limits.lims[1]-arr_limits.lims[0])+arr_limits.lims[0];
		/* need to transfer this to real value in the FITS file */
		/* using the inverse scaling , zero */
#ifdef DEBUG
		printf("cutoff for %lfx100 %% is %lf\n",cutoff,fbuff->arr_dims.tol);
#endif
    /* apply the rate function to each row of the table */
#ifdef DEBUG
    printf("Calling iterator function...%d\n", status);
#endif


    rows_per_loop = 0;  /* use default optimum number of rows */
    offset = 0;         /* process all the rows */
		fbuff->arr_dims.datatype=datatype;
    fits_iterate_data(n_cols, cols, offset, rows_per_loop,
                      zero_image_float, (void*)&fbuff->arr_dims, &status);

		/* sanity check: if no good pixels are found, include whole
		 * image 
		 */
		if(fbuff->arr_dims.lpix[0]==fbuff->arr_dims.d[0] || fbuff->arr_dims.lpix[1]==fbuff->arr_dims.d[1]
			||fbuff->arr_dims.hpix[0]==-1 || fbuff->arr_dims.hpix[1]==-1) {
					printf("No pixels found\n");
     fbuff->arr_dims.hpix[0]=fbuff->arr_dims.d[0]-1;
     fbuff->arr_dims.hpix[1]=fbuff->arr_dims.d[1]-1;
     fbuff->arr_dims.lpix[0]=fbuff->arr_dims.lpix[1]=0;
		}
    } else if (cutoff<0) { /* use the given grid*/
     fbuff->arr_dims.hpix[0]=xhigh-1;
     fbuff->arr_dims.hpix[1]=yhigh-1;
     fbuff->arr_dims.lpix[0]=xlow-1;
     fbuff->arr_dims.lpix[1]=ylow-1;
    } else {
      /* use whole image */
      fbuff->arr_dims.hpix[0]=fbuff->arr_dims.d[0]-1;
      fbuff->arr_dims.hpix[1]=fbuff->arr_dims.d[1]-1;
      fbuff->arr_dims.lpix[0]=fbuff->arr_dims.lpix[1]=0;
     }
 
		/* correct the coordinates for 1 indexing */
     fbuff->arr_dims.hpix[0]++;
     fbuff->arr_dims.hpix[1]++;
     fbuff->arr_dims.lpix[0]++;
     fbuff->arr_dims.lpix[1]++;

		 /* only work with stokes I */
     fbuff->arr_dims.hpix[2]=fbuff->arr_dims.hpix[3]=fbuff->arr_dims.lpix[2]=fbuff->arr_dims.lpix[3]=1;

#ifdef DEBUG
	  printf("(%ld %ld %ld %ld) ",fbuff->arr_dims.lpix[0],
									fbuff->arr_dims.lpix[1],	fbuff->arr_dims.lpix[2], fbuff->arr_dims.lpix[3]);
	  printf(" to (%ld %ld %ld %ld)\n",fbuff->arr_dims.hpix[0],
									fbuff->arr_dims.hpix[1],	fbuff->arr_dims.hpix[2], fbuff->arr_dims.hpix[3]);
#endif
	  /******* create new array **********/	
		new_naxis[0]=fbuff->arr_dims.hpix[0]-fbuff->arr_dims.lpix[0]+1;
		new_naxis[1]=fbuff->arr_dims.hpix[1]-fbuff->arr_dims.lpix[1]+1;
		new_naxis[2]=fbuff->arr_dims.hpix[2]-fbuff->arr_dims.lpix[2]+1;
		new_naxis[3]=fbuff->arr_dims.hpix[3]-fbuff->arr_dims.lpix[3]+1;
		/* calculate total number of pixels */
		totalpix=new_naxis[0]*new_naxis[1]*1*1; /* consider only one plane fron freq, and stokes axes because RA,Dec will not change */
#ifdef DEBUG
		printf("selecting %ld pixels\n",totalpix);
#endif
		
		if ((*myarr=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}


		/* read the subset increment=[1,1,1,..]*/
		fits_read_subset(fbuff->fptr, TDOUBLE, fbuff->arr_dims.lpix, fbuff->arr_dims.hpix, increment,
									 &nullval, *myarr, &null_flag, &status);

    /*****************************************************************/
    /* if use_mask>0, try to read mask file */
    if (use_mask>0) {
      if ((maskname=(char*)calloc(strlen(filename)+1+strlen(".MASK"),sizeof(char)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
      }
      strcpy(maskname,filename);
      /* if image is test.fits, mask file MUST be test.MASK.fits */
      ii=strlen(filename);
      maskname[ii-4]='M';
      maskname[ii-3]='A';
      maskname[ii-2]='S';
      maskname[ii-1]='K';
      maskname[ii]='.';
      maskname[ii+1]='f';
      maskname[ii+2]='i';
      maskname[ii+3]='t';
      maskname[ii+4]='s';
      maskname[ii+5]='\0';

     printf("MASK= %s\n",maskname);
      /* check if mask exists */
      maskd=0;
      maskd=fopen(maskname,"r");
      if (!maskd) {
       /* error in finding mask, so skip */
       fprintf(stderr,"%s: %d: no mask file %s found\n",__FILE__,__LINE__,maskname);
      } else {
       fclose(maskd);
       fits_open_file(&maskptr, maskname, READONLY, &status);
     	 if ((maskarr=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
			  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			  return 1;
		   }
		   /* read the subset increment=[1,1,1,..]*/
		   fits_read_subset(maskptr, TDOUBLE, fbuff->arr_dims.lpix, fbuff->arr_dims.hpix, increment,
			  &nullval, maskarr, &null_flag, &status);
       fits_close_file(maskptr, &status); 

      /* AND the mask with data */
      for (ii=0; ii<totalpix; ++ii) {
         (*myarr)[ii] *=maskarr[ii];
         /* also cound mask pixels >0 */
         if (maskarr[ii]>0) {
           *nmaskpix=*nmaskpix+1;
         }
      }
      
      free(maskarr);

      } 

      free(maskname);
    }   
    /*****************************************************************/

		/* ******************BEGIN create grid for the cells using WCS */
		if (use_wcs) {
#ifdef DEBUG
    printf("found axis %d using %d\n",fbuff->wcs->naxis,naxis);
#endif
		/* allocate memory for pixel/world coordinate arrays */
		ncoord=totalpix; /* consider only one plane fron freq, and stokes axes because RA,Dec will not change */
  	if ((pixelc=(double*)calloc((size_t)ncoord*4,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
  	if ((imgc=(double*)calloc((size_t)ncoord*4,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
  	if ((worldc=(double*)calloc((size_t)ncoord*4,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((phic=(double*)calloc((size_t)ncoord,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((thetac=(double*)calloc((size_t)ncoord,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((statc=(int*)calloc((size_t)ncoord,sizeof(int)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}

		/* fill up the pixel coordinate array */
    kk=0;
    for (ii=fbuff->arr_dims.lpix[0];ii<=fbuff->arr_dims.hpix[0];ii++)
     for (jj=fbuff->arr_dims.lpix[1];jj<=fbuff->arr_dims.hpix[1];jj++) {
						 pixelc[kk+0]=(double)ii;
						 pixelc[kk+1]=(double)jj;
						 pixelc[kk+2]=(double)1.0;
						 pixelc[kk+3]=(double)1.0;
						 kk+=4;
		 }
		/* now kk has passed the last pixel */
#ifdef DEBUG
		printf("total %d, created %d\n",ncoord,kk);
#endif
		if ((status = wcsp2s(fbuff->wcs, ncoord, naxis, pixelc, imgc, phic, thetac,
			 worldc, statc))) {
			 fprintf(stderr,"wcsp2s ERROR %2d\n", status);
			 /* Handle Invalid pixel coordinates. */
			 if (status == 8) status = 0;
	  }

    /* cleanup the array, clip values below threshold, also find the max value for center */
    if (clipmin!=clipmax && (clipmin!=0.0 || clipmax !=0.0)) {
    for (ii=0; ii<ncoord; ii++) {
       /* check for NaN and Inf values and remove them */
       if (isnan((*myarr)[ii])) {
          (*myarr)[ii]=0.0;
       } else if (isinf((*myarr)[ii])) {
          (*myarr)[ii]=0.0;
       } else {
        if ((*myarr)[ii] < clipmin) {
          (*myarr)[ii]=0.0; /*NOTE: put min value as 0.0 */
       }
        if ((*myarr)[ii]>clipmax) {
          (*myarr)[ii]=clipmax;
        }
       }
    }
    }

    /* the coordinate grid is in row major order, while the data is in column major order,
     so convert kk to row major order value */
//    jj=(kk+1)/new_naxis[0]-1;
//    ii=kk-(jj+1)*new_naxis[0]-1;
//    ii=new_naxis[0]-ii-1;
//    printf("PIX (%d,%d) axis [%ld,%ld]\n",ii,jj,new_naxis[0],new_naxis[1]);
//    /* recalculate kk in row major order */
//    kk=4*(ii*new_naxis[1]+jj);
//    printf("MAX pix=%d val=%lf\n",kk,maxval);
    if (p>0.0 || q>0.0) {
     /* shift center to pixels (p,q) */
    cpixelc[0]=p;
    cpixelc[1]=q;
    printf("Using center given (%lf,%lf)\n",p,q);
    cpixelc[2]=cpixelc[3]=1.0;
		if ((status = wcsp2s(fbuff->wcs, 1, naxis, cpixelc, cimgc, cphic, cthetac,
			 cworldc, cstatc))) {
			 fprintf(stderr,"wcsp2s ERROR %2d\n", status);
			 /* Handle Invalid pixel coordinates. */
			 if (status == 8) status = 0;
	  }
    l0=cimgc[0];
    m0=cimgc[1];

    cworldc[0]*=24.0/360.0;
    ra_h=(int)((int)cworldc[0]%24);
    ra_m=(int)((cworldc[0]-ra_h)*60.0);
    ra_s=(cworldc[0]-ra_h-ra_m/60.0)*3600.0;
    dec_d=(int)((int)cworldc[1]%180);
    dec_m=(int)((cworldc[1]-dec_d)*60.0);
    dec_s=(cworldc[1]-dec_d-dec_m/60.0)*3600.0;
    } else {
     /* no shift in pixel center, use origin */
 
    /* find center coordinates, handle even numbers correctly */
		/* even, use the pixel to the right as the center */
		/* odd, use middle pixel */
		/* find corresponding kk value =(l_c-1)*M +m_c-1 */
		 kk=4*((new_naxis[0]/2)*new_naxis[1]+(new_naxis[1]/2)); 

#ifdef DEBUG
		printf("found center %d\n",kk);
#endif
    /* find the phase centre in RA [0:360] ,Dec [0:180] (degrees) */
		l0=imgc[kk];
		m0=imgc[kk+1];
    /* RA DEC */
    worldc[kk]*=24.0/360.0;
    ra_h=(int)((int)worldc[kk]%24);
    ra_m=(int)((worldc[kk]-ra_h)*60.0);
    ra_s=(worldc[kk]-ra_h-ra_m/60.0)*3600.0;
    dec_d=(int)((int)worldc[kk+1]%180);
    dec_m=(int)((worldc[kk+1]-dec_d)*60.0);
    dec_s=(worldc[kk+1]-dec_d-dec_m/60.0)*3600.0;
 
    printf("RA DEC= %lf, %lf\n",worldc[kk],worldc[kk+1]);
   }

   /* if ra_h<0, change min sec to positive */
   if (ra_h<0) {
     ra_m=-ra_m;
     ra_s=-ra_s;
   }
   /* if dec_d<0, change min sec to positive */
   if (dec_d<0) {
     dec_m=-dec_m;
     dec_s=-dec_s;
   }
   printf("RA DEC= (%d:%d:%lf), (%d,%d,%lf)\n",ra_h,ra_m,ra_s,dec_d,dec_m,dec_s);
    cen->ra_h=ra_h;
    cen->ra_m=ra_m;
    cen->ra_s=ra_s;
    cen->dec_d=dec_d;
    cen->dec_m=dec_m;
    cen->dec_s=dec_s;

		}

		/* now calculate new (l,m) values for the phi,theta values with the new phase centre */

	if ((*lgrid=(double*)calloc((size_t)new_naxis[0],sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
	if ((*mgrid=(double*)calloc((size_t)new_naxis[1],sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}

   /* WCS can fail for very small images, so do a sanity check */
		if (use_wcs) {
		for (ii=0;ii<new_naxis[0];ii++) {
					//(*lgrid)[new_naxis[0]-ii-1]=(imgc[ii*4*new_naxis[1]]-l0)*M_PI/180.0;
					(*lgrid)[ii]=-(imgc[ii*4*new_naxis[1]]-l0)*M_PI/180.0;
		}
		/* different strides */
		for (ii=0;ii<new_naxis[1];ii++) {
					(*mgrid)[ii]=(imgc[ii*4+1]-m0)*M_PI/180.0;
		}

    /* sanity check the grid */
		for (ii=1;ii<new_naxis[0];ii++) {
      del_lm=(*lgrid)[ii]-(*lgrid)[ii-1];
      if (del_lm<=0) {
			  fprintf(stderr,"%s: %d: WCS L grid incorrect, creating by hand\n",__FILE__,__LINE__);
        create_grid=1;
        break;

      }
    }

    for (ii=1;ii<new_naxis[1];ii++) {
      del_lm=(*mgrid)[ii]-(*mgrid)[ii-1];
      if (del_lm<=0) {
			  fprintf(stderr,"%s: %d: WCS M grid incorrect\n",__FILE__,__LINE__);
        create_grid=1;
        break;
      }
    }

		}

/*
    for (ii=0;ii<new_naxis[0];ii++) {
     printf("L[%d]=%lf\n",ii,(*lgrid)[ii]);
    }
    for (ii=0;ii<new_naxis[1];ii++) {
     printf("M[%d]=%lf\n",ii,(*mgrid)[ii]);
    }
*/


    if (!use_wcs || create_grid) {
			/* no WCS, fill the grid by hand */
			/* in range -1 to +1 */
			l0=(double)new_naxis[0]/2;
			m0=(double)new_naxis[1]/2;
		  for (ii=0;ii<new_naxis[0];ii++) {
					(*lgrid)[ii]=-l0+ii;
		  }
	  	for (ii=0;ii<new_naxis[1];ii++) {
					(*mgrid)[ii]=-m0+ii;
	  	}
		}


    if (status) 
        fits_report_error(stderr, status);  /* print out error messages */

		free(header);
		
		if (use_wcs) {
		free(pixelc);
		free(imgc);
		free(worldc);
		free(phic);
		free(thetac);
		free(statc);
    }
    return(status);
}

/* filename: file name
 * myarr: 4D data array of truncated image
 * new_naxis: array of dimensions of each axis
 * lgrid: grid points in l axis
 * mgrid: grid points in m axis
 * ignore_wcs: sometimes too small images screw up WCS, set 1 to use manual grid
   cen : position (RA (h,m,s) Dec (h,m,s), input 
   Note: the grid will be shifted so that this will be the center
   xlow,xhigh,ylow,yhigh: will use these to read a sub image, if all not equal to zero
 */
extern int 
read_fits_file_recon(const char *filename,double**myarr, long int *new_naxis, double **lgrid, double **mgrid, io_buff *fbuff, int ignore_wcs, position cen, int xlow, int xhigh, int ylow, int yhigh){

    int status;
		int naxis;
		int bitpix;

		int ii,jj,kk;
		int datatype=0;
		long int totalpix;
		double bscale,bzero;
		long int increment[4]={1,1,1,1};
		int null_flag=0;

		/* stuctures from WCSLIB */
		char *header;
		int ncard,nreject,nwcs;
		//extern const char *wcshdr_errmsg[];
		int ncoord;
		double *pixelc, *imgc, *worldc, *phic, *thetac;
		int *statc;
		double l0,m0;
    double del_lm; 
    int create_grid=0;

		int stat[NWCSFIX];

    double ra0,dec0;
    
		int use_wcs=1;
    if (ignore_wcs) use_wcs=0;
		
    status = 0; 
#ifdef DEBUG
    printf("File =%s\n",filename);
#endif
    fits_open_file(&fbuff->fptr, filename, READWRITE, &status); /* open file */

/* WCSLIB et al. */
		/* read FITS header */
		if ((status = fits_hdr2str(fbuff->fptr, 1, NULL, 0, &header, &ncard, &status))) {
		 fits_report_error(stderr, status);
		 return 1;
		}

#ifdef DEBUG
		printf("header %s\n",header); 
#endif
/* try to Parse the primary header of the FITS file. */
    if ((status = wcspih(header, ncard, WCSHDR_all, 2, &nreject, &nwcs, &fbuff->wcs))) {
	      fprintf(stderr, "wcspih ERROR %d, ignoring WCS\n", status);
				use_wcs=0;
				/* to avoid compiler warnings do this : */
				pixelc=imgc=worldc=phic=thetac=0;
				statc=0;
				l0=m0=0;
		}

		status=0;
		if (use_wcs) {
		/* Fix non-standard WCS keyvalues. */
		if ((status = wcsfix(7, 0, fbuff->wcs, stat))) {
		  printf("wcsfix ERROR, status returns: (");
			  for (ii = 0; ii < NWCSFIX; ii++) {
					printf(ii ? ", %d" : "%d", stat[ii]);
				}
				printf(")\n\n");
	  }

		if ((status = wcsset(fbuff->wcs))) {
		  fprintf(stderr, "wcsset ERROR %d:\n", status);
		  return 1;
		}

#ifdef DEBUG
	  /* Print the struct. */
	  if ((status = wcsprt(fbuff->wcs))) return status;
#endif
		}

		/* find proper scaling if any */
    fits_read_key(fbuff->fptr,TDOUBLE,"BSCALE",&bscale,0,&status);
    if (status) { status=0; bscale=1.0; }
    fits_read_key(fbuff->fptr,TDOUBLE,"BZERO",&bzero,0,&status);
    if (status) { status=0; bzero=0.0; }
#ifdef DEBUG
		printf("Settin scale=%lf zero=%lf\n",bscale,bzero);
#endif
		//bscale=1.0; bzero=0.0;
    fits_set_bscale(fbuff->fptr,  bscale, bzero, &status);



		fits_get_img_dim(fbuff->fptr, &naxis, &status);
#ifdef DEBUG
		printf("Axis=%d\n",naxis);
#endif
    /* fix zero length axes */
		if (naxis<4) naxis=4;

		if ((fbuff->arr_dims.d=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((fbuff->arr_dims.lpix=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((fbuff->arr_dims.hpix=(long int*)calloc((size_t)naxis,sizeof(long int)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		/* get axis dimensions */
		fits_get_img_size(fbuff->fptr, naxis, fbuff->arr_dims.d, &status);
    /* reset zero length axes to 1 */
    if (fbuff->arr_dims.d[2]==0) { fbuff->arr_dims.d[2]=1; }
    if (fbuff->arr_dims.d[3]==0) { fbuff->arr_dims.d[3]=1; }
		fbuff->arr_dims.naxis=naxis;
		fbuff->arr_dims.tol=1.0;
		/* get data type */
		fits_get_img_type(fbuff->fptr, &bitpix, &status);
		if(bitpix==BYTE_IMG) {
#ifdef DEBUG
			printf("Type Bytes\n");
#endif
			datatype=TBYTE;
		}else if(bitpix==SHORT_IMG) {
#ifdef DEBUG
			printf("Type Short Int\n");
#endif
			datatype=TSHORT;
		}else if(bitpix==LONG_IMG) {
#ifdef DEBUG
			printf("Type Long Int\n");
#endif
			datatype=TLONG;
		}else if(bitpix==FLOAT_IMG) {
#ifdef DEBUG
			printf("Type Float\n");
#endif
			datatype=TFLOAT;
		}else if(bitpix==DOUBLE_IMG) {
#ifdef DEBUG
			printf("Type Double\n");
#endif
			datatype=TDOUBLE;
		}


    if (xlow || xhigh || ylow || yhigh) { /* use the given grid*/
     fbuff->arr_dims.hpix[0]=xhigh-1;
     fbuff->arr_dims.hpix[1]=yhigh-1;
     fbuff->arr_dims.lpix[0]=xlow-1;
     fbuff->arr_dims.lpix[1]=ylow-1;
    } else {
      /* use whole image */
      fbuff->arr_dims.hpix[0]=fbuff->arr_dims.d[0]-1;
      fbuff->arr_dims.hpix[1]=fbuff->arr_dims.d[1]-1;
      fbuff->arr_dims.lpix[0]=fbuff->arr_dims.lpix[1]=0;
     }
 
		/* correct the coordinates for 1 indexing */
     fbuff->arr_dims.hpix[0]++;
     fbuff->arr_dims.hpix[1]++;
     fbuff->arr_dims.hpix[2]++;
     fbuff->arr_dims.hpix[3]++;
     fbuff->arr_dims.lpix[0]++;
     fbuff->arr_dims.lpix[1]++;
     fbuff->arr_dims.lpix[2]++;
     fbuff->arr_dims.lpix[3]++;

		 /* only work with stokes I */
     fbuff->arr_dims.hpix[2]=fbuff->arr_dims.hpix[3]=fbuff->arr_dims.lpix[2]=fbuff->arr_dims.lpix[3]=1;

#ifdef DEBUG
	  printf("(%ld %ld %ld %ld) ",fbuff->arr_dims.lpix[0],
									fbuff->arr_dims.lpix[1],	fbuff->arr_dims.lpix[2], fbuff->arr_dims.lpix[3]);
	  printf(" to (%ld %ld %ld %ld)\n",fbuff->arr_dims.hpix[0],
									fbuff->arr_dims.hpix[1],	fbuff->arr_dims.hpix[2], fbuff->arr_dims.hpix[3]);
#endif
	  /******* create new array **********/	
		new_naxis[0]=fbuff->arr_dims.hpix[0]-fbuff->arr_dims.lpix[0]+1;
		new_naxis[1]=fbuff->arr_dims.hpix[1]-fbuff->arr_dims.lpix[1]+1;
		new_naxis[2]=fbuff->arr_dims.hpix[2]-fbuff->arr_dims.lpix[2]+1;
		new_naxis[3]=fbuff->arr_dims.hpix[3]-fbuff->arr_dims.lpix[3]+1;
		/* calculate total number of pixels */
    totalpix=((fbuff->arr_dims.hpix[0]-fbuff->arr_dims.lpix[0]+1)
     *(fbuff->arr_dims.hpix[1]-fbuff->arr_dims.lpix[1]+1)
     *(fbuff->arr_dims.hpix[2]-fbuff->arr_dims.lpix[2]+1)
     *(fbuff->arr_dims.hpix[3]-fbuff->arr_dims.lpix[3]+1));

#ifdef DEBUG
		printf("selecting %ld pixels\n",totalpix);
#endif
		
		if ((*myarr=(double*)calloc((size_t)totalpix,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}


		/* read the subset increment=[1,1,1,..]*/
		fits_read_subset(fbuff->fptr, TDOUBLE, fbuff->arr_dims.lpix, fbuff->arr_dims.hpix, increment,
									 0, *myarr, &null_flag, &status);

    

		/* ******************BEGIN create grid for the cells using WCS */
		if (use_wcs) {
#ifdef DEBUG
    printf("found axis %d using %d\n",fbuff->wcs->naxis,naxis);
#endif
		/* allocate memory for pixel/world coordinate arrays */
		ncoord=new_naxis[0]*new_naxis[1]*1*1; /* consider only one plane fron freq, and stokes axes because RA,Dec will not change */
  	if ((pixelc=(double*)calloc((size_t)ncoord*4,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
  	if ((imgc=(double*)calloc((size_t)ncoord*4,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
  	if ((worldc=(double*)calloc((size_t)ncoord*4,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((phic=(double*)calloc((size_t)ncoord,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((thetac=(double*)calloc((size_t)ncoord,sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
		if ((statc=(int*)calloc((size_t)ncoord,sizeof(int)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}

    /* calculate the RA,Dec to be shifted */
    ra0=(cen.ra_h+cen.ra_m/60.0+cen.ra_s/3600.0)*360/24.0;
    dec0=(cen.dec_d+cen.dec_m/60.0+cen.dec_s/3600.0);
    worldc[0]=ra0;
    worldc[1]=dec0;
    worldc[2]=1.0;
    worldc[3]=1.0;
    if ((status = wcss2p(fbuff->wcs, 1, naxis, worldc, phic, thetac, imgc, pixelc, statc))) {
       fprintf(stderr,"wcsp2s ERROR %2d\n", status);
       /* Handle Invalid pixel coordinates. */
       if (status == 8) status = 0;
    }
		l0=fbuff->arr_dims.d[0]/2-pixelc[0];
		m0=fbuff->arr_dims.d[1]/2-pixelc[1];
    /* somehow, RA pixel coords also shift like RA itself, right to left */
    if(xlow && xhigh) {
      fbuff->wcs->crpix[0]=xhigh-(pixelc[0]-xlow);
    } else { /* full image */ /* FIXME: get offset right */
      fbuff->wcs->crpix[0]=fbuff->arr_dims.d[0]-pixelc[0]+1;
    }
    fbuff->wcs->crpix[1]=pixelc[1];
    printf("Model RA DEC= %lf, %lf\n",worldc[0],worldc[1]);
    printf("Model Pixel l,m=(%lf,%lf)\n",pixelc[0],pixelc[1]);

 
		/* fill up the pixel coordinate array */
    kk=0;
    for (ii=fbuff->arr_dims.lpix[0];ii<=fbuff->arr_dims.hpix[0];ii++)
     for (jj=fbuff->arr_dims.lpix[1];jj<=fbuff->arr_dims.hpix[1];jj++) {
						 pixelc[kk+0]=(double)ii;//+l0;
						 pixelc[kk+1]=(double)jj;//+m0;
						 pixelc[kk+2]=(double)1.0;
						 pixelc[kk+3]=(double)1.0;
						 kk+=4;
		 }
		/* now kk has passed the last pixel */
#ifdef DEBUG
		printf("total %d, created %d\n",ncoord,kk);
#endif
		if ((status = wcsp2s(fbuff->wcs, ncoord, naxis, pixelc, imgc, phic, thetac,
			 worldc, statc))) {
			 fprintf(stderr,"wcsp2s ERROR %2d\n", status);
			 /* Handle Invalid pixel coordinates. */
			 if (status == 8) status = 0;
	  }

    printf("RA DEC= %lf, %lf\n",ra0,dec0);
    printf("l,m=(%lf,%lf)\n",l0,m0);

		}

		/* now calculate new (l,m) values for the phi,theta values with the new phase centre */

	if ((*lgrid=(double*)calloc((size_t)new_naxis[0],sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}
	if ((*mgrid=(double*)calloc((size_t)new_naxis[1],sizeof(double)))==0) {
			fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
			return 1;
		}

   /* WCS can fail for very small images, so do a sanity check */
		if (use_wcs) {
		for (ii=0;ii<new_naxis[0];ii++) {
					(*lgrid)[new_naxis[0]-ii-1]=(imgc[ii*4*new_naxis[1]])*M_PI/180.0;
		}
		/* different strides */
		for (ii=0;ii<new_naxis[1];ii++) {
					(*mgrid)[ii]=(imgc[ii*4+1])*M_PI/180.0;
		}

/*
		for (ii=0;ii<new_naxis[0];ii++) {
					printf("%d %lf\n",ii,(*lgrid)[ii]);
		}
		for (ii=0;ii<new_naxis[1];ii++) {
					printf("%d %lf\n",ii,(*mgrid)[ii]);
		}
*/

    /* sanity check the grid */
		for (ii=1;ii<new_naxis[0];ii++) {
      del_lm=(*lgrid)[ii]-(*lgrid)[ii-1];
      if (del_lm<=0) {
			  fprintf(stderr,"%s: %d: WCS L grid incorrect, creating by hand\n",__FILE__,__LINE__);
        create_grid=1;
        break;

      }
    }

    for (ii=1;ii<new_naxis[1];ii++) {
      del_lm=(*mgrid)[ii]-(*mgrid)[ii-1];
      if (del_lm<=0) {
			  fprintf(stderr,"%s: %d: WCS M grid incorrect\n",__FILE__,__LINE__);
        create_grid=1;
        break;
      }
    }

		}
    if (!use_wcs || create_grid) {
			/* no WCS, fill the grid by hand */
			/* in range -1 to +1 */
			l0=(double)new_naxis[0]/2;
			m0=(double)new_naxis[1]/2;
		  for (ii=0;ii<new_naxis[0];ii++) {
					(*lgrid)[ii]=-l0+ii;
		  }
	  	for (ii=0;ii<new_naxis[1];ii++) {
					(*mgrid)[ii]=-m0+ii;
	  	}
		}


    if (status) 
        fits_report_error(stderr, status);  /* print out error messages */

		free(header);
		
		if (use_wcs) {
		free(pixelc);
		free(imgc);
		free(worldc);
		free(phic);
		free(thetac);
		free(statc);
    }
    return(status);
}


/* simple read of FITS file, if the PSF is given as the FITS file*/
/* filename: fits file name
   psf: file will be stored in this Npx x Npy vector
*/
int 
read_fits_file_psf(const char *filename,double **psf, int *Npx, int *Npy){
  double *px, *py;
  position psf_cen;
  io_buff psf_filep;

  long int naxis[4]={0,0,0,0};

  px=py=0;
  /* read PSF fits file, ignore WCS utils (arg 7 is 1)*/
  read_fits_file_recon(filename, psf, naxis, &px, &py, &psf_filep, 1, psf_cen, 0,0,0,0);

  *Npx=naxis[0];
  *Npy=naxis[1];
  free(px);
  free(py);
  close_fits_file(psf_filep);

  return 0;
}

int close_fits_file(io_buff fbuff) {
   int status=0;


//		fits_write_subset(fbuff.fptr, TDOUBLE, fbuff.arr_dims.lpix, fbuff.arr_dims.hpix, myarr, &status);



    fits_close_file(fbuff.fptr, &status);      /* all done */

    if (status) 
        fits_report_error(stderr, status);  /* print out error messages */

		free(fbuff.arr_dims.d);
		free(fbuff.arr_dims.lpix);
		free(fbuff.arr_dims.hpix);
		wcsfree(fbuff.wcs);
    free(fbuff.wcs);
	
		return(status);
}


/* filename: original FITS file name
 * cutoff: cutoff to truncate the image
 * newfile:  new FITS file name
 * myarr: 4D data array tobe written to newfile
 */
int write_fits_file1(const char *newfile, double *myarr, io_buff fbuff) {
    fitsfile *newfptr;

    int status;
		int naxis;
		int bitpix,nkeys,ii;

    char card[81];

		double bscale,bzero;
    
    status = 0; 
#ifdef DEBUG
    printf("NEW File =%s\n",newfile);
#endif

    fits_create_file(&newfptr,newfile, &status);
    bitpix=DOUBLE_IMG;
    naxis=4;
    fits_create_img(newfptr, bitpix, naxis, fbuff.arr_dims.d, &status);
    /* copy all the user keywords (not the structural keywords) */
    fits_get_hdrspace(fbuff.fptr, &nkeys, NULL, &status);

    for (ii = 1; ii <= nkeys; ii++) {
      fits_read_record(fbuff.fptr, ii, card, &status);
      if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
         fits_write_record(newfptr, card, &status);
    }

   //fits_copy_file(fbuff.fptr,newfptr,1,1,1,&status);
   fits_copy_hdu(fbuff.fptr,newfptr,0,&status);

		bscale=1.0; bzero=0.0;
    fits_set_bscale(newfptr,  bscale, bzero, &status);
 		fits_write_subset(newfptr, TDOUBLE, fbuff.arr_dims.lpix, fbuff.arr_dims.hpix, myarr, &status);

    fits_close_file(newfptr, &status);      /* all done */

    if (status)
        fits_report_error(stderr, status);  /* print out error messages */

    
    return(status);
}


int write_fits_file(const char *oldfile,const char *newfile, double *myarr, io_buff fbuff) {
    fitsfile *newfptr;

    int status;
		int naxis;
		int bitpix,ii;
    long int naxes[4]={1,1,1,1};

    char *command;

		double bscale,bzero;
    
    status = 0; 
#ifdef DEBUG
    printf("NEW File =%s\n",newfile);
#endif
    if ((command=(char*)calloc(strlen("cp -f ")+strlen(oldfile)+strlen(" ")+strlen(newfile)+1,sizeof(char)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
    }

    strcpy(command,"cp -f ");
    strcat(command,oldfile);
    strcat(command," ");
    strcat(command,newfile);
    printf("running %s\n",command);
    ii=my_system(command);
    if (ii<0) {
      fprintf(stderr,"%s: %d: exec failed\n",__FILE__,__LINE__);
      return 1;
    }
    free(command);

    fits_open_file(&newfptr, newfile, READWRITE, &status); /* open file */
    fits_get_img_param(newfptr, 4, &bitpix, &naxis, naxes, &status);

		bscale=1.0; bzero=0.0;
    fits_set_bscale(newfptr,  bscale, bzero, &status);

 		fits_write_subset(newfptr, TDOUBLE, fbuff.arr_dims.lpix, fbuff.arr_dims.hpix, myarr, &status);


    fits_close_file(newfptr, &status);      /* all done */

    if (status)
        fits_report_error(stderr, status);  /* print out error messages */

    
    return(status);
}
