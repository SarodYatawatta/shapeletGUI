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

//#define DEBUG
/*
 * the bulk of the work is done here;
 * filename: input FITS file
 * cutoff: cutoff [0,1)
 * x,Nx: Nx by 1 grid of x axis
 * y,Ny: Ny by 1 grid of y axis
 * beta: scale (if -1, the value is calculated, else this value is used)
 * M: no of modes (if -1, the value is calculated, else this value is used)
 * n0: actual no of modes in each axis
 * p,q: shift center to this pixel vales (starting from 1), if zero, ignored
 * clipmin,clipmax: if not zero, clip values out this range before decomposition
 * b: image pixel values
 * av: decomposed coefficient array
 * z: reconstructed pixel values
 cen: position (RA,Dec) h,m,s, d, m ,s
   convolve_psf: if >0, will attempt to get the PSF from the FITS file and convolve modes with it before decomposition
 psf_filename: if not 0, will try to read the PSF from this FITS file

  use_mask: if 1, look for fitfile.MASK.fits and use it as a mask NOTE: mask file is assumed to have same dimensions as the image fits file
 */
int decompose_fits_file(char* filename, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double p, double q, double clipmin, double clipmax, double **b, double **av, double **z, position *cen, int convolve_psf, char *psf_filename,  int use_mask) {
	double *Av;

	int modes;

	long int naxis[4]={0,0,0,0};
	io_buff filep;
  double fits_bmaj,fits_bmin,fits_bpa;
  double freq;

  double deltax,deltay,pixels_beam,tflux;


	int i;
  char *newfilename=0;
  double *res;

#ifdef DEBUG
	FILE *dbg;
#endif

  double *psf;
  int Npx,Npy;
  int Nm; /* MASK pixels, if any */

	/* read fits file */
	read_fits_file(filename,cutoff,b,naxis,x, y, &filep,0, cen,0,0,0,0,p,q, clipmin, clipmax, use_mask, &Nm, &fits_bmaj, &fits_bmin, &fits_bpa, &deltax, &deltay, &freq);
	*Nx=naxis[0]; *Ny=naxis[1];


	/* sanity check for cutoff, M and beta */
	if (*M!=-1) {
		if (*M<4) {
      fprintf(stderr,"Overriding modes from %d to %d\n",*M,-1);
			*M=-1;
		}
	}
	if (*beta!=-1.0) {
		if (*beta<1e-6) {
      fprintf(stderr,"Overriding scale from %lf to %d\n",*beta,-1);
			*beta=-1.0;
		}
	}


  calculate_mode_vectors(*x, *Nx, *y, *Ny, M, beta, &Av, n0);

  /* calculate foot print of psf (pi*a*b)/(deltax*deltay) no of pixels*/
  if (fits_bmaj!=-1 && fits_bmin!=-1) {
   pixels_beam=M_PI*(fits_bmaj)*(fits_bmin)/(deltax*deltay);
   if (pixels_beam<0) { pixels_beam=-pixels_beam; }
  } else {
   pixels_beam=1.0;
  }
  printf("foot print %lf pix\n",pixels_beam);


  if (convolve_psf) {
    if (!psf_filename) {
     //fits_bmaj=fits_bmin=7e-5;
     if (fits_bmaj!=-1 && fits_bmin!=-1) {
       convolve_with_psf(*x,*Nx,*y,*Ny,Av,*n0,fits_bmaj,fits_bmin,fits_bpa);
     } else {
      fprintf(stderr,"Unable to extract PSF information from the FITS file. BMAJ, BMIN, BPA not found. Not convolving.\n");
     }
    } else { /* read PSF from psf_filename */
      printf("using PSF FITS file= (%s)\n",psf_filename);
      read_fits_file_psf(psf_filename,&psf, &Npx, &Npy);
      convolve_with_psf_fits(*x,*Nx,*y,*Ny,Av,*n0,psf,Npx,Npy);
      free(psf);
    }
  }

 /* re-normalize image so that total_flux/(total_pixels/pixels_per_beam)=total_flux/beam=1 */
 /* Note: use absolute value of flux, to handle -ve sources */
  /* FIXME: use BLAS for this */
  tflux=0.0;
  for (i=0; i<(*Nx)*(*Ny); ++i) {
    tflux+=fabs((*b)[i]);
  }
  printf("Normalizing for total flux=%lf\n",tflux);
  /* if mask is used, normalize with the mask pixels */
  /* also scale by 1/sqrt(beta) */
  if (use_mask && Nm>0) {
   tflux=1.0/tflux*((double)(Nm)/pixels_beam)*(1.0/sqrt(*beta));
  } else {
   tflux=1.0/tflux*((double)((*Nx)*(*Ny))/pixels_beam)*(1.0/sqrt(*beta));
  }
  /* FIXME: use BLAS for this */
  //for (i=0; i<(*Nx)*(*Ny); ++i) {
  //  (*b)[i]*=tflux;
  //}
  dscal((*Nx)*(*Ny),tflux,*b);

#ifndef DEBUG
	printf("Image dimension=%d by %d\n",*Nx,*Ny);
	printf("Matrix dimension=%d by %d\n",(*Nx)*(*Ny),(*n0)*(*n0));
	printf("scale=%lf, modes=%d <(%d)^2\n",*beta,*M,*n0);
#endif

	/* allocate memory for LAPACK */
  if ((*av=(double*)calloc((size_t)((*n0)*(*n0)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
#ifdef DEBUG
  for(i=0;i<*n0*(*n0); i++) {
		printf("%d: %lf\n",i,(*av)[i]);
	}
#endif

    lsq_lapack(Av,*b,*av, (*Nx)*(*Ny), (*n0)*(*n0));

#ifdef DEBUG
	printf("solution\n");
  for(i=0;i<*n0*(*n0); i++) {
		printf("%d: %lf\n",i,(*av)[i]);
	}
#endif
	/* reconstruct the image using the solution */
  if ((*z=(double*)calloc((size_t)(*Nx*(*Ny)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

	modes=(*n0)*(*n0);
#ifdef DEBUG
	printf("A=\n");
  for(i=0;i<*Nx*(*Ny);i++) {
		printf("%d: %lf\n",i,Av[i]);
	}
#endif
  for(i=0;i<modes; i++) {
		/* y=y+A(:,i)*a[i] */
    daxpy(*Nx*(*Ny), &(Av[i*(*Nx)*(*Ny)]), (*av)[i], *z);
	}
	
#ifdef DEBUG
	printf("solution\n");
  for(i=0;i<*Nx*(*Ny); i++) {
		printf("%d: %lf ?? %lf\n",i,(*z)[i],(*b)[i]);
	}
#endif

#ifdef DEBUG
	/* print debug to file */
	if ((dbg=fopen("debug","w"))==0) {
	  fprintf(stderr,"%s: %d: cannot open file\n",__FILE__,__LINE__);
	  exit(1);
	}
  for(i=0;i<*Nx*(*Ny); i++) {
		fprintf(dbg,"%d %lf  %lf\n",i,(*z)[i],(*b)[i]);
	}
	fclose(dbg);
  /* print mode coefficients to another file */
	if ((dbg=fopen("mode_debug","w"))==0) {
	  fprintf(stderr,"%s: %d: cannot open file\n",__FILE__,__LINE__);
	  exit(1);
	}
  for(i=0;i<(*n0)*(*n0); i++) {
		fprintf(dbg,"%d %lf\n",i,(*av)[i]);
	}
	fclose(dbg);
#endif


 /*also save the reconstructed image as a FITS file */
  if ((res=(double*)calloc((size_t)(*Nx*(*Ny)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}


  if ((newfilename=(char*)calloc(strlen(filename)+1+strlen(".residual"),sizeof(char)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  strcpy(newfilename,filename);
  /* add suffix */
  i=strlen(newfilename);
  newfilename[i-5]='.';
  newfilename[i-4]='r';
  newfilename[i-3]='e';
  newfilename[i-2]='s';
  newfilename[i-1]='i';
  newfilename[i-0]='d';
  newfilename[i+1]='u';
  newfilename[i+2]='a';
  newfilename[i+3]='l';
  newfilename[i+4]='.';
  newfilename[i+5]='f';
  newfilename[i+6]='i';
  newfilename[i+7]='t';
  newfilename[i+8]='s';
  newfilename[i+9]='\0';
  /* calculate residual */
  /* res=b */
  dcopy(*Nx*(*Ny),*b,res);
  /* res=res - z */
  daxpy(*Nx*(*Ny), *z, -1.0, res);
  /* scale back residual to match correct flux */
  dscal((*Nx)*(*Ny),1.0/tflux,res);
  write_fits_file(filename,newfilename,res,filep);

	
	/*close FITS file */
	close_fits_file(filep);

  free(res);
	free(Av);
  free(newfilename);
	return 0;
}

/*
 * read FITS file and calculate mode vectors, do not do a decomposition
 * cutoff: cutoff [0,1)
 * x,Nx: Nx by 1 grid of x axis
 * y,Ny: Ny by 1 grid of y axis
 * beta: scale (if -1, the value is calculated, else this value is used)
 * M: no of modes (if -1, the value is calculated, else this value is used)
 * n0: actual no of modes in each axis
 * p,q: shift center to this pixel vales (starting from 1), if zero, ignored
 * b: image pixel values
 * Av: matrix with columns as each mode vector: size (Nx . Ny) x (n0 . n0)
 cen: position (RA,Dec) h,m,s, d, m ,s

   convolve_psf: if >0, will attempt to get the PSF from the FITS file and convolve modes with it before decomposition
 psf_filename: if not 0, will try to read the PSF from this FITS file
 */
int calcmode_fits_file(char* filename, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double p, double q, double **b, double **Av, position *cen, int convolve_psf , char *psf_filename) {

	long int naxis[4]={0,0,0,0};
	io_buff filep;
  double fits_bmaj,fits_bmin,fits_bpa, deltax, deltay;
  double freq;

  double *psf;
  int Npx,Npy;
  int Nm;


	/* read fits file */
	read_fits_file(filename,cutoff,b,naxis,x, y, &filep,0, cen,0,0,0,0,p,q,0.0, 0.0, 0, &Nm, &fits_bmaj, &fits_bmin, &fits_bpa, &deltax, &deltay, &freq);
	*Nx=naxis[0]; *Ny=naxis[1];


	/* sanity check for cutoff, M and beta */
	if (*M!=-1) {
		if (*M<4) {
      fprintf(stderr,"overriding modes from %d to %d\n",*M,-1);
			*M=-1;
		}
	}
	if (*beta!=-1) {
		if (*beta<MIN_TOL) {
      fprintf(stderr,"overriding scale from %lf to %d\n",*beta,-1);
			*beta=-1;
		}
	}


  calculate_mode_vectors(*x, *Nx, *y, *Ny, M, beta, Av, n0);

  if (convolve_psf) {
   if (!psf_filename) {
    if (fits_bmaj!=-1 && fits_bmin!=-1) {
     convolve_with_psf(*x,*Nx,*y,*Ny,*Av,*n0,fits_bmaj,fits_bmin,fits_bpa);
    }
   } else { /* read PSF from psf_filename */
     printf("using PSF FITS file= (%s)\n",psf_filename);
     read_fits_file_psf(psf_filename,&psf, &Npx, &Npy);
     convolve_with_psf_fits(*x,*Nx,*y,*Ny,*Av,*n0,psf,Npx,Npy);
     free(psf);
   }
  }

#ifdef DEBUG
	printf("Image dimension=%d by %d\n",*Nx,*Ny);
	printf("Matrix dimension=%d by %d\n",*Nx*(*Ny),*n0*(*n0));
	printf("scale=%lf, modes=%d <(%d)^2\n",*beta,*M,*n0);
#endif

 close_fits_file(filep);
 return 0;
}



/*
 * decompose with the given linear transform applied.
 * filename: input FITS file
 * cutoff: cutoff [0,1)
 * x,Nx: Nx by 1 grid of x axis
 * y,Ny: Ny by 1 grid of y axis
 * beta: scale
 * M: no of modes
 * n0: actual no of modes in each axis
 * a0,b0,theta,p,q: parameters in linear transform
 * a0/b0 : major/minor axes (better to keep one at 1, because beta can also change)
 * theta: position angle
 * p,q: central pixel value (absolute, from image)
 * clipmin,clipmax: if not zero, clip values out this range before decomposition
 * b: image pixel values
 * av: decomposed coefficient array
 * z: reconstructed pixel values
 cen: position (RA,Dec) h,m,s, d, m ,s
   convolve_psf: if >0, will attempt to get the PSF from the FITS file and convolve modes with it before decomposition
 psf_filename: if not 0, will try to read the PSF from this FITS file

  use_mask: if 1, look for fitfile.MASK.fits and use it as a mask NOTE: mask file is assumed to have same dimensions as the image fits file
 */
int
decompose_fits_file_tf(char* filename, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double a0, double b0, double theta, double p, double q, double clipmin, double clipmax, double **b, double **av, double **z, position *cen, int convolve_psf , char *psf_filename, int use_mask) {
	double *Av;

	int modes;

	long int naxis[4]={0,0,0,0};
	io_buff filep;

	int i;
  char *newfilename=0;
  double *res;
  double fits_bmaj,fits_bmin,fits_bpa,deltax,deltay;
  double freq;

#ifdef DEBUG
	FILE *dbg;
#endif

  double *psf;
  int Npx,Npy;
  int Nm;


	/* read fits file */
	read_fits_file(filename,cutoff,b,naxis,x, y, &filep,0, cen,0,0,0,0,p,q,clipmin, clipmax, use_mask, &Nm,  &fits_bmaj, &fits_bmin, &fits_bpa, &deltax, &deltay, &freq);
	*Nx=naxis[0]; *Ny=naxis[1];


	/* sanity check for cutoff, M and beta */
	if (*M!=-1) {
		if (*M<4) {
      fprintf(stderr,"overriding modes from %d to %d\n",*M,-1);
			*M=-1;
      *n0=-1;
		} else { 
     /* calculate n0 */
     *n0=sqrt((double)*M)+1;
    }
	} else {
      *n0=-1;
  }

	if (*beta!=-1) {
		if (*beta<MIN_TOL) {
      fprintf(stderr,"overriding scale from %lf to %d\n",*beta,-1);
			*beta=-1;
		}
	}

  calculate_mode_vectors_tf(*x, *Nx, *y, *Ny, a0, b0, theta, n0, beta, &Av);

  if (convolve_psf) {
   if (!psf_filename) {
    if (fits_bmaj!=-1 && fits_bmin!=-1) {
     convolve_with_psf(*x,*Nx,*y,*Ny,Av,*n0,fits_bmaj,fits_bmin,fits_bpa);
    }
   } else { /* read PSF from psf_filename */
     printf("using PSF FITS file= (%s)\n",psf_filename);
     read_fits_file_psf(psf_filename,&psf, &Npx, &Npy);
     convolve_with_psf_fits(*x,*Nx,*y,*Ny,Av,*n0,psf,Npx,Npy);
     free(psf);
   }
  }



#ifndef DEBUG
	printf("Image dimension=%d by %d\n",*Nx,*Ny);
	printf("Matrix dimension=%d by %d\n",*Nx*(*Ny),*n0*(*n0));
	printf("scale=%lf, modes=%d <(%d)^2\n",*beta,*M,*n0);
#endif

	/* allocate memory for LAPACK */
  if ((*av=(double*)calloc((size_t)((*n0)*(*n0)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
#ifdef DEBUG
  for(i=0;i<*n0*(*n0); i++) {
		printf("%d: %lf\n",i,*av[i]);
	}
#endif

    lsq_lapack(Av,*b,*av, (*Nx)*(*Ny), (*n0)*(*n0));

#ifdef DEBUG
	printf("solution\n");
  for(i=0;i<*n0*(*n0); i++) {
		printf("%d: %lf\n",i,*av[i]);
	}
#endif
	/* reconstruct the image using the solution */
  if ((*z=(double*)calloc((size_t)(*Nx*(*Ny)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}

	modes=(*n0)*(*n0);
#ifdef DEBUG
	printf("A=\n");
  for(i=0;i<*Nx*(*Ny);i++) {
		printf("%d: %lf\n",i,Av[i]);
	}
#endif
  for(i=0;i<modes; i++) {
		/* y=y+A(:,i)*a[i] */
    daxpy(*Nx*(*Ny), &(Av[i*(*Nx)*(*Ny)]), (*av)[i], *z);
	}
	
#ifdef DEBUG
	printf("solution\n");
  for(i=0;i<*Nx*(*Ny); i++) {
		printf("%d: %lf ?? %lf\n",i,(*z)[i],(*b)[i]);
	}
#endif

#ifdef DEBUG
	/* print debug to file */
	if ((dbg=fopen("debug","w"))==0) {
	  fprintf(stderr,"%s: %d: cannot open file\n",__FILE__,__LINE__);
	  exit(1);
	}
  for(i=0;i<*Nx*(*Ny); i++) {
		fprintf(dbg,"%d %lf  %lf\n",i,(*z)[i],(*b)[i]);
	}
	fclose(dbg);
  /* print mode coefficients to another file */
	if ((dbg=fopen("mode_debug","w"))==0) {
	  fprintf(stderr,"%s: %d: cannot open file\n",__FILE__,__LINE__);
	  exit(1);
	}
  for(i=0;i<(*n0)*(*n0); i++) {
		fprintf(dbg,"%d %lf\n",i,(*av)[i]);
	}
	fclose(dbg);
#endif

 /*also save the reconstructed image as a FITS file */
  if ((res=(double*)calloc((size_t)(*Nx*(*Ny)),sizeof(double)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}


  if ((newfilename=(char*)calloc(strlen(filename)+1+strlen(".residual"),sizeof(char)))==0) {
	  fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
	  exit(1);
	}
  strcpy(newfilename,filename);
  /* add suffix */
  i=strlen(newfilename);
  newfilename[i-5]='.';
  newfilename[i-4]='r';
  newfilename[i-3]='e';
  newfilename[i-2]='s';
  newfilename[i-1]='i';
  newfilename[i-0]='d';
  newfilename[i+1]='u';
  newfilename[i+2]='a';
  newfilename[i+3]='l';
  newfilename[i+4]='.';
  newfilename[i+5]='f';
  newfilename[i+6]='i';
  newfilename[i+7]='t';
  newfilename[i+8]='s';
  newfilename[i+9]='\0';
  /* calculate residual */
  /* res=b */
  dcopy(*Nx*(*Ny),*b,res);
  /* res=res - z */
  daxpy(*Nx*(*Ny), *z, -1.0, res);
  write_fits_file(filename,newfilename,res,filep);

	
	/*close FITS file */
	close_fits_file(filep);

	free(Av);
	return 0;


}


/*
 * support linear transforms
 * read FITS file and calculate mode vectors, do not do a decomposition
 * cutoff: cutoff [0,1)
 * x,Nx: Nx by 1 grid of x axis
 * y,Ny: Ny by 1 grid of y axis
 * beta: scale (if -1, the value is calculated, else this value is used)
 * M: no of modes (if -1, the value is calculated, else this value is used)
 * n0: actual no of modes in each axis
 * b: image pixel values
 * Av: matrix with columns as each mode vector: size (Nx . Ny) x (n0 . n0)
 * a0,b0,theta,p,q: params of linear transform
 cen: position (RA,Dec) h,m,s, d, m ,s
 convolve_psf: if >0, will attempt to get the PSF from the FITS file and convolve modes with it before decomposition
 psf_filename: if not 0, will try to read the PSF from this FITS file
 */
extern int
calcmode_fits_file_tf(char* filename, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double a0, double b0, double theta, double p, double q, double **b, double **Av, position *cen, int convolve_psf , char *psf_filename) {
	long int naxis[4]={0,0,0,0};
	io_buff filep;

  double fits_bmaj,fits_bmin,fits_bpa,deltax,deltay,freq;
  int status=0;

  double *psf;
  int Npx,Npy;
  int Nm;


	/* read fits file */
	read_fits_file(filename,cutoff,b,naxis,x, y, &filep,0, cen,0,0,0,0,p,q, 0.0, 0.0, 0, &Nm, &fits_bmaj, &fits_bmin, &fits_bpa, &deltax, &deltay, &freq);
	*Nx=naxis[0]; *Ny=naxis[1];


  printf("Modes TF\n");
	/* sanity check for cutoff, M and beta */
	if (*M!=-1) {
		if (*M<4) {
      fprintf(stderr,"overriding modes from %d to %d\n",*M,-1);
			*M=-1;
		}
	}
	if (*beta!=-1) {
		if (*beta<1e-6) {
      fprintf(stderr,"overriding scale from %lf to %d\n",*beta,-1);
			*beta=-1;
		}
	}

  *n0=-1; 

  calculate_mode_vectors_tf(*x, *Nx, *y, *Ny, a0, b0, theta, n0, beta, Av);
  fits_bmaj=fits_bmin=-1; /* assume no key present */
  if (convolve_psf) {
   if (!psf_filename) {
    /* try to read psf params from file */
    fits_read_key(filep.fptr,TDOUBLE,"BMAJ",&fits_bmaj,0,&status);
    /* recover error from missing key */
    if (status) {
     status=0;
     fits_bmaj=fits_bmin=-1; /* no key present */
    } else {
     fits_read_key(filep.fptr,TDOUBLE,"BMIN",&fits_bmin,0,&status);
     if (status) {
      status=0;
      fits_bmaj=fits_bmin=-1; /* no key present */
     } else {
      fits_read_key(filep.fptr,TDOUBLE,"BPA",&fits_bpa,0,&status);
      if (status) {
       status=0;
       fits_bmaj=fits_bmin=-1; /* no key present */
      } else { /* convert to radians */
       fits_bpa=(fits_bpa)/180.0*M_PI;
       printf("beam= (%lf,%lf,%lf)\n",fits_bmaj,fits_bmin,fits_bpa);
       fits_bmaj=fits_bmaj/360.0*M_PI;
       fits_bmin=fits_bmin/360.0*M_PI;
     }
    }
    }
    if (fits_bmaj!=-1 && fits_bmin!=-1) {
     convolve_with_psf(*x,*Nx,*y,*Ny,*Av,*n0,fits_bmaj,fits_bmin,fits_bpa);
    }
   } else { /* read PSF from psf_filename */
     printf("using PSF FITS file= (%s)\n",psf_filename);
     read_fits_file_psf(psf_filename,&psf, &Npx, &Npy);
     convolve_with_psf_fits(*x,*Nx,*y,*Ny,*Av,*n0,psf,Npx,Npy);
     free(psf);
   }
  }




	close_fits_file(filep);
#ifdef DEBUG
	printf("Image dimension=%d by %d\n",*Nx,*Ny);
	printf("Matrix dimension=%d by %d\n",*Nx*(*Ny),*n0*(*n0));
	printf("scale=%lf, modes=%d <(%d)^2\n",*beta,*M,*n0);
#endif

 return 0;
}


