/*
   Copyright (C) 2022 Sarod Yatawatta <sarod@users.sf.net>  
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


#ifndef SHAPELET_H
#define SHAPELET_H
#ifdef __cplusplus
        extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <wcs.h>
#include <prj.h>
#include <wcshdr.h>
#include <wcsfix.h>

#ifndef TOL
#define TOL 1e-12
#endif
#ifndef MIN_TOL
#define MIN_TOL 1e-24
#endif
#ifndef MAX_TOL
#define MAX_TOL 1e+24
#endif


#ifndef ABS
#define ABS(x) \
  ((x)>=0? (x):-(x))
#endif

#ifndef IS_ZERO
#define IS_ZERO(x)\
        (ABS(x)<= TOL)
#endif

/**********************************************************
 * fits_utils.c
 **********************************************************/
typedef struct nlimits_ {
		int naxis;
		long int *d;
		long int *lpix;
		long int *hpix;
		double tol;
		int datatype;
} nlims; /* struct to store array dimensions */

typedef struct drange_ {
		double lims[2]; /* min, max values */
		int datatype;
} drange; /* struct to store min, max limits of image, needs the data type used in the pixels as well */

typedef struct __io_buff__ {
		fitsfile *fptr;
		nlims arr_dims;
		struct wcsprm *wcs;
} io_buff;


/* structure to store RA,Dec position */
typedef struct position__ {
  int ra_h;
  int ra_m;
  double ra_s;
  int dec_d;
  int dec_m;
  double dec_s;
} position; 

extern int 
zero_image_float(long totalrows, long offset, long firstrow, long nrows,
				int ncols, iteratorCol *cols, void *user_struct);
extern int 
get_min_max(long totalrows, long offset, long firstrow, long nrows,
				int ncols, iteratorCol *cols, void *user_struct);

/* filename: file name
 * cutoff: cutoff to truncate the image, -1 to ignore this, 1.0 to read whole image
 * myarr: 4D data array of truncated image (output)
 * new_naxis: array of dimensions of each axis (output)
 * lgrid: grid points in l axis (output)
 * mgrid: grid points in m axis (output)
 * fbuff: (output)
 * ignore_wcs: sometimes too small images screw up WCS, set 1 to use manual grid
   cen : position (RA (h,m,s) Dec (h,m,s)  (output)
   xlow,xhigh,ylow,yhigh: if cutoff==-1, then will use these to read a sub image
   p,q: if not zero, will shift center to these pixels (p,q)
* clipmin,clipmax: if not zero, clip values out this range before decomposition
  use_mask: if 1, look for fitfile.MASK.fits and use it as a mask NOTE: mask file is assumed to have same dimensions as the image fits file
  nmaskpix: total pixels in the mask (output)
  bmaj,bmin,bpa: PSF (beam in radians) if the header includes the info (output)
  deltax, deltay : pixel deltas in radians (output)
  freq: (central) frequency in Hz (output)
*/
extern int 
read_fits_file(const char *filename,double cutoff, double **myarr, long int *new_naxis, double **lgrid, double **mgrid, io_buff *fbuff, int ignore_wcs, position *cen, int xlow, int xhigh, int ylow, int yhigh, double p, double q, double clipmin, double clipmax, int use_mask, int *nmaskpix, double *bmaj, double *bmin, double *bpa, double *deltax, double *deltay, double *freq);

/* filename: file name
 * myarr: 4D data array of truncated image
 * new_naxis: array of dimensions of each axis
 * lgrid: grid points in l axis
 * mgrid: grid points in m axis
 * ignore_wcs: sometimes too small images screw up WCS, set 1 to use manual grid
cen : position (RA (h,m,s) Dec (h,m,s), Note: input 
   Note: the grid will be shifted so that this will be the center
   xlow,xhigh,ylow,yhigh: will use these to read a sub image, if all not equal to zero
 */
extern int 
read_fits_file_recon(const char *filename,double**myarr, long int *new_naxis, double **lgrid, double **mgrid, io_buff *fbuff, int ignore_wcs, position cen, int xlow, int xhigh, int ylow, int yhigh);


/* 
 *   newfile:  new FITS file name
 *   myarr: 4D data array tobe written to newfile
 *     */
extern int 
write_fits_file(const char *oldfine, const char *newfile, double *myarr, io_buff fbuff);

extern int 
close_fits_file(io_buff fbuff);

/* simple read of FITS file, if the PSF is given as the FITS file*/
/* filename: fits file name
   psf: file will be stored in this Npx x Npy vector
*/
extern int 
read_fits_file_psf(const char *filename,double **psf, int *Npx, int *Npy);


/****************************************
 * hermite.c
 ****************************************/
/* evaluate Hermite polynomial value using recursion
 */
extern double
H_e(double x, int n);

/********************************************
 * shapelet_lm.c
 *******************************************/
/** calculate the mode vectors, for regular grid
 * in: x,y: arrays of the grid points, sorted!
 * out:        
 *      M: number of modes
 *      beta: scale factor
 *      Av: array of mode vectors size Nx.Ny times n0.n0, in column major order
 *      n0: number of modes in each dimension
 *
 */
extern int
calculate_mode_vectors(double *x, int Nx, double *y, int Ny, int *M, double
								                *beta, double **Av, int *n0);

/* calculate mode vectors for each (x,y) point given by the arrays x, y
 * of equal length. can be unregular grid
 *
 * in: x,y: arrays of the grid points, unsorted
 *      N: number of grid points
 *      beta: scale factor
 *      n0: number of modes in each dimension
 * out:        
 *      Av: array of mode vectors size N times n0.n0, in column major order
 *
 */
extern int
calculate_mode_vectors_bi(double *x, double *y, int N,  double beta, int n0, double **Av);

/* no sorting etc, just simple calculations */
extern int
calculate_mode_vectors_simple(double *x, double *y, int N,  double beta, int n0, double **Av);

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
 *     Av: size Nx*Ny time n0*n0 matrix in column major order
 *     n0: no. of modes in each dimension
 *     beta: scale
 *  if n0<0 or beta<0, default values will be calculated
 */
extern int
calculate_mode_vectors_tf(double *x, int Nx, double *y, int Ny,
                   double a, double b, double theta,
								        int *n0, double *beta, double **Av);

/*************************************************
 * shapelet_uv.c
 ************************************************/
/** calculate the UV mode vectors
 * in: u,v: arrays of the grid points in UV domain
 *      M: number of modes
 *      beta: scale factor
 *      n0: number of modes in each dimension
 * out:
 *      Av: array of mode vectors size Nu.Nv times n0.n0, in column major order
 *      cplx: array of integers, size n0*n0, if 1 this mode is imaginary, else real
 *
 */
extern int
calculate_uv_mode_vectors(double *u, int Nu, double *v, int Nv, double beta, int n0, double **Av, int **cplx);

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
extern int
calculate_uv_mode_vectors_scalar(double *u, int Nu, double *v, int Nv, double beta, int n0, double **Av, int **cplx);

/* calculate mode vectors for each (u,v) point given by the arrays u, v
 * of equal length. can be unregular grid
 *
 * in: u,v: arrays of the grid points, unsorted
 *      N: number of grid points
 *      beta: scale factor
 *      n0: number of modes in each dimension
 * out:        
 *      Av: array of mode vectors size N times n0.n0, in column major order
 *      cplx: array of integers, size n0*n0, if 1 this mode is imaginary, else real
 *
 */
extern int
calculate_uv_mode_vectors_bi(double *u, double *v, int N,  double beta, int n0, double **Av, int **cplx);



/* calculate mode vectors for the regular grid given by the u,v arrays after 
 * performing the given linear transform. i.e.,
 * |X| =    |a  0| |cos(th)  sin(th)| { |u| - |p| }
 * |Y|      |0  b| |-sin(th) cos(th)| { |v| - |q| }
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
 *     Av: size Nu*Nv time n0*n0 matrix in column major order
 *     n0: no. of modes in each dimension
 *     beta: scale
 *      cplx: array of integers, size n0*n0, if 1 this mode is imaginary, else real
 */
extern int
calculate_uv_mode_vectors_tf(double *u, int Nu, double *v, int Nv,
                   double a, double b, double theta, double p, double q,
								        int n0, double beta, double **Av, int **cplx);



/*******************************************************
 * lapack.c
 *******************************************************/
/* solve the linear least squares problem using LAPACK */
/* min_x |Ax-b|_2 norm */
/* A: N by M, N> M 
 * b: N by 1 vector 
 * x: M by 1 vector */
extern int 
lsq_lapack(double *Av,double *b,double *x, int N, int M);


/* solve the linear least squares problem using FISTA */
/* min_x |Ax-b|_2 norm + mu |x|_1 + lambda |x|_2^2 */
/* A: N by M, N> M
 * b: N by 1 vector
 * x: M by 1 vector
 * mu: L1 penalty, lambda L2 penalty */
extern int
elasticnet_fista(double *Av,double *b,double *x, int N, int M, double lambda, double mu, int maxiter);

/*******************************************************
 * myblas.c
 *******************************************************/
/* y = a.x + y */
extern void
my_daxpy(int N, double *x, double a, double *y);
extern void
my_saxpy(int N, float *x, float a, float *y);
/* y = x */
extern void
my_dcopy(int N, double *x, int Nx, double *y, int Ny); 
extern void
my_scopy(int N, float *x, int Nx, float *y, int Ny);
/* scale */
extern void
my_dscal(int N, double a, double *x);
extern void
my_sscal(int N, float a, float *x);
/* norm || ||_2 */
extern double
my_dnrm2(int N, double *x);
extern float
my_snrm2(int N, float *x);
/* max |x| id, start from 1... */
extern int
my_idamax(int N, double *x, int incx);
/* min |x| id, start from 1... */
extern int
my_idamin(int N, double *x, int incx);

/* BLAS DGEMM C = alpha*op(A)*op(B)+ beta*C */
/* op(A) : M x K, op(B) : K x N, C: MxN */
extern void
my_dgemm(char transa, char transb, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc);

/* A=U S VT, so V needs NOT to be transposed */
extern int
my_dgesvd(char JOBU, char JOBVT, int M, int N, double *A, int LDA, double *S,
   double *U, int LDU, double *VT, int LDVT, double *WORK, int LWORK);

/* BLAS SGEMV  y = alpha*op(A)*x+ beta*y : op 'T' or 'N' */
extern void
my_dgemv(char trans, int M, int N, double alpha, double *A, int lda, double *x, int incx,  double beta, double *y, int incy);
extern void
my_sgemv(char trans, int M, int N, float alpha, float *A, int lda, float *x, int incx,  float beta, float *y, int incy);
/*******************************************************
 * decom_fits.c
 ******************************************************/
/*
 * the bulk of the work is done here;
 * filename: input FITS file
 * cutoff: cutoff [0,1)
 * x,Nx: Nx by 1 grid of x axis
 * y,Ny: Ny by 1 grid of y axis
 * beta: scale
 * M: no of modes
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
extern int 
decompose_fits_file(char* filename, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double p, double q, double clipmin, double clipmax, double **b, double **av, double **z, position *cen, int convolve_psf, char *psf_filename, int use_mask); 

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
extern int 
decompose_fits_file_tf(char* filename, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double a0, double b0, double theta, double p, double q, double clipmin, double clipmax, double **b, double **av, double **z, position *cen, int convolve_psf, char *psf_filename, int use_mask); 



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
extern int 
calcmode_fits_file(char* filename, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double p, double q, double **b, double **Av, position *cen, int convolve_psf, char *psf_filename); 


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
calcmode_fits_file_tf(char* filename, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double a0, double b0, double theta, double p, double q, double **b, double **Av, position *cen, int convolve_psf, char *psf_filename); 


/**************************************************************
 fft.c 
**************************************************************/

/* evaluate 2D fft using fftw
   data: column major data, fortran array
   Nx, Ny: data dimension Nx x Ny (rows, columns)
   Nud, Nvd: FFT dimension Nud x Nvd  Nud> Nx and Nvd>Ny (proper zero padding)
   freal, fimag: copied values, size Nu x Nv (smaller than Nud x Nvd), C arrays
   In a nutshell:
   Nx x Ny image -> zero pad to get Nud x Nvd image -> FFT to get Nud x Nvd values 
   -> extract inner Nu x Nv value -> copy back to output (columns,rows)
   colmajor: if 1, return the result in column major order, else row major order
 */
extern int
eval_fft(double *data, int Nx, int Ny, double *freal, double *fimag, int Nud, int Nvd, int Nu, int Nv, int colmajor);

/** convolves the mode vectors with the given PSF  (using the FFT)
 *      x,y: arrays of the grid points, sorted!
 *      Av: array of mode vectors size Nx.Ny times n0.n0, in column major order
 *      n0: number of modes in each dimension
 *      bmaj (radians),bmin (radians),bpa (degrees): Gaussian PSF parameters 
 */
extern int
convolve_with_psf(double *x, int Nx, double *y, int Ny,
								  double *Av, int n0, double bmaj, double bmin, double bpa);

/** convolves the mode vectors with the given PSF fits array (using the FFT)
 *      x,y: arrays of the grid points, sorted!
 *      M: number of modes
 *      Av: array of mode vectors size Nx.Ny times n0.n0, in column major order
 *      n0: number of modes in each dimension
 *      PSF parameters 
 *      psf: array in column major order, size Npx x Npy
 *      constraint Npx<= Nx and Npy <= Ny
 */
extern int
convolve_with_psf_fits(double *x, int Nx, double *y, int Ny,
                  double *Av, int n0, double *psf, int Npx, int Npy);

/**************************************************************
 multi_fits.c 
**************************************************************/
/* read multiple FITS files in the directory (select suitable files)
 * fitsdir: directory of image files
 * cutoff: cutoff to truncate the image, -1 to ignore this, 1.0 to read whole image
 * myarr: 4D data array of truncated image
 * new_naxis: array of dimensions of each axis
 * lgrid: grid points in l axis
 * mgrid: grid points in m axis
 * ignore_wcs: sometimes too small images screw up WCS, set 1 to use manual grid
 *  cen : position (RA (h,m,s) Dec (h,m,s)
 *  xlow,xhigh,ylow,yhigh: if cutoff==-1, then will use these to read a sub image
 p,q: if not zero, will shift center to these pixels (p,q)
* clipmin,clipmax: if not zero, clip values out this range before decomposition
* Nf : total frequencies (= suitable FITS files)
* freqs: array of frequencies Nfx1
* bmaj,bmin,bpa: PSF info Nfx1
*/
extern int 
read_fits_dir(const char *fitsdir, double cutoff, double**myarr, long int *new_naxis, double **lgrid, double **mgrid, io_buff *fbuff, int ignore_wcs, position *cen, int xlow, int xhigh, int ylow, int yhigh,double p, double q, double clipmin, double clipmax, int *Nf, double **freqs, double **bmaj, double **bmin, double **bpa);

/*
 * similar to decompose_fits(), instead of one FITS file, a directory
 * fitsdir: directory
 * cutoff: cutoff [0,1)
 * x,Nx: Nx by 1 grid of x axis
 * y,Ny: Ny by 1 grid of y axis
 * beta: scale
 * M: no of modes
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
extern int 
decompose_fits_dir(const char *fitsdir, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double p, double q, double clipmin, double clipmax, int *Nf, double **freqs, double **b, double **av, double **z, position *cen, int convolve_psf, char *psf_filename, int use_mask);

/**************************************************************
 apc.c 
**************************************************************/
extern int 
apc_decompose_fits_file(char* filename, double cutoff, int *Nx, int *Ny, double *beta, int *M, int *n0, double **img, double **av, double **z, position *cen); 

#ifdef __cplusplus
     } /* extern "C" */
#endif

#endif /* SHAPELET_H */
