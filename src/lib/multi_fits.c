/* FST - a Fast Shapelet Transformer
 *
   Copyright (C) 2006-2022 Sarod Yatawatta <sarod@users.sf.net>  
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


#include <dirent.h>
#include "shapelet.h"

/* for scandir() */
static int
select_file(const struct dirent *dent)
{
  if (!strcmp(dent->d_name,"..")) return 0;
  if (strlen(dent->d_name)>1) {
   /* only select files ending with .fits or .FITS */
   if (strstr(dent->d_name,".fits")) {
    return 1;
   }
   if (strstr(dent->d_name,".FITS")) {
    return 1;
   }
  }

  return 0;
}

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
int read_fits_dir(const char *fitsdir, double cutoff, double**myarr, long int *new_naxis, double **lgrid, double **mgrid, io_buff *fbuff, int ignore_wcs, position *cen, int xlow, int xhigh, int ylow, int yhigh,double p, double q, double clipmin, double clipmax, int *Nf, double **freqs, double **bmaj, double **bmin, double **bpa) {

  long int naxis[4]={0,0,0,0};
  io_buff fbuff0;
  double *deltax, *deltay;
  double *x, *y;
  double *pixval;

 /* for handling files in directory */
 char *fullname=0;
 struct dirent **namelist;

 /* scan given directory */
 *Nf=scandir(fitsdir, &namelist, select_file, alphasort);
 if (*Nf<=0) {
     fprintf(stderr,"%s: %d: invalid directory\n",__FILE__,__LINE__);
     return 1;
 }
 /* allocate memory for freqs and PSF */
 if ((*freqs=(double*)calloc((size_t)(*Nf),sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
  if ((*bmaj=(double*)calloc((size_t)*Nf,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((*bmin=(double*)calloc((size_t)*Nf,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((*bpa=(double*)calloc((size_t)*Nf,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((deltax=(double*)calloc((size_t)*Nf,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((deltay=(double*)calloc((size_t)*Nf,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }

 int totalpix=0;
 /* iterate over FITS files */
 int Nm=0;
 /**************************************************************/
 for (int cnt=0; cnt< *Nf; ++cnt) {
       /* create full path */
       if ((fullname=(char*)calloc((size_t)(strlen(fitsdir)+strlen(namelist[cnt]->d_name)+2),sizeof(char)))==0) {
         fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
         exit(1);
       }
       strcpy(fullname,fitsdir);
       /* append extra '/' to catch errors */
       strcat(fullname,"/");
       strcat(fullname,namelist[cnt]->d_name);
       printf("processing file %s\n",fullname);
       /* xlow=xhigh=ylow=yhigh=0 here */
       read_fits_file(fullname,cutoff, &pixval,naxis, &x, &y, &fbuff0,0, cen,0,0,0,0,p,q, clipmin, clipmax, 0, &Nm, &(*bmaj)[cnt], &(*bmin)[cnt], &(*bpa)[cnt], &deltax[cnt], &deltay[cnt], &(*freqs)[cnt]);
       free(fullname);
      if (!cnt) { /* allocate memory for all files, using first file size */
        totalpix=naxis[0]*naxis[1];
        if ((*myarr=(double*)calloc((size_t)totalpix*(*Nf),sizeof(double)))==0) {
         fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
         exit(1);
       }
      } else {
        /* check all other files have same size */
        if (totalpix!=naxis[0]*naxis[1]){
         fprintf(stderr,"%s: %d: FITS files do not match in size\n",__FILE__,__LINE__);
         exit(1);
        }
      }
      /* copy to common buffer */
      memcpy(&(*myarr)[cnt*totalpix],pixval,sizeof(double)*totalpix);

      free(x);
      free(y);
      free(pixval);
      close_fits_file(fbuff0);
 }
 /* copy axis parameters */
 new_naxis[0]=naxis[0];
 new_naxis[1]=naxis[1];
 new_naxis[2]=naxis[2];
 new_naxis[3]=naxis[3];
 /**************************************************************/

 for (int cnt=0; cnt< *Nf; ++cnt) {
   free(namelist[cnt]);
 }
 free(namelist);
 free(deltax);
 free(deltay);
 return 0;
}


int decompose_fits_dir(const char *fitsdir, double cutoff, double **x, int *Nx, double **y, int *Ny, double *beta, int *M, int *n0, double p, double q, double clipmin, double clipmax, int *Nf, double **freqs, double **b, double **av, double **z, position *cen, int convolve_psf, char *psf_filename, int use_mask) {



  return 0;
}
