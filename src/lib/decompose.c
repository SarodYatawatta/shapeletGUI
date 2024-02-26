/* FST - a Fast Shapelet Transformer
 *
   Copyright (C) 2006 Sarod Yatawatta <sarod@users.sf.net>  
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
#include <unistd.h>

#include "shapelet.h"
//#define DEBUG

void
print_help(void) {
   fprintf(stderr,"Usage:\n");
   fprintf(stderr,"-f infile.fits -m modes -b beta -o outfile.fits -j jobs -s out.txt\n");
   fprintf(stderr,"-f : input FITS image\n");
   fprintf(stderr,"-j : number of subtasks\n");
   fprintf(stderr,"-t : number of threads\n");
   fprintf(stderr,"-m : approximate number of shapelet basis functions\n");
   fprintf(stderr,"-b : shapelet basis scale\n");
   fprintf(stderr,"-o : (if given) write model to output FITS image\n");
   fprintf(stderr,"-s : (if given) write model coefficients to output text file\n");
}


/* for getopt() */
extern char *optarg;
extern int optind, opterr, optopt;


int main(int argc, char* argv[]) {
  int c;
  char *infile=0;
  char *outfile=0;
  char *outmodes=0;
  int M=10;
  int J=10;
  int Nt=4;
  double beta=0.01;

  while ((c=getopt(argc,argv,"b:f:j:m:o:s:t:h"))!=-1) {
    switch(c) {
    case 'f':
      if (optarg) {
        infile=(char*)calloc((size_t)strlen((char*)optarg)+1,sizeof(char));
        if ( infile== 0 ) {
         fprintf(stderr,"%s: %d: no free memory",__FILE__,__LINE__);
         exit(1);
        }
        strcpy(infile,(char*)optarg);
      }
      break;
    case 'o':
      if (optarg) {
        outfile=(char*)calloc((size_t)strlen((char*)optarg)+1,sizeof(char));
        if ( outfile== 0 ) {
         fprintf(stderr,"%s: %d: no free memory",__FILE__,__LINE__);
         exit(1);
        }
        strcpy(outfile,(char*)optarg);
      }
      break;
    case 's':
      if (optarg) {
        outmodes=(char*)calloc((size_t)strlen((char*)optarg)+1,sizeof(char));
        if ( outmodes== 0 ) {
         fprintf(stderr,"%s: %d: no free memory",__FILE__,__LINE__);
         exit(1);
        }
        strcpy(outmodes,(char*)optarg);
      }
      break;
    case 'm':
      if (optarg) {
        M=atoi(optarg);
      }
      break;
    case 'j':
      if (optarg) {
        J=atoi(optarg);
      }
      break;
    case 't':
      if (optarg) {
        Nt=atoi(optarg);
      }
      break;
    case 'b':
      if (optarg) {
        beta=atof(optarg);
      }
      break;
    default:
      print_help();
      break;
    }
  } 

  if (!infile) {
    print_help();
    exit(0);
  }

  int Nx,Ny;
  int n0;
  double *image=0;
  double *coeff=0;
  double *model=0;
  position center;
  apc_decompose_fits_file(infile,1.0,&Nx, &Ny, &beta,&M,&n0,&image,&coeff,&model,&center,outfile,J,Nt);


  if (outmodes) {
   save_decomposition(outmodes,beta,n0,coeff,center);
  }
  if(infile) free(infile);
  if(outfile) free(outfile);
  if(outmodes) free(outmodes);
  if(image) free(image);
  if(coeff) free(coeff);
  if(model) free(model);
  return 0;
}
