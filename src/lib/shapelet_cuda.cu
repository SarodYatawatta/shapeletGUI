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

#include <cuda_runtime.h>
#include "shapelet.h"

#define DEFAULT_TH_PER_BK 64

/* Hermite polynomial, non recursive version */
__device__ float
H_e(float x, int n) {
  if(n==0) return 1.0f;
  if(n==1) return 2.0f*x;
  /* else iterate */
  float Hn_1,Hn,Hnp1;
  Hn_1=1.0f;
  Hn=2.0f*x;
  int ci;
  for (ci=1; ci<n; ci++) {
    Hnp1=2.0f*x*Hn-2.0f*((float)ci)*Hn_1;
    Hn_1=Hn;
    Hn=Hnp1;
  }

  return Hn;
}

__global__ void 
kernel_calculate_shapelet_lm(float *Ad,float *xd,float *yd,float *fact,float beta,int N,int n0, int startpix, int endpix) {
  /* pixel 0..N, n+startpix */
  unsigned int n = threadIdx.x + blockDim.x*blockIdx.x + startpix;
  /* mode 0...n0^2 */
  unsigned int mode = threadIdx.y + blockDim.y*blockIdx.y;
  /* separate mode to n1,n2 */
  unsigned int n1=mode%n0;
  unsigned int n2=mode/n0;

  if (n<=endpix && n1<n0 && n2<n0) {
   float xx=xd[n]/beta;
   float yy=yd[n]/beta;

   Ad[n+mode*N]=H_e(xx,n1)/sqrtf(powf(2.0f,(float)n1+1)*fact[n1])*expf(-0.5f*xx*xx)
    *H_e(yy,n2)/sqrtf(powf(2.0f,(float)n2+1)*fact[n2])*expf(-0.5f*yy*yy);

  }
}

extern "C" {

static void
checkCudaError(cudaError_t err, const char *file, int line)
{
    if(!err)
        return;
    fprintf(stderr,"GPU (CUDA): %s %s %d\n", cudaGetErrorString(err),file,line);
    exit(EXIT_FAILURE);
}

int 
calculate_mode_vectors_cuda(double *x, double *y, int N,  double beta, int n0, double **Av) {

  cudaError_t error;
  error = cudaGetLastError();

  /* set up factorial array */
  float *fact;
  if ((fact=(float*)calloc((size_t)(n0),sizeof(float)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  fact[0]=1.0f;
  for (int ci=1; ci<(n0); ci++) {
    fact[ci]=float(ci)*fact[ci-1];
  }

  float *Ad=0;
  float *factd=0;
  float *xd=0;
  float *yd=0;
  /* setup device memory */
  error=cudaMalloc((void**)&Ad, (size_t)(N*(n0)*(n0))*sizeof(float));
  checkCudaError(error,__FILE__,__LINE__);
  error=cudaMalloc((void**)&factd, (size_t)(n0)*sizeof(float));
  checkCudaError(error,__FILE__,__LINE__);
  error=cudaMemcpy(factd,fact,(size_t)(n0)*sizeof(float),cudaMemcpyHostToDevice);
  checkCudaError(error,__FILE__,__LINE__);
  error=cudaMalloc((void**)&xd, (size_t)(N)*sizeof(float));
  checkCudaError(error,__FILE__,__LINE__);
  error=cudaMalloc((void**)&yd, (size_t)(N)*sizeof(float));
  checkCudaError(error,__FILE__,__LINE__);
  float *xf,*yf;
  if ((xf=(float*)calloc((size_t)(N),sizeof(float)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((yf=(float*)calloc((size_t)(N),sizeof(float)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  for (int ci=0; ci<N; ci++) {
    xf[ci]=(float)x[ci];
    yf[ci]=(float)y[ci];
  }
  error=cudaMemcpy(xd,xf,(size_t)(N)*sizeof(float),cudaMemcpyHostToDevice);
  checkCudaError(error,__FILE__,__LINE__);
  error=cudaMemcpy(yd,yf,(size_t)(N)*sizeof(float),cudaMemcpyHostToDevice);
  checkCudaError(error,__FILE__,__LINE__);

  free(xf);
  free(yf);

  int ThreadsPerBlock=DEFAULT_TH_PER_BK;
  /* if pixels x modes / threads > 16384, split the calculations into smaller chunks */
  int MAX_BLOCKS=16384;
  if (n0*n0*N/ThreadsPerBlock < MAX_BLOCKS) {
    /* 2D grid of threads: x dim-> pixel (N), y dim-> shapelet mode (n0*n0) */
    dim3 grid(1, 1, 1);
    grid.x = (int)ceilf(N / (float)ThreadsPerBlock);
    grid.y = n0*n0;
    int startpix=0;
    int endpix=N-1;

    kernel_calculate_shapelet_lm<<<grid,ThreadsPerBlock>>>(Ad,xd,yd,factd,(float)beta,N,n0,startpix,endpix);
  } else {
    /* divide N pixels into blocks so that
       pixels * modes / threads < 16384 */
    int n_runs=(int)ceilf(ceilf((N*n0*n0)/(float)ThreadsPerBlock)/(float)MAX_BLOCKS);
    int pix_per_run=N/n_runs;
    for (int run=0; run<n_runs; run++) {
      int startpix=run*pix_per_run;
      int endpix=(run+1)*pix_per_run-1;
      if (endpix>N-1) { endpix=N-1; }
      dim3 grid(1, 1, 1);
      grid.x = (int)ceilf((endpix-startpix+1) / (float)ThreadsPerBlock);
      grid.y = n0*n0;

      kernel_calculate_shapelet_lm<<<grid,ThreadsPerBlock>>>(Ad,xd,yd,factd,(float)beta,N,n0,startpix,endpix);
    }
  }
  cudaDeviceSynchronize();
  error = cudaGetLastError();
  checkCudaError(error,__FILE__,__LINE__);

  free(fact);
  error=cudaFree(factd);
  checkCudaError(error,__FILE__,__LINE__);
  error=cudaFree(xd);
  checkCudaError(error,__FILE__,__LINE__);
  error=cudaFree(yd);
  checkCudaError(error,__FILE__,__LINE__);

  float *Af=0;
  if ((Af=(float*)calloc((size_t)(N*(n0)*(n0)),sizeof(float)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  error=cudaMemcpy(Af,Ad,(size_t)(N*(n0)*(n0))*sizeof(float),cudaMemcpyDeviceToHost);
  checkCudaError(error,__FILE__,__LINE__);
  error=cudaFree(Ad);
  checkCudaError(error,__FILE__,__LINE__);

  if ((*Av=(double*)calloc((size_t)(N*(n0)*(n0)),sizeof(double)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  for (int ci=0; ci<N*n0*n0; ci++) {
    (*Av)[ci]=(double)Af[ci];
  }

  free(Af);
  return 0;

}

}
