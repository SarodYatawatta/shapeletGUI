/*
 *
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

#include "mygraphicsview.h"
#include <QGraphicsPixmapItem>
#include <QMessageBox>

#include <iostream>
int MyGraphicsView::modes() const
{
return modes_;
}

void MyGraphicsView::setModes(int modes)
{
modes_ = modes;
}

double MyGraphicsView::scale() const
{
return scale_;
}

void MyGraphicsView::setScale(double scale)
{
scale_ = scale;
}

double MyGraphicsView::cutoff() const
{
return cutoff_;
}

void MyGraphicsView::setCutoff(double cutoff)
{
cutoff_ = cutoff;
}

double MyGraphicsView::rotation() const
{
return rotation_;
}

void MyGraphicsView::setRotation(double rotation)
{
rotation_ = rotation;
}

double MyGraphicsView::xscale() const
{
return xscale_;
}

void MyGraphicsView::setXscale(double xscale)
{
xscale_ = xscale;
}

double MyGraphicsView::yscale() const
{
return yscale_;
}

void MyGraphicsView::setYscale(double yscale)
{
yscale_ = yscale;
}

bool MyGraphicsView::tf() const
{
return tf_;
}

void MyGraphicsView::setTf(bool tf)
{
tf_ = tf;
}

QString MyGraphicsView::fileName() const
{
    return file_name_;
}

void MyGraphicsView::setFileName(QString file_name)
{
    file_name_ = file_name;
}

QString MyGraphicsView::dirName() const
{
    return dir_name_;
}

void MyGraphicsView::setDirName(QString dir_name)
{
    dir_name_ = dir_name;
}

double MyGraphicsView::xoff() const
{
    return xoff_;
}

void MyGraphicsView::setXoff(double xoff)
{
    xoff_ = xoff;
}

double MyGraphicsView::yoff() const
{
    return yoff_;
}

void MyGraphicsView::setYoff(double yoff)
{
    yoff_ = yoff;
}

QRgb MyGraphicsView::getRGB(double z, double maxval) {
 if (maxval==0) {
    return qRgb(1,1,1);
 }

 double cl=z/maxval;
 //std::cout<<"z="<<z<<"/ "<<maxval<<"="<<cl<<std::endl;
 if (cl<0.25) {
  return qRgb(0,int(cl*256*4),255);
 } else if (cl<0.5) {
  return qRgb(0,255,int((2-cl*4)*256));
 } else if (cl<0.75) {
  return qRgb(int((cl*4-2)*256),255,0);
 }
 // else
  return qRgb(255,int((4-cl*4)*256),0);
}

QImage *MyGraphicsView::createArrayImage(double *data, int Nx, int Ny, double *minval, double *maxval, bool fortran_array) {
  int idmax=idamax(Nx*Ny,data,1)-1; /* 0.. */
  *maxval=data[idmax];
  double *data1=new double[Nx*Ny];
  dcopy(Nx*Ny,data,data1);
  /* subtract maxval to get the -ve high value */
  daxpy(Nx*Ny,data,-1.0*(*maxval),data1);
  int idmin=idamin(Nx*Ny,data1,1)-1;
  *minval=data[idmin];
  std::cout<<"Max "<<idmax<<" "<<*maxval<< " Min "<<idmin<<" "<<*minval<<std::endl;
  //create Image
  QImage *qim;
  if (fortran_array) {
   qim=new QImage(Nx,Ny,QImage::Format_RGB32);
  } else {
   qim=new QImage(Ny,Nx,QImage::Format_RGB32);
  }
  qim->fill(qRgb(255,255,255));
  //fill pixel values -- reverse x
  std::cout<<"Fort "<<fortran_array<<std::endl;
  if (fortran_array) {
   for (int ii=0; ii<Nx; ii++) {
    for (int jj=0; jj<Ny; jj++) {
      //qim->setPixel(ii,jj,getRGB(ia(Ny-jj-1,ii)-*minval,*maxval-*minval));
      qim->setPixel(ii,jj,getRGB(data[(Ny-1-jj)*Nx+ii]-*minval,*maxval-*minval));
    }
   }
  } else {
   for (int ii=0; ii<Ny; ii++) {
    for (int jj=0; jj<Nx; jj++) {
      qim->setPixel(ii,jj,getRGB(data[ii*Ny+jj]-*minval,*maxval-*minval));
    }
   }
  }

  delete[] data1;
  return qim;

}

QImage *MyGraphicsView::createDiffArrayImage(double *data1, double *data2, int Nx, int Ny, double *minval, double *maxval) {
  // data= data1-data2
  double *data=new double[Nx*Ny];
  dcopy(Nx*Ny,data1,data);
  daxpy(Nx*Ny,data2,-1.0,data);
  int idmax=idamax(Nx*Ny,data,1)-1; /* 0.. */
  *maxval=data[idmax];
  // data=data2-data1
  dcopy(Nx*Ny,data2,data);
  daxpy(Nx*Ny,data1,-1.0,data);
  int idmin=idamax(Nx*Ny,data,1)-1;
  *minval=data[idmin];
  if(*minval > *maxval) { double tmp=*maxval; *maxval=*minval; *minval=tmp; }
  std::cout<<"Max "<<*maxval<< " Min "<<*minval<<std::endl;
  //create Image
  QImage *qim=new QImage(Nx,Ny,QImage::Format_RGB32);
  qim->fill(qRgb(255,255,255));
  //fill pixel values
  for (int ii=0; ii<Nx; ii++) {
    for (int jj=0; jj<Ny; jj++) {
      qim->setPixel(ii,jj,getRGB(data[(Ny-1-jj)*Nx+ii]-*minval,*maxval-*minval));
    }
  }

  delete[] data;
  return qim;

}

QString MyGraphicsView::saveName() const
{
    return save_name_;
}

void MyGraphicsView::setSaveName(const QString &save_name)
{
    save_name_ = save_name;
}

bool MyGraphicsView::convolve_psf() const
{
    return convolve_psf_;
}

void MyGraphicsView::setConvolve_psf(bool convolve_psf)
{
    convolve_psf_ = convolve_psf;
}

double MyGraphicsView::clipmax() const
{
    return clipmax_;
}

void MyGraphicsView::setClipmax(double clipmax)
{
    clipmax_ = clipmax;
}

double MyGraphicsView::clipmin() const
{
    return clipmin_;
}

void MyGraphicsView::setClipmin(double clipmin)
{
    clipmin_ = clipmin;
}


int MyGraphicsView::readFITSFile(void)
{
    
    long int naxis[4]={0,0,0,0};
    io_buff filep;
    int ignore_wcs=0;
   int use_mask=0;
   int Nm=0;
   int xlow=0;
   int xhigh=0;
   int ylow=0;
   int yhigh=0;
 
   double beam_maj, beam_min, beam_pa, freq, deltax, deltay;

   clearMemory();
   scene->clear();
   read_fits_file(this->fileName().toLocal8Bit().data(),this->cutoff(),&(this->pix_),naxis,&(this->x_),&(this->y_),&filep,ignore_wcs,&(this->cen_),xlow,xhigh,ylow,yhigh,this->xoff(),this->yoff(),this->clipmin(),this->clipmax(),use_mask, &Nm, &beam_maj, &beam_min, &beam_pa, &deltax, &deltay, &freq);
   close_fits_file(filep);

   double minval;
   double maxval;
   QImage *qim=createArrayImage(this->pix_,static_cast<int>(naxis[0]),static_cast<int>(naxis[1]),&minval,&maxval,true);
   //scale to match canvas size
   QImage qimc=qim->scaled(CANVAS_WIDTH/2, CANVAS_HEIGHT/2, Qt::IgnoreAspectRatio,Qt::SmoothTransformation);
   delete qim;
   std::cout<<"Scene="<<scene->height()<<","<<scene->width()<<std::endl;
   std::cout<<"Image="<<qimc.height()<<","<<qimc.width()<<std::endl;
   QGraphicsPixmapItem *itm=scene->addPixmap(QPixmap::fromImage(qimc));
   itm->setPos(CANVAS_WIDTH/2, CANVAS_HEIGHT/2);
   itm->setZValue(0.0);
   itm->setToolTip(this->fileName());
   this->show();

   return 0;
}


int MyGraphicsView::decompose(void)
{
  if (this->fileName()==nullptr) {
   QMessageBox msg;
   msg.setText(tr("No input given. Open a FITS file first."));
   msg.exec();
   return 1;
  }
  int Nx,Ny;
  int M=this->modes();
  double beta=this->scale();
  int n0=-1;

  clearMemory();
  scene->clear();
  decompose_fits_file(this->fileName().toLocal8Bit().data(),this->cutoff(),&(this->x_),&Nx,&(this->y_),&Ny,&beta,&M,&n0,this->xoff(),this->yoff(),this->clipmin(),this->clipmax(),&(this->pix_),&(this->av_),&(this->z_),&(this->cen_),this->convolve_psf(),nullptr,0);

  this->setScale(beta);
  this->setModes(M);
  this->n0_=n0;

  double minval;
  double maxval;
  QImage *qim=createArrayImage(this->pix_,Nx,Ny,&minval,&maxval,true);
  //scale to match canvas size
  QImage qimc=qim->scaled(CANVAS_WIDTH/2, CANVAS_HEIGHT/2, Qt::IgnoreAspectRatio,Qt::SmoothTransformation);
  delete qim;
  QGraphicsPixmapItem *itm=scene->addPixmap(QPixmap::fromImage(qimc));
  itm->setPos(CANVAS_WIDTH/2, CANVAS_HEIGHT/2);
  itm->setZValue(0.0);
  itm->setToolTip(this->fileName()+" min "+QString::number(minval)+" max "+QString::number(maxval));

  QImage *qim_mod=createArrayImage(this->z_,Nx,Ny,&minval,&maxval,true);
  //scale to match canvas size
  QImage qimc_mod=qim_mod->scaled(CANVAS_WIDTH/2, CANVAS_HEIGHT/2, Qt::IgnoreAspectRatio,Qt::SmoothTransformation);
  delete qim_mod;
  QGraphicsPixmapItem *itm_mod=scene->addPixmap(QPixmap::fromImage(qimc_mod));
  itm_mod->setPos(CANVAS_WIDTH, CANVAS_HEIGHT/2);
  itm_mod->setZValue(0.0);
  itm_mod->setToolTip(tr("Model min ")+QString::number(minval)+" max "+QString::number(maxval));

  QImage *qim_coef=createArrayImage(this->av_,n0,n0,&minval,&maxval,true);
  //scale to match canvas size
  QImage qimc_coef=qim_coef->scaled(CANVAS_WIDTH/2, CANVAS_HEIGHT/2, Qt::IgnoreAspectRatio,Qt::FastTransformation);
  delete qim_coef;
  QGraphicsPixmapItem *itm_coef=scene->addPixmap(QPixmap::fromImage(qimc_coef));
  itm_coef->setPos(CANVAS_WIDTH/2, CANVAS_HEIGHT);
  itm_coef->setZValue(0.0);
  itm_coef->setToolTip(tr("Coefficients min ")+QString::number(minval)+" max "+QString::number(maxval));

  QImage *qim_res=createDiffArrayImage(this->pix_,this->z_,Nx,Ny,&minval,&maxval);
  //scale to match canvas size
  QImage qimc_res=qim_res->scaled(CANVAS_WIDTH/2, CANVAS_HEIGHT/2, Qt::IgnoreAspectRatio,Qt::FastTransformation);
  delete qim_res;
  QGraphicsPixmapItem *itm_res=scene->addPixmap(QPixmap::fromImage(qimc_res));
  itm_res->setPos(CANVAS_WIDTH, CANVAS_HEIGHT);
  itm_res->setZValue(0.0);
  itm_res->setToolTip(tr("Residual min ")+QString::number(minval)+" max "+QString::number(maxval));


  this->show();
  return 0;
}


int MyGraphicsView::saveDecomp_ascii(const char* filename, double beta, int n0, double *modes, position& cen) {

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
  fprintf(fp,"L %lf %lf %lf\n",this->xscale(),this->yscale(),(this->rotation()+90.0)*M_PI/180.0);
  /* last lines additional info, save info on any linear transform used */
  fprintf(fp,"#\n#\n");
  fprintf(fp,"#a=%lf b=%lf theta=%lf p=%lf q=%lf\n",this->xscale(),this->yscale(),(this->rotation()+90.0)*M_PI/180.0,this->xoff(),this->yoff());
  /* save filename, original beta */
  fprintf(fp,"#file=%s beta=%lf ",filename,beta);
  if (this->convolve_psf_) {
   fprintf(fp," convolved with PSF\n");
  } else {
   fprintf(fp,"\n");
  }
  /* as help save full line to be included in sky model */
  /*  name h m s d m s I Q U V spectral_index RM extent_X(rad) extent_Y(rad) pos_angle(rad) freq0 */
  fprintf(fp,"# LSM format:\n");
  fprintf(fp,"## %s %d %d %lf %d %d %lf 1 0 0 0 0 0 %lf %lf %lf 1000000.0\n",filename,cen.ra_h,cen.ra_m,cen.ra_s, cen.dec_d,cen.dec_m,cen.dec_s,this->xscale(),this->yscale(),(this->rotation())*M_PI/180.0);
  fclose(fp);
  return 0;
}

int MyGraphicsView::getNf() const
{
    return Nf_;
}

void MyGraphicsView::setNf(int Nf)
{
    Nf_ = Nf;
}


int MyGraphicsView::saveDecomp(void)
{
    // check to see if we have a valid result to save
    if (this->av_==nullptr) {
   QMessageBox msg;
   msg.setText(tr("No valid shapelet decomposition to save. First open a FITS file and run decomposition."));
   msg.exec();
   return 1;
  }
  if (this->saveName().endsWith(".modes")) {
   saveDecomp_ascii(this->saveName().toLocal8Bit().data(), this->scale(), this->n0_, this->av_, this->cen_);
  } else {
    QString fullname=this->saveName()+".modes";
    saveDecomp_ascii(fullname.toLocal8Bit().data(), this->scale(), this->n0_, this->av_, this->cen_);
  }
  return 0;
}


int MyGraphicsView::readFITSDir(void)
{
    
   long int naxis[4]={0,0,0,0};
   io_buff filep;
   int ignore_wcs=0;
   int xlow=0;
   int xhigh=0;
   int ylow=0;
   int yhigh=0;
 
   clearMemory();
   scene->clear();
   read_fits_dir(this->dirName().toLocal8Bit().data(),this->cutoff(),&(this->pix_),naxis,&(this->x_),&(this->y_),&filep,ignore_wcs,&(this->cen_),xlow,xhigh,ylow,yhigh,this->xoff(),this->yoff(),this->clipmin(),this->clipmax(),&(this->Nf_),&(this->freqs_), &(this->bmaj_), &(this->bmin_), &(this->bpa_));

   return 0;
}
