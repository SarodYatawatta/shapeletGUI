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
 
   if (this->pix_) { free(this->pix_); }
   if (this->x_) { free(this->x_); }
   if (this->y_) { free(this->y_); }

   read_fits_file(this->fileName().toLocal8Bit().data(),this->cutoff(),&(this->pix_),naxis,&(this->x_),&(this->y_),&filep,ignore_wcs,&(this->cen_),xlow,xhigh,ylow,yhigh,this->xoff(),this->yoff(),this->clipmin(),this->clipmax(),use_mask, &Nm);
   close_fits_file(filep);

   double minval;
   double maxval;
   QImage *qim=createArrayImage(this->pix_,static_cast<int>(naxis[0]),static_cast<int>(naxis[1]),&minval,&maxval,true);
   //scale to match canvas size
   QImage qimc=qim->scaled(CANVAS_WIDTH/2, CANVAS_HEIGHT/2, Qt::IgnoreAspectRatio,Qt::SmoothTransformation);
   delete qim;
   scene->clear();
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

  if (this->pix_) { free(this->pix_); }
  if (this->x_) { free(this->x_); }
  if (this->y_) { free(this->y_); }
  if (this->av_) { free(this->av_); }
  if (this->z_) { free(this->z_); }

  decompose_fits_file(this->fileName().toLocal8Bit().data(),this->cutoff(),&(this->x_),&Nx,&(this->y_),&Ny,&beta,&M,&n0,this->xoff(),this->yoff(),this->clipmin(),this->clipmax(),&(this->pix_),&(this->av_),&(this->z_),&(this->cen_),this->convolve_psf(),nullptr,0);

  this->setScale(beta);
  this->setModes(M);
  return 0;
}
