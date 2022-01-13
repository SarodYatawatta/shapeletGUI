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

#ifndef MYGRAPHICSVIEW_H
#define MYGRAPHICSVIEW_H

#include <QWidget>
#include <QGraphicsView>
#include <QGraphicsScene>

#include "shapelet.h"

#ifndef CANVAS_WIDTH
#define CANVAS_WIDTH 600
#endif
#ifndef CANVAS_HEIGHT
#define CANVAS_HEIGHT 600
#endif

class MyGraphicsView : public QGraphicsView
{
public:
 MyGraphicsView(QWidget *parent=0) : QGraphicsView(parent )
 {
   scene=new QGraphicsScene(this);
   // set default scene centered at half width,height with width,height
   scene->setSceneRect(CANVAS_WIDTH/2.0,CANVAS_HEIGHT/2.0,CANVAS_WIDTH,CANVAS_HEIGHT);
   this->setScene(scene);

   // set default options
   this->setModes(-1);
   this->setScale(-1.0);
   this->setCutoff(0.9);
   this->setXscale(1.0);
   this->setYscale(1.0);
   this->setRotation(0.0);
   this->setXoff(0.0);
   this->setYoff(0.0);
   this->setTf(false);
   this->setClipmin(0.0);
   this->setClipmax(0.0);
   this->setConvolve_psf(true);

   this->pix_=nullptr;
   this->x_=nullptr;
   this->y_=nullptr;
   this->av_=nullptr;
   this->z_=nullptr;
 }
 ~MyGraphicsView() {
   delete scene;
   if (this->pix_) { free(this->pix_); }
   if (this->x_) { free(this->x_); }
   if (this->y_) { free(this->y_); }
   if (this->av_) { free(this->av_); }
   if (this->z_) { free(this->z_); }
 }

 int modes() const;
 void setModes(int modes);

 double scale() const;
 void setScale(double scale);

 double cutoff() const;
 void setCutoff(double cutoff);

 double rotation() const;
 void setRotation(double rotation);

 double xscale() const;
 void setXscale(double xscale);

 double yscale() const;
 void setYscale(double yscale);

 bool tf() const;
 void setTf(bool tf);


 QString fileName() const;
 void setFileName(QString file_name);

 QString dirName() const;
 void setDirName(QString dir_name);

 double xoff() const;
 void setXoff(double xoff);

 double yoff() const;
 void setYoff(double yoff);

 QGraphicsScene *scene;

 int readFITSFile();

 int decompose();

 int saveDecomp();

 double clipmin() const;
 void setClipmin(double clipmin);

 double clipmax() const;
 void setClipmax(double clipmax);

 bool convolve_psf() const;
 void setConvolve_psf(bool convolve_psf);

 QString saveName() const;
 void setSaveName(const QString &save_name);

protected:
 QRgb getRGB(double z, double maxval);
 QImage* createArrayImage(double *data, int Nx, int Ny, double *minval, double *maxval, bool fortran_array=false);
 int saveDecomp_ascii(const char* filename, double beta, int n0, double *modes, position& cen);
private:

 // variables for decomposition
 int modes_; // upper limit for no. of modes
 int n0_; // actual number of modes, modes>n0*n0
 double scale_; // scale=beta
 double cutoff_;
 // linear transform parameters
 double rotation_;
 double xscale_;
 double yscale_;
 bool tf_;
 double xoff_;
 double yoff_;
 double clipmin_;
 double clipmax_;
 // convolve with PSF
 int convolve_psf_;

 QString file_name_;
 QString dir_name_;
 QString save_name_;

 // variables for FITS IO
 position cen_;
 // storage (will be allocated while reading file)
 double *pix_; // pixel values
 double *x_; // l grid values
 double *y_; // m grid values
 // array of solutions
 double *av_; // size n0*n0
 // reconstructed image
 double *z_; // size equal to image
};



#endif /* MYGRAPHICSVIEW_H */
