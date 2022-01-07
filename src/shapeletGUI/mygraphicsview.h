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
class MyGraphicsView : public QGraphicsView
{
public:
 MyGraphicsView(QWidget *parent=0) : QGraphicsView(parent )
 {
   scene=new QGraphicsScene(this);
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
 }
 ~MyGraphicsView() {delete scene;}

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

 QGraphicsScene *scene;
 double xoff() const;
 void setXoff(double xoff);

 double yoff() const;
 void setYoff(double yoff);

private:

 int modes_;
 double scale_;
 double cutoff_;
 double rotation_;
 double xscale_;
 double yscale_;
 bool tf_;
 double xoff_;
 double yoff_;

 QString file_name_;
 QString dir_name_;
};



#endif /* MYGRAPHICSVIEW_H */
