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
