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

#ifndef SHAPELETMAIN_H
#define SHAPELETMAIN_H

#include <QMainWindow>
#include <QString>

QT_BEGIN_NAMESPACE
namespace Ui { class ShapeletMain; }
QT_END_NAMESPACE

class ShapeletMain : public QMainWindow
{
    Q_OBJECT

public:
    ShapeletMain(QWidget *parent = nullptr);
    ~ShapeletMain();

private slots:
    void on_actionOpen_FITS_triggered();

    void on_actionExit_triggered();

    void on_actionRun_single_triggered();

    void on_actionRun_multifrequency_triggered();

    void on_actionSettings_triggered();

    void on_actionOpen_Directory_triggered();

    void on_actionAbout_triggered();

private:
    Ui::ShapeletMain *ui;
};
#endif // SHAPELETMAIN_H
