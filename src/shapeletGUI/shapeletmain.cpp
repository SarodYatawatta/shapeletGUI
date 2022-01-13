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

#include "shapeletmain.h"
#include "ui_shapeletmain.h"
#include "optionsdialog.h"
#include "textdialog.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QGraphicsPixmapItem>
#include <iostream>

ShapeletMain::ShapeletMain(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::ShapeletMain)
{
    ui->setupUi(this);
}

ShapeletMain::~ShapeletMain()
{
    delete ui;
}


void ShapeletMain::on_actionOpen_FITS_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Open FITS"), ".", tr("FITS Files (*.fits *.FITS)"));
    if (fileName==nullptr) {
      QMessageBox msg;
      msg.setText(tr("Cannot open FITS file"));
      msg.exec();
    } else {
std::cout<<"Opening "<<fileName.toLocal8Bit().data()<<std::endl;
    ui->graphicsView->setFileName(fileName);
    ui->graphicsView->readFITSFile();
    }
}

void ShapeletMain::on_actionExit_triggered()
{
    this->close();
}

void ShapeletMain::on_actionRun_single_triggered()
{
  //QGraphicsPixmapItem *item=new QGraphicsPixmapItem(QPixmap("/rec/modelFITS.png"));
  //item->setScale(1.0);
  //item->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
  //ui->graphicsView->scene->addItem(item);
  //item->setPos(100,100);
  //ui->graphicsView->scene->setSceneRect(0,0,400,400);
  //ui->graphicsView->scene->addRect(QRectF(0, 0, 100, 100));
  //ui->graphicsView->show();

  ui->graphicsView->decompose();
}

void ShapeletMain::on_actionRun_multifrequency_triggered()
{
   ui->graphicsView->scene->clear();

std::cout<<"Opening "<< ui->graphicsView->fileName().toLocal8Bit().data()<<std::endl;
}

void ShapeletMain::on_actionSettings_triggered()
{
 OptionsDialog *opt=new OptionsDialog(this,
    ui->graphicsView->modes(),
    ui->graphicsView->scale(),
    ui->graphicsView->cutoff(),
    ui->graphicsView->rotation(),
    ui->graphicsView->xscale(),
    ui->graphicsView->yscale(),
    ui->graphicsView->tf(),
    ui->graphicsView->xoff(),
    ui->graphicsView->yoff()
    );

 if (opt->exec()== QDialog::Accepted) {
    ui->graphicsView->setModes(opt->modes());
    ui->graphicsView->setScale(opt->scale());
    ui->graphicsView->setCutoff(opt->cutoff());
    ui->graphicsView->setRotation(opt->rotation());
    ui->graphicsView->setXscale(opt->xscale());
    ui->graphicsView->setYscale(opt->yscale());
    ui->graphicsView->setTf(opt->tf());
    ui->graphicsView->setXoff(opt->xoff());
    ui->graphicsView->setYoff(opt->yoff());
 }
}

void ShapeletMain::on_actionOpen_Directory_triggered()
{
    QString dirName = QFileDialog::getExistingDirectory(this,
        tr("Open Directory"), ".", QFileDialog::ShowDirsOnly|QFileDialog::DontResolveSymlinks);
    if (dirName==nullptr) {
      QMessageBox msg;
      msg.setText(tr("Cannot open directory"));
      msg.exec();
    } else {
        ui->graphicsView->setDirName(dirName);
    }
}

void ShapeletMain::on_actionAbout_triggered()
{
    TextDialog *abt=new TextDialog(this, QString("<b>ShapeletGUI</b>"));
    abt->exec();
}
