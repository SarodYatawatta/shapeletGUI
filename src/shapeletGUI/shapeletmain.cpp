#include "shapeletmain.h"
#include "ui_shapeletmain.h"
#include "optionsdialog.h"
#include "textdialog.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QGraphicsPixmapItem>

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

    }
}

void ShapeletMain::on_actionExit_triggered()
{
    this->close();
}

void ShapeletMain::on_actionRun_single_triggered()
{
  QGraphicsPixmapItem *item=new QGraphicsPixmapItem(QPixmap("/rec/modelFITS.png"));
  item->setScale(1.0);
  //item->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
  ui->graphicsView->scene->addItem(item);
  item->setPos(100,100);
  ui->graphicsView->scene->setSceneRect(0,0,400,400);
  ui->graphicsView->scene->addRect(QRectF(0, 0, 100, 100));
  ui->graphicsView->show();
}

void ShapeletMain::on_actionRun_multifrequency_triggered()
{
   ui->graphicsView->scene->clear();
}

void ShapeletMain::on_actionSettings_triggered()
{
 OptionsDialog *opt=new OptionsDialog(this, 100, 2.0, 0.9, 0.0, 1.0, 1.0, false);

 if (opt->exec()== QDialog::Accepted) {
    ui->graphicsView->setModes(opt->modes());
    ui->graphicsView->setCutoff(opt->cutoff());
    ui->graphicsView->setScale(opt->scale());
    ui->graphicsView->setXscale(opt->xscale());
    ui->graphicsView->setYscale(opt->yscale());
    ui->graphicsView->setRotation(opt->rotation());
    ui->graphicsView->setTf(opt->tf());
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
