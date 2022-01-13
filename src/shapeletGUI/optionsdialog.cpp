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

#include "optionsdialog.h"
#include "ui_optionsdialog.h"

OptionsDialog::OptionsDialog(QWidget *parent, int modes, double scale, 
    double cutoff, double rotation,
 double xscale, double yscale, bool tf, double xoff, double yoff) :
    QDialog(parent),
    ui(new Ui::OptionsDialog),
    modes_(modes), scale_(scale), cutoff_(cutoff),
    rotation_(rotation), xscale_(xscale), yscale_(yscale), tf_(tf),
    xoff_(xoff), yoff_(yoff)
{
    ui->setupUi(this);
    QLocale::setDefault(QLocale(QLocale::English, QLocale::UnitedStates));

    ui->lineEdit_modes->setText(QString::number(modes_));
    ui->lineEdit_cutoff->setText(QString::number(cutoff_));
    ui->lineEdit_scale->setText(QString::number(scale_));
    ui->lineEdit_rotation->setText(QString::number(rotation_));
    ui->lineEdit_xscale->setText(QString::number(xscale_));
    ui->lineEdit_yscale->setText(QString::number(yscale_));
    ui->checkBox_tf->setChecked(tf_);
    ui->lineEdit_xoff->setText(QString::number(xoff_));
    ui->lineEdit_yoff->setText(QString::number(yoff_));
}

OptionsDialog::~OptionsDialog()
{
    delete ui;
}

double OptionsDialog::cutoff() const
{
    return cutoff_;
}

void OptionsDialog::setCutoff(double cutoff)
{
    cutoff_ = cutoff;
}

double OptionsDialog::scale() const
{
    return scale_;
}

void OptionsDialog::setScale(double scale)
{
    scale_ = scale;
}

int OptionsDialog::modes() const
{
    return modes_;
}

void OptionsDialog::setModes(int modes)
{
    modes_ = modes;
}

double OptionsDialog::yscale() const
{
    return yscale_;
}

void OptionsDialog::setYscale(double yscale)
{
    yscale_ = yscale;
}

double OptionsDialog::xscale() const
{
    return xscale_;
}

void OptionsDialog::setXscale(double xscale)
{
    xscale_ = xscale;
}

double OptionsDialog::rotation() const
{
    return rotation_;
}

void OptionsDialog::setRotation(double rotation)
{
    rotation_ = rotation;
}

double OptionsDialog::yoff() const
{
    return yoff_;
}

void OptionsDialog::setYoff(double yoff)
{
    yoff_ = yoff;
}

double OptionsDialog::xoff() const
{
    return xoff_;
}

void OptionsDialog::setXoff(double xoff)
{
    xoff_ = xoff;
}

bool OptionsDialog::tf() const
{
    return tf_;
}

void OptionsDialog::setTf(bool tf)
{
    tf_ = tf;
}

void OptionsDialog::on_lineEdit_modes_textChanged(const QString &arg1)
{
 bool ok;
 int modes=arg1.toInt(&ok);
 if (ok && modes<= 0) modes=-1;
 if (ok) setModes(modes);
}

void OptionsDialog::on_lineEdit_cutoff_textChanged(const QString &arg1)
{
    bool ok;
    double tempval=arg1.toDouble(&ok);
    if (ok && tempval>0.0) setCutoff(tempval);
}

void OptionsDialog::on_lineEdit_rotation_textChanged(const QString &arg1)
{
    bool ok;
    double tempval=arg1.toDouble(&ok);
    if (ok) setRotation(tempval);
}

void OptionsDialog::on_lineEdit_xscale_textChanged(const QString &arg1)
{
    bool ok;
    double tempval=arg1.toDouble(&ok);
    if (ok && tempval>0.0) setXscale(tempval);
}

void OptionsDialog::on_lineEdit_yscale_textChanged(const QString &arg1)
{
    bool ok;
    double tempval=arg1.toDouble(&ok);
    if (ok && tempval>0.0) setYscale(tempval);
}

void OptionsDialog::on_checkBox_tf_toggled(bool checked)
{
    setTf(checked);
}

void OptionsDialog::on_lineEdit_xoff_textChanged(const QString &arg1)
{
    bool ok;
    double tempval=arg1.toDouble(&ok);
    if (ok) setXoff(tempval);
}

void OptionsDialog::on_lineEdit_yoff_textChanged(const QString &arg1)
{
    bool ok;
    double tempval=arg1.toDouble(&ok);
    if (ok) setYoff(tempval);
}

void OptionsDialog::on_lineEdit_scale_textChanged(const QString &arg1)
{
    bool ok;
    double tempval=arg1.toDouble(&ok);
    if (ok) setScale(tempval);
}
