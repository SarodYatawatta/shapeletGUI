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

#ifndef OPTIONSDIALOG_H
#define OPTIONSDIALOG_H

#include <QDialog>

namespace Ui {
class OptionsDialog;
}

class OptionsDialog : public QDialog
{
    Q_OBJECT

public:
    explicit OptionsDialog(QWidget *parent = nullptr, int modes=10, double scale=1.0,
                           double cutoff=0.9, double rotation=0.0, double xscale=1.0,
                           double yscale=1.0, bool tf=false, double xoff=0.0, double yoff=0.0);
    ~OptionsDialog();

    int modes() const;
    void setModes(const int modes);

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

    double xoff() const;
    void setXoff(double xoff);

    double yoff() const;
    void setYoff(double yoff);

private slots:
    void on_lineEdit_modes_textChanged(const QString &arg1);

    void on_lineEdit_cutoff_textChanged(const QString &arg1);

    void on_lineEdit_rotation_textChanged(const QString &arg1);

    void on_lineEdit_xscale_textChanged(const QString &arg1);

    void on_lineEdit_yscale_textChanged(const QString &arg1);

    void on_checkBox_tf_toggled(bool checked);

    void on_lineEdit_xoff_textChanged(const QString &arg1);

    void on_lineEdit_yoff_textChanged(const QString &arg1);

    void on_lineEdit_scale_textChanged(const QString &arg1);

private:
    Ui::OptionsDialog *ui;
    int modes_;
    double scale_;
    double cutoff_;
    double rotation_;
    double xscale_;
    double yscale_;
    bool tf_;
    double xoff_;
    double yoff_;
};

#endif // OPTIONSDIALOG_H
