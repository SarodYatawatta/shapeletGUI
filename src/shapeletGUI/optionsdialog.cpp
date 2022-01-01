#include "optionsdialog.h"
#include "ui_optionsdialog.h"

OptionsDialog::OptionsDialog(QWidget *parent, int modes, double scale, double cutoff, double rotation, double xscale, double yscale, bool tf) :
    QDialog(parent),
    ui(new Ui::OptionsDialog),
    modes_(modes), scale_(scale), cutoff_(cutoff),
    rotation_(rotation), xscale_(xscale), yscale_(yscale), tf_(tf)
{
    ui->setupUi(this);
    QLocale::setDefault(QLocale(QLocale::English, QLocale::UnitedStates));

    ui->lineEdit_modes->setText(QString::number(modes_));
    ui->lineEdit_cutoff->setText(QString::number(cutoff_));
    ui->lineEdit_beta->setText(QString::number(scale_));
    ui->lineEdit_rotation->setText(QString::number(rotation_));
    ui->lineEdit_xscale->setText(QString::number(xscale_));
    ui->lineEdit_yscale->setText(QString::number(yscale_));
    ui->checkBox_tf->setChecked(tf_);
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

void OptionsDialog::on_lineEdit_modes_textChanged(const QString &arg1)
{
 bool ok;
 int modes=arg1.toInt(&ok);
 if (ok && modes>0) setModes(modes);
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

void OptionsDialog::on_checkBox_tf_toggled(bool checked)
{
    setTf(checked);
}

bool OptionsDialog::tf() const
{
    return tf_;
}

void OptionsDialog::setTf(bool tf)
{
    tf_ = tf;
}
