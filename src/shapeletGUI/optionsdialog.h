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
                           double yscale=1.0, bool tf=false);
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

private slots:
    void on_lineEdit_modes_textChanged(const QString &arg1);

    void on_lineEdit_cutoff_textChanged(const QString &arg1);

    void on_lineEdit_rotation_textChanged(const QString &arg1);

    void on_lineEdit_xscale_textChanged(const QString &arg1);

    void on_lineEdit_yscale_textChanged(const QString &arg1);

    void on_checkBox_tf_toggled(bool checked);

private:
    Ui::OptionsDialog *ui;
    int modes_;
    double scale_;
    double cutoff_;
    double rotation_;
    double xscale_;
    double yscale_;
    bool tf_;
};

#endif // OPTIONSDIALOG_H
