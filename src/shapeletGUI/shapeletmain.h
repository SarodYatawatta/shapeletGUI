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
