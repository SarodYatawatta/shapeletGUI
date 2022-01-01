#include "shapeletmain.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    ShapeletMain w;
    w.show();
    return a.exec();
}
