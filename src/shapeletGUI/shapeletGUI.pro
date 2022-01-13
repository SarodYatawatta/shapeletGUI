QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    mygraphicsview.cpp \
    optionsdialog.cpp \
    shapeletmain.cpp \
    textdialog.cpp

HEADERS += \
    ../lib/shapelet.h \
    mygraphicsview.h \
    optionsdialog.h \
    shapeletmain.h \
    textdialog.h

FORMS += \
    optionsdialog.ui \
    shapeletmain.ui \
    textdialog.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    resources.qrc

unix: LIBS += -L$$PWD/../lib/ -L/usr/lib -lshapelet -lm -lwcs -lcfitsio -lopenblas -lpthread -lgfortran  -lfftw3_threads -lfftw3 -lm

INCLUDEPATH += . $$PWD/../ /usr/include/wcslib
DEPENDPATH += $$PWD/../

unix: PRE_TARGETDEPS += $$PWD/../lib/libshapelet.a

INCLUDEPATH += $$PWD/../lib
DEPENDPATH += $$PWD/../lib

unix: PRE_TARGETDEPS += $$PWD/../lib/libshapelet.a
