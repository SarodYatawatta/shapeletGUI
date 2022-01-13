do 13 jan 2022 13:34:00 CET
# ShapeletGUI Installation

## Pre-requisites:
Libraries: cfitsio, wcslib, fftw, openblas and GUI: QT development tools.
In Ubuntu 20.04: libcfitsio-dev, wcslib-dev, libopenblas-dev, libfftw3-dev and QT development packages including qmake.

## Download and unpack the source:
Get the source from [here](https://github.com/SarodYatawatta/shapeletGUI).

## Compile the library
  * change to src/lib 
  * edit Makefile
  *  change WCSI to point to wcslib include files (wcs.h)
  *  change WCSL to point to wcslib library (libwcs.so)
  * also change INCLUDES and LIBPATH to find cfitsio headers (fitsio2.h) and lib (libcfitsio.so)
  * type "make libshapelet.a"

## Compile the GUI
 * change to src/gui
 * edit shapeletGUI.pro file
 * change INCLUDEPATH and LIBS to include wcslib and cfitsio (same as step 2)
 * type "qmake"
 * type "make clean"
 * type "make"



## Start using the program "shapeletGUI"
