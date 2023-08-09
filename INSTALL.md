wo  9 aug 2023 12:32:11 CEST
# ShapeletGUI Installation

## Pre-requisites:
Libraries: cfitsio, wcslib, fftw, openblas and GUI: QT development tools.
In Ubuntu 20.04: libcfitsio-dev, wcslib-dev, libopenblas-dev, libfftw3-dev and QT development packages including qmake.

## Download and unpack the source:
Get the source from [here](https://github.com/SarodYatawatta/shapeletGUI).

## Compile the binary (cmake)
Using cmake it is simple as the following:

```
 mkdir build && cd build
 cmake ..
 make
```

## Compile the library (without cmake)
  * change to src/lib 
  * edit Makefile
  *  change WCSI to point to wcslib include files (wcs.h)
  *  change WCSL to point to wcslib library (libwcs.so)
  * also change INCLUDES and LIBPATH to find cfitsio headers (fitsio2.h) and lib (libcfitsio.so)
  * type "make libshapelet.a"

## Compile the GUI (without cmake)
 * change to src/shapeletGUI
 * edit shapeletGUI.pro file
 * change INCLUDEPATH and LIBS to include wcslib and cfitsio (same as step 2)
 * type "qmake"
 * type "make clean"
 * type "make"



## Start using the program "shapeletGUI"
