0)Pre-requisites:
 * libcfitsio: get it from http://heasarc.nasa.gov/fitsio/fitsio.html
 * wcslib: get it from http://www.atnf.csiro.au/people/mcalabre/WCS/
 * LAPACK
 * QT, Qmake
Build both and install.
1)Unpack
   gzip -dc fst-*.*.*.tar.gz | tar xvf -
2)Compile the library
  * change to src/lib 
  * edit Makefile
    change WCSI to point to wcslib include files (wcs.h)
    change WCSL to point to wcslib library (libwcs.so)
    also change INCLUDES and LIBPATH to find cfitsio headers (fitsio2.h) and lib (libcfitsio.so)
  * type "make libshapelet.a"

3)Compile the GUI
 * change to src/gui
 * edit gui.pro file
   change INCLUDEPATH and LIBS to include wcslib and cfitsio (same as step 2)
 * type "qmake"
   type "make clean"
   type "make"



4)Start using the program "gui"
