CC=gcc
CXX=g++
CFLAGS= -O3 -Wall -g -DSTANDALONE_FITS
CLIBS= -lm -lcfitsio -lfftw3_threads -lfftw3
WCSI= -I/usr/include/wcslib
WCSL= -L/usr/lib
WCSLIBS= -lwcs


#LAPACK=-L/usr/lib/atlas/sse2 -llapack -lblas
#LAPACK=-L/usr/local/GotoBLAS2/lib -lgoto2 -lpthread -lgfortran
#LAPACK=-L/usr/local/OpenBLAS/lib -lopenblas -lpthread -lgfortran
LAPACK=-L/usr/lib -lopenblas -lpthread -lgfortran

INCLUDES= -I. -I/usr/include $(WCSI) $(LMINCLUDE) $(L1INCLUDE)
LIBPATH= -L/usr/lib $(WCSL) $(LMLIB) $(L1LIB)

RANLIB=ranlib
OBJECTS=hermite.o lapack.o shapelet_lm.o shapelet_uv.o \
       decom_fits.o fits_utils.o fft.o
default:libshapelet.a
fits_utils.o:fits_utils.c shapelet.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
shapelet_lm.o:shapelet_lm.c shapelet.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
shapelet_uv.o:shapelet_uv.c shapelet.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
decom_fits.o:decom_fits.c shapelet.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
lapack.o:lapack.c shapelet.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
hermite.o:hermite.c shapelet.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
fft.o:fft.c shapelet.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

libshapelet.a:$(OBJECTS)
	ar rv libshapelet.a $(OBJECTS); \
	$(RANLIB) libshapelet.a;
