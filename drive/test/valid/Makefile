FC = ifort

CC = icc

FFLAGS = -c -m32
FFLAGS = -O3 -funroll-loops -c -m32 -g -pg
FFLAGS = -O3 -c -g -pg -m32 -openmp -diag-enable openmp -openmp-report1 -Bstatic -v -static -static-libgcc
FFLAGS = -O3 -c -m32 -openmp -Bstatic -static -static-libgcc -parallel -unroll -ipo
# -prof-use

CCFLAGS = -ansi -g  -pg -c -m32 -static -static-libgcc
CCFLAGS = -ansi -c -m32 -openmp -static-libgcc -parallel -unroll -ipo -Wall -Wremarks -diag-disable 981
#,1419,2259

FCL = -lm -g77libs -m32 -pg
FCL = -lm -m32 -openmp
FCL = -lm -m32 -nofor_main -openmp -Bstatic -static-libgcc -static -parallel -unroll -ipo
# -prof-use


all: Drive_God_lin

# make sussix object file from fortran source
sussix4drivexxNoO.o: sussix4drivexxNoO.f Makefile
	$(FC)  $(FFLAGS)  sussix4drivexxNoO.f

# make drive object file from c source
Drive_God_lin.o: Drive_God_lin.c Makefile
	$(CC)  $(CCFLAGS)  Drive_God_lin.c

# link to binary
Drive_God_lin:  Drive_God_lin.o  sussix4drivexxNoO.o   Makefile
	$(FC) -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o   $(FCL)


Drive_God_lin_dev:  Drive_God_lin_dev.o  sussix4drivexxNoO.o   Makefile
	$(FC) -o Drive_God_lin_dev Drive_God_lin_dev.o  sussix4drivexxNoO.o   $(FCL)

Drive_God_lin_dev.o: Drive_God_lin_dev.c Makefile
	$(CC)  $(CCFLAGS)  Drive_God_lin_dev.c


clean:
	rm -f *.o Drive_God_lin

