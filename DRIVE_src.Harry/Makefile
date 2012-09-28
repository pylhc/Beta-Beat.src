
FC = g77
FC = ifort
 
FFLAGS = -c -m32
FFLAGS = -O3 -funroll-loops -c -m32 -g -pg
FFLAGS = -O3 -c -g -pg -m32 -openmp -diag-enable openmp -openmp-report1 -Bstatic -v -static -static-libgcc
FFLAGS = -O3 -c -m32 -openmp -Bstatic -static -static-libgcc

CFLAGS = -ansi -g  -pg -c -m32 -static -static-libgcc
CFLAGS = -ansi -c -m32 -fopenmp -static-libgcc

FCL = -lm -g77libs -m32 -pg
FCL = -lm -m32 -openmp
FCL =  -lm -m32 -nofor_main -openmp -Bstatic -static-libgcc -static


#Drive_God:  Drive_God.o  sussix4drivexxNoO.o   Makefile
#	g77 -o Drive_God Drive_God.o  sussix4drivexxNoO.o   -lm -g77libs


all: Drive_God_lin

sussix4drivexxNoO.o: sussix4drivexxNoO.f Makefile
	$(FC)  $(FFLAGS)  sussix4drivexxNoO.f

#Drive_God.o: Drive_God.c Makefile
#	gcc  -ansi -g  -c Drive_God.c


Drive_God_lin:  Drive_God_lin.o  sussix4drivexxNoO.o   Makefile
	$(FC) -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o   $(FCL)
#	gcc -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o   $(FCL)
#	ifort $(FCL) -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o

Drive_God_lin.o: Drive_God_lin.c Makefile
	gcc  $(CFLAGS)  Drive_God_lin.c



Drive_God_lin_dev:  Drive_God_lin_dev.o  sussix4drivexxNoO.o   Makefile
	$(FC) -o Drive_God_lin_dev Drive_God_lin_dev.o  sussix4drivexxNoO.o   $(FCL)
#	gcc -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o   $(FCL)
#	ifort $(FCL) -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o

Drive_God_lin_dev.o: Drive_God_lin_dev.c Makefile
	gcc  $(CFLAGS)  Drive_God_lin_dev.c






sussix4drivexxNoO.4000.o: sussix4drivexxNoO.4000.f Makefile
	$(FC)  $(FFLAGS)  sussix4drivexxNoO.4000.f

Drive_God_lin_4000.o: Drive_God_lin_4000.c Makefile
	gcc  $(CFLAGS) Drive_God_lin_4000.c

Drive_God_lin_4000:  Drive_God_lin_4000.o  sussix4drivexxNoO.4000.o   Makefile
	$(FC) -o Drive_God_lin_4000 Drive_God_lin_4000.o  sussix4drivexxNoO.4000.o  $(FCL)

clean: 
	rm -f *.o Drive_God_lin
