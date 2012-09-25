
FC = g77
FC = ifort
 
FFLAGS = -c -m32
FFLAGS = -O3 -funroll-loops -c -m32
FFLAGS = -O3 -c -Bstatic -v -static -static-libgcc

CFLAGS = -ansi -g  -c -m32 -static-libgcc
CFLAGS = -ansi -g  -c -m32 -static -static-libgcc

FCL = -lm -m32
FCL = -lm -g77libs -m32
FCL =  -lm -m32 -nofor_main -Bstatic -v -static-libgcc -static 


#Drive_God:  Drive_God.o  sussix4drivexxNoO.o   Makefile
#	g77 -o Drive_God Drive_God.o  sussix4drivexxNoO.o   -lm -g77libs


all: Drive_God_lin

sussix4drivexxNoO.o: sussix4drivexxNoO.f Makefile
	$(FC)  $(FFLAGS)  sussix4drivexxNoO.f

#Drive_God.o: Drive_God.c Makefile
#	gcc  -ansi -g  -c Drive_God.c


Drive_God_lin:  Drive_God_lin.o  sussix4drivexxNoO.o   Makefile
#	$(FC) -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o   $(FCL)
#	gcc -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o   $(FCL)
	ifort $(FCL) -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o

Drive_God_lin.o: Drive_God_lin.c Makefile
	gcc  $(CFLAGS)  Drive_God_lin.c


sussix4drivexxNoO.4000.o: sussix4drivexxNoO.4000.f Makefile
	$(FC)  $(FFLAGS)  sussix4drivexxNoO.4000.f

Drive_God_lin_4000.o: Drive_God_lin_4000.c Makefile
	gcc  $(CFLAGS) Drive_God_lin_4000.c

Drive_God_lin_4000:  Drive_God_lin_4000.o  sussix4drivexxNoO.4000.o   Makefile
	$(FC) -o Drive_God_lin_4000 Drive_God_lin_4000.o  sussix4drivexxNoO.4000.o  $(FCL)

clean: 
	rm -f *.o Drive_God_lin
