#Drive_God:  Drive_God.o  sussix4drivexxNoO.o   Makefile
#	g77 -o Drive_God Drive_God.o  sussix4drivexxNoO.o   -lm -g77libs


all: Drive_God_lin

sussix4drivexxNoO.o: sussix4drivexxNoO.f Makefile
	g77  -c  sussix4drivexxNoO.f

#Drive_God.o: Drive_God.c Makefile
#	gcc  -ansi -g  -c Drive_God.c


Drive_God_lin:  Drive_God_lin.o  sussix4drivexxNoO.o   Makefile
	g77 -o Drive_God_lin Drive_God_lin.o  sussix4drivexxNoO.o   -lm -g77libs

Drive_God_lin.o: Drive_God_lin.c Makefile
	gcc  -ansi -g  -c Drive_God_lin.c


sussix4drivexxNoO_Massimo.o: sussix4drivexxNoO_Massimo.f Makefile
	g77  -c  sussix4drivexxNoO_Massimo.f

Drive_God_lin_4000.o: Drive_God_lin_4000.c Makefile
	gcc  -ansi -g  -c Drive_God_lin_4000.c

Drive_God_lin_4000:  Drive_God_lin_4000.o  sussix4drivexxNoO_Massimo.o   Makefile
	g77 -o Drive_God_lin_4000 Drive_God_lin_4000.o  sussix4drivexxNoO_Massimo.o   -lm -g77libs
