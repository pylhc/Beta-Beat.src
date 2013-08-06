#Windows Makefile, to eventually be merged with linuxMakefile
#For GUI, must rename Drive_God_lin.exe to Drive_God_lin in linux

FC = ifort

ifdef SystemRoot
	CC = icl
	RM = del -f *.obj
	FFLAGS = /O3 /Qopenmp /Qparallel /Qunroll /Qipo-c 
	CCFLAGS = /c /Qopenmp /Qparallel /Qunroll /Qipo-c /Wall
	LINK = xilink -qipo /out:Drive_God_lin_win.exe
	FCL = 
	DRIVE_OBJ = Drive_God_lin.obj
	SUSSIX_OBJ = sussix4drivexxNoO.obj
else
	CC = icc
	RM = rm -f *.o
	FFLAGS = -O3 -c -m32 -openmp -Bstatic -static -static-libgcc -parallel -unroll -ipo
	CCFLAGS = -ansi -c -m32 -openmp -static-libgcc -parallel -unroll -ipo -Wall -Wremarks -diag-disable 981
	LINK = ifort -o Drive_God_lin 
	FCL = -lm -m32 -nofor_main -openmp -Bstatic -static-libgcc -static -parallel -unroll -ipo -lstdc++
	DRIVE_OBJ = Drive_God_lin.o
	SUSSIX_OBJ = sussix4drivexxNoO.o
endif	


all: Drive_God_lin

# make sussix object file from fortran source
$(SUSSIX_OBJ):
	$(FC) $(FFLAGS) sussix4drivexxNoO.f 
	
# make drive object file from c source	
$(DRIVE_OBJ):
	$(CC) $(CCFLAGS) Drive_God_lin.cpp 

# link to binary
Drive_God_lin: $(DRIVE_OBJ) $(SUSSIX_OBJ)
	$(LINK) $(DRIVE_OBJ) $(SUSSIX_OBJ) $(FCL)

	
clean:
	$(RM) Drive_God_lin


 
#From old Makefile: Wasn't being used, may be left over from a previous version. 

#Drive_God_lin_dev:  Drive_God_lin_dev.o  sussix4drivexxNoO.o   Makefile
#$(FC) -o Drive_God_lin_dev Drive_God_lin_dev.o  sussix4drivexxNoO.o   $(FCL)

#Drive_God_lin_dev.o: Drive_God_lin_dev.c Makefile
#$(CC)  $(CCFLAGS)  Drive_God_lin_dev.c
