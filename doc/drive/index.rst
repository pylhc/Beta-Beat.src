Drive documentation
===================
The executable *Drive.God_lin* computes a refined Fourier transform and writes frequencies, phases
and amplitudes for all BPMs to an output file. The output files, \*.linx and \*.liny, are in the
Table File System (TFS) format.
This is a parallelized version which uses OpenMP. `Here <https://computing.llnl.gov/tutorials/openMP/>`_ is good tutorial.

*Drive_God_lin.cpp* is a wrapper, written in C++, for *SUSSIX(sussix4drivexxNoO.f)*, written in FORTRAN.

We use make and the Intel compilers icc(Linux)/icl(Windows) and ifort to build the binaries.


Drive folder
------------
The folder Beta-Beat.src/drive contains the following items:
 - **test**: a directory containing a Python Unit test
 - **__init__.py**: Indicates that this folder is a Python package (needed for the subpackage *test*)
 - **Drive_God_lin**: The executable created by running make in this directory on a Linux machine
 - **Drive_God_lin_win.exe**: The executable for windows. Created by running make on a Windows machine
 - **Drive_God_lin.cpp**: The parallelised C++-code with #pragma omp constructs.
 - **Drive_God_lin.o**: The compiled C++-code on Linux machines
 - **Drive_God_lin.obj**: The compiled C++-code on Windows machines
 - **Makefile**: *build*-file which holds instructions to build drive
 - **README**: Text file containing a link to this website
 - **sussix4drivexxNoO.f**: The parallelised Fortran with !$OMP constructs
 - **sussix4drivexxNoO.o**: The compiled Fortran for Linux machines
 - **sussix4drivexxNoO.obj**: The compiled Fortran for Windows machines

Building drive
--------------
To build drive one needs the Intel compilers for C++ and Fortran and the utility *make* which is on
almost every UNIX system available. For windows `MinGW <http://www.mingw.org/>`_ provides *make*.

The Intel compilers are commercial but CERN has licences. The section below shows how to obtain it.

To build drive execute ``make``. In case you have ancient build files use ``make clean`` before.
*make* will build according to your Os either *Drive_God_lin*(Linux) or *Drive_God_lin_win.exe*(Windows).

Linux - Intel compilers
```````````````````````
To use the Intel compiler on CERN lxplus one needs to execute:
::

	``source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.csh ia32``

``icc --version`` and ``ifort --version`` should print version information after that.

Hint: One have to run the csh console instead of bash console otherwise
there will be a syntax error. To switch to the csh console enter ``csh``.

For the number of cores to use enter before running drive ``setenv OMP_NUM_THREADS n`` where ``n`` is from 1 (use for
validation runs) to the number available.
Tests show the optimum is about 12 - probably due to i/o issues.


Windows - Intel compilers
`````````````````````````
The following tutorial shows how to obtain the Intel compilers on windows:
https://twiki.cern.ch/twiki/bin/view/Openlab/IntelTools#Obtaining_and_installing_Int_AN1

``icl /help`` and ``ifort /help`` should print the manual.

Probably one has to add paths to missing header files to the *Path* system variable.
If the compiler misses header files search for them on your system. The header files are probably
in a similar directory to ``C:\Program Files (x86)\Intel\Composer XE 2013\bin\ia32`` or
``C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin``.

Running drive
-------------
To run drive you have to state the path to the output directory as first argument:

::

	./Drive_God_lin path/to/output/

The output directory have to contain *Drive.inp* and *DrivingTerms*.

There is often a small numerical difference in the 6th or 7th significant figures of the Q2RMS
values in line 2 of the linx and liny files between the single and
multicore executions. This is put down to the different order in which the
tune values are summed during parallel execution to which the RMS values are
sensitive and is not considered important.

A minor problem is that the sussix_v4.inp file created by Drive_God_lin.cpp defines NLINE= 0 which,
using the Intel Fortran compiler, turns the following read of L, M, K into a zero-trip
loop and all values after are taken out of sequence. This does not affect
execution so has been left as it is but if NLINE is given a value then
ISME= 1 would create a smear file so this has been changed to ISME= 0 at
line 1337 in Drive_God_lin.c.

Drive.inp
`````````
Drive.inp is a text file which contains settings for Drive.
Here is an example:

::

  KICK=1
  CASE(1[H], 0[V])=1
  KPER(KICK PERCE.)=0.5
  TUNE X=0.31
  TUNE Y=0.32
  PICKUP START=0
  PICKUP END=538
  ISTUN=0.01
  LABEL RUN (1[yes])=0
  WINDOWa1=0.04
  WINDOWa2=0.1
  WINDOWb1=0.4
  NATURAL X=0.30
  NATURAL Y=0.33

DrivingTerms
````````````
The text file DrivingTerms is a settings file for Drive and contains a single line includig the path to the source sdds file, the start turn and the end turn.

**(path to datafile) (int start turn) (int end turn)**
::

	Beam1@Turn@2012_07_09@15_47_56_059_0.sdds.new.new 1 2100

PyUnit test
-----------
test_output.py runs Drive_God_lin for every input folder in *drive/test/data/input/* and outs the output into
a subfolder of *drive/test/data/to_check/*. Afterwards ndiff.py will be used to compare the produced output
with precalculated, assumed files in *drive/test/data/valid*. If one needs to create new valid files the
executable in *drive/test/valid/* can be used. *drive/test/valid/Drive_God_lin* is a previous version of drive.
