Directory /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/DRIVE_src.Harry
contains a parallelised version of Drive_God_lin using the gcc and
Intel Fortran IA32 32-bit Openmp compiler options. 
A good Openmp tutorial is at: https://computing.llnl.gov/tutorials/openMP/
To use the Intel compiler on CERN lxplus needs to execute:

source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.csh ia32

For the number of cores to use enter: setenv OMP_NUM_THREADS n
where n is from 1 (use for validation runs) to the number available.
Tests show the optimum is about 12 - probably due to i/o issues.

This version requires the number of turns to be processed to be given in
the DrivingTerms input file. If this is more than actually in the BPM
data then extra values of zero will be filled and used.

Files in this directory are:

README                   these instructions
Drive_God_lin.c          the parallelised C-code with #pragma omp constructs.
sussix4drivexxNoO.f      the parallelised Fortran with !$OMP constructs.
Makefile                 the make file for building the executable (needs cleaning!)
Drive_God_lin            the executable created by running make in this directory
Drive_God_lin.o          the compiled C-code
sussix4drivexxNoO.o      the compiled Fortran
TEST                     a directory containing the TEST data from E.Macleod
The TEST directory contains required input files and also
ALLBPMs_linx.orig and ALLBPMs_liny.orig
to compare with the TEST output from the original non-parallel version.

The Makefile uses the static linking option and does not need the Intel compiler
environment on lxplus but may do on other machines.

To run the test case in single and multicore modes:

cd TEST
setenv OMP_NUM_THREADS 1
time ../Drive_God_lin . > ! testout
which on lxplus436 gives:
43.715u 0.220s 0:44.58 98.5%
then
setenv OMP_NUM_THREADS 12
time ../Drive_God_lin . > ! testout
62.255u 0.456s 0:06.26 1001.5%
so real time from 44 to 6 seconds and effectively about 10 cores being used.
Note that the stdout in the testout file is in a random order in the multicore
case while the ALLBPMs files are sorted into the original order. There is often
a small numerical difference in the 6th or 7th significant figures of the Q2RMS
values in line 2 of the ALLBPMs linx and liny files between the single and 
multicore executions. This is put down to the different order in which the
tune values are summed during parallel execution to which the RMS values are
sensitive and is not considered important.

A detail is that a minor bug whereby the last bpm data was processed twice
wasting a little time but not affecting the results has been fixed at line
581 in Drive_God_lin.c. Another minor problem is that the sussix_v4.inp
file created by Drive_God_lin.c defines NLINE= 0 which, using the Intel
Fortran compiler, turns the following read of L, M, K into a zero-trip
loop and all values after are taken out of sequence. This does not affect
execution so has been left as it is but if NLINE is given a value then
ISME= 1 would create a smear file so this has been changed to ISME= 0 at
line 1337 in Drive_God_lin.c.

A change was made on 29.03.2012 to fix a bug when executing with more than one
thread. The global constant WINDOW values were mistakenly included in the
Drive_God_lin.c threadprivate pragma meaning they were undefined (zeroes) in
secondary threads. This resulted in the noise1, co and co2 values being computed
as zero in all secondary threads.
