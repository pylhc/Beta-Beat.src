/* Parallelised version of Drive_God_lin.c of 15/05/2011 using Openmp.
   NB this version uses as the number of turns the minimum of the number
   actually read or the number given in the DrivingTerms file.


   Change 05/08/2013:-Removed any rejections in BPMstatus function
   -Changed formatLinFile to be OS compatible, makefile also modified, drive can now
   be compiled using 'make' on both linux and windows
   -Removed some unused and unnecessary code
   -Forced input tune to be in [-0.5,0.5] by taking closes member in this interval
    -Updated code to use c++ string manipulations     -asherman
    12/08: Added two classes and an array of structs for better data organization

   Version <x3> 20121106 (tbach)
   - removed all (?) unused code
   - changed options reading to allow trailing spaces
   - introduced more options and more columns in output for natural tune
   - cleaned source code

   Version <x2> 20121031 (rwasef/tbach)
   - introduced more columns in output, increased variable size for this (rwasef)
   - removed iformat variable and not used code
   - rewritten option file parsing, order is not important anymore
   - removed all unused variables

   Version <x> 20121001 (tbach):
   - fixed comment line reading (by rtomas)
   - changed variable names to more meaningful and readable names
   - fixed errors from static code analysis
   - removed outcommented code
   - removed unused variables and functions
   - formatted the source
   - changed error messages to more helpful content for the user

   Change 10/08/2012 increase size of character strings
   dataFilePath[500], noisefile[500] from 300 to 500. Long datafile name was causing
   a segv when it was being read.
   Increased size of maxturns define to handle more turns.
   Fortran sourcecode changed, too.

   Change 29/03/2012 at line 43: remove window variables from the
   threadprivate pragma. These are global constants which were undefined
   other than in the primary thread so noise1, co and co2 were being
   calculated as zero in all secondary threads.

   Change 29/09/2011 at lines 715 and 724 to find lines with any
   bpm name to sort in order by looking for a " rather than a name string.
   Has matching sussix4drivexxNoO.f      H.Renshall & E.Maclean
   */

 /*Determines if OS is windows or not and declares external fortran function.
Linux and windows have different requirements for names of subroutines in fortran called from C.
Code works on assumption that other operating systems will work with the linux version, but if this
is not the case then it is easy to add another OS to be used*/
 #ifdef _WIN32
extern "C" { void SUSSIX4DRIVENOISE (double *, double*, double*, double*, double*, double*, double*, double*, char*); }
#define OS "Windows32"
#else
extern "C" { void sussix4drivenoise_(double *, double*, double*, double*, double*, double*, double*, double*, char*); }
#define OS "linux"
#endif

//C++ libraries
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>

 //C libraries
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <ctype.h>
#define MAXPICK 1100
#define MAXTURNS 10000 /* Critical parameter in sussixfordiveNoO.f, same parameter name !*/
#define MAXTURNS4 40000 /*Always four times  MAXTURNS*/
#define MAXRUNS 100
#define NATTUNE_DEFAULT -100

#if !defined(LOG_INFO)
#define LOG_INFO 0 /*set to 1 to enable log prints (tbach) */
#endif

void get_and_check_file_paths(std::string *, std::string *, std::string *, std::string *,const char *);
void get_inp_data(std::string);
void check_inp_data();
int readDrivingTerms(std::istream&, int*, std::string *);
void writeSussixInput(std::string, const int, const double, const double, const double);
bool BPMstatus(const int, const int);
void formatLinFile(std::string, const int, const double, const double, const int, const double, const double,int);
bool cannotOpenFile(std::string,char);
bool inputFileCheck(std::string);
bool outputFileCheck(std::string);

class Input_Data{
public:
    double  istun, kper, tunex, tuney, windowa1, windowa2, windowb1, windowb2, nattunex, nattuney;
    int     turns, pickstart, pickend, kcase, kick, labelrun;
    void init_input_values(){ //initializes all input data to 0 except for natural tunes which have default value defined
        istun = kper = tunex = tuney = windowa1 = windowa2 = windowb1 = windowb2 = 0.0;
        turns = pickstart = pickend = kcase = kick = labelrun = 0;
        nattunex = nattuney = NATTUNE_DEFAULT;
    }
    void check_tune(double tune, char plane){ //checks that magnitude of tune is less than 0.5, changes it by integer value until it is
        int change;
        if(tune > 0.5 || tune < -0.5){
            change=0;
            while(tune > 0.5){
                tune -= 1;
                change--;
            }
            while(tune < -0.5){
                tune += 1;
                change++;
            }
            std::cout << "tune_" << plane <<" input increased by " << change << ". Should be less than or equal to 0.5 in magnitude. New value is: " << tune << std::endl;
        }
    }
}InpData;

class Calculated_Data{
public:
    double tunesumx, tunesumy, tune2sumx, tune2sumy, nattunexsum, nattunex2sum, nattuneysum, nattuney2sum;
    int    tunecountx, tunecounty, nattunexcount, nattuneycount;
    void init_calc_values(){ //initializes all values to be calculated to zero
        tunesumx = tunesumy = tune2sumx = tune2sumy = nattunexsum = nattunex2sum = nattuneysum = nattuney2sum = 0.0;
        tunecountx = tunecounty = nattunexcount = nattuneycount = 0;
    }
}CalcData;

    double calculatednattuney, calculatednattunex, calculatednatampy, calculatednatampx, co, co2, maxamp,
           noise1, maxfreq, maxmin, maxpeak, noiseAve;

    double allampsx[300], allampsy[300], allfreqsx[300], allfreqsy[300], amplitude[19],
               doubleToSend[MAXTURNS4 + 4], phase[19], tune[2];

    int nslines=0, Nturns=0;

    struct BPM{ /*Structure for each BPM- has name, position, plane, if it's pickedup, and it's tbtdata*/
    std::string bpmname;
    double bpmpos;
    int plane; /* 0 is horizontal, 1 is vertical*/
    bool pickedup;
    double tbtdata[MAXTURNS];
    }BPMs[MAXPICK];

#pragma omp threadprivate(amplitude, doubleToSend, tune, phase,\
        noise1, noiseAve, maxpeak, maxfreq, maxmin, co, co2,\
        allfreqsx, allampsx, allfreqsy, allampsy)

int main(int argc, char **argv)
{

    char string1[1000];

    int counth, countv, maxcounthv, start, Nbpms,i, j, bpmCounter, columnCounter, horizontalBpmCounter, verticalBpmCounter, flag, loopcounter=0;

    std::ofstream linxFile, linyFile, noiseFile, spectrumFile;
    std::ifstream driveInputFile, drivingTermsFile, dataFile;
    std::string temp_str, dataFilePath, bpmFileName, workingDirectoryPath, sussixInputFilePath,
    driveInputFilePath, drivingTermsFilePath, noiseFilePath, linxFilePath, linyFilePath, spectrumFilePath;

    #ifdef _WIN32 /*Changes minor formatting difference in windows regarding the output of a number in scientific notation.*/
        _set_output_format(_TWO_DIGIT_EXPONENT);
    #endif

    omp_set_dynamic(0);

    //To output scientific notation
    std::cout << std::setiosflags (std::ios::scientific);

    get_and_check_file_paths(&workingDirectoryPath, &drivingTermsFilePath, &driveInputFilePath, &sussixInputFilePath, argv[1]);
    std::cout << "READ THIS:" << workingDirectoryPath << std::endl;

    /* set all options to defaults, it could happen that they are not included in the file (tbach) */
    InpData.init_input_values();

    //Reads input data from Drive.inp into global InpData class
    get_inp_data(driveInputFilePath);

    //Checks tune_x/y, kcase, kick, labelrun, pickend, and pickend
    check_inp_data();


    drivingTermsFile.open(drivingTermsFilePath.c_str());
    if(cannotOpenFile(drivingTermsFilePath,'i')){
        std::cout << "Leaving drive due to error" << std::endl;
        exit(EXIT_FAILURE);
    }

    /* From drivingTermsFilePath assign dataFilePath, assign turns. */
    while (readDrivingTerms(drivingTermsFile, &(InpData.turns), &dataFilePath)) {
        /* set all values to be calculated to default values */
        CalcData.init_calc_values();

        dataFile.open(dataFilePath.c_str());
        /* Check the file dataFilePath */
        if(cannotOpenFile(dataFilePath,'i'))
            continue;
        loopcounter++;
        std::cout << "Data file: " << dataFilePath << std::endl;

        bpmFileName = dataFilePath.substr(dataFilePath.rfind('/',dataFilePath.size()-2));

        std::cout << "bpmFileName: " << bpmFileName << std::endl;

        noiseFilePath = workingDirectoryPath+'/'+bpmFileName+"_noise";
        linxFilePath = workingDirectoryPath+'/'+bpmFileName+"_linx";
        linyFilePath = workingDirectoryPath+'/'+bpmFileName+"_liny";
        bpmFileName +="_bpm";


        linxFile.open(linxFilePath.c_str());
        linyFile.open(linyFilePath.c_str());
        if (InpData.labelrun == 1) {
            noiseFilePath = workingDirectoryPath+'/'+bpmFileName+"_noise";
            noiseFile.open(noiseFilePath.c_str());
            if(cannotOpenFile(noiseFilePath,'o')){
                std::cerr << "Leaving drive due to error" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        if(cannotOpenFile(linxFilePath,'o') || cannotOpenFile(linyFilePath,'o')){
            std::cerr << "Leaving drive due to error" << std::endl;
            exit(EXIT_FAILURE);
        }

        linxFile << "* NAME S    BINDEX SLABEL TUNEX MUX  AMPX NOISE PK2PK AMP01 PHASE01 CO   CORMS AMP_20 PHASE_20 AMP02 PHASE02 AMP_30 PHASE_30 AMP_1_1 PHASE_1_1 AMP2_2 PHASE2_2 AMP0_2 PHASE0_2 NATTUNEX NATAMPX\n";
        linxFile << std::scientific << "$ %s  %le %le   %le   %le  %le %le %le  %le  %le  %le    %le %le  %le   %le     %le  %le    %le   %le     %le    %le      %le   %le     %le   %le     %le     %le\n";
        linyFile << "* NAME S    BINDEX SLABEL TUNEY MUY  AMPY NOISE PK2PK AMP10 PHASE10 CO   CORMS AMP_1_1 PHASE_1_1 AMP_20 PHASE_20 AMP1_1 PHASE1_1 AMP0_2 PHASE0_2 AMP0_3 PHASE0_3 NATTUNEY NATAMPY\n";
        linyFile << std::scientific << "$ %s  %le %le   %le   %le  %le %le %le  %le  %le  %le    %le %le  %le    %le      %le   %le     %le   %le     %le   %le     %le   %le     %le       %le\n";

        flag = 0;
        for (i = 0; i < MAXPICK; i++)
            BPMs[i].pickedup = false;

        /* start data file reading, puts the tbt data each BPM struct into  with all the data from the pick-ups */
        bpmCounter = 0;
        columnCounter = 0;
        horizontalBpmCounter = -1;
        verticalBpmCounter = MAXPICK / 2 - 1;
        i = 0;
        string1[0] = (char)dataFile.get();
        while (string1[0] == '#') {       /* then it is a comment line (tbach) */
            while (dataFile.get() != '\n');       /* read until the end of the line (tbach) */
            string1[0] = (char)dataFile.get();        /* read the first char of the new line (tbach) */
        }
        /* after this, we have skipped all the comment lines, and s[0] is the first character of a new line which is not a "#" (tbach) */
        if (LOG_INFO)
            printf("BPM file content:\n");
		int currentPlane = -1;
        while (string1[0] != EOF) {
            if (string1[0] == '\n') {
                ++bpmCounter;
                if (LOG_INFO)
                    printf("\n");
                columnCounter = 0;
            }
            if (isspace((int)string1[0]) && flag == 1)
                flag = 0;
            if (!isspace((int)string1[0]) && flag == 0) {
                while (!isspace((int)string1[i]) && string1[i] != EOF) {
                    ++i;
                    string1[i] = (char)dataFile.get();
                    if (i > 100) {
                        string1[i + 1] = '\0';
                        fprintf(stderr, "Found a value which has more than 100 characters, exit parsing.\n"
                            "This is most probably a malformatted file. bpmCounter=%d columnCounter=%d string1=%s\n", bpmCounter, columnCounter, string1);
                        exit(EXIT_FAILURE);
                    }
                }
                string1[i + 1] = string1[i];
                string1[i] = '\0';
                if (LOG_INFO)
                    printf("%s ", string1);
                if (columnCounter >= MAXTURNS) {
                    fprintf(stderr, "Found >= %d Turns, this turn size is not supported. Reduce amount of turns. bpmCounter:%d\n", MAXTURNS - 3, bpmCounter); /* 0,1,2 is plane, name and location (tbach) */
                    exit(EXIT_FAILURE);
                }
                if (bpmCounter >= MAXPICK) {
                    fprintf(stderr, "Found >= %d BPMs, this size is not supported. Reduce amount of BPMs. columnCounter:%d\n", MAXPICK, columnCounter);
                    exit(EXIT_FAILURE);
                }
                if (columnCounter == 0) {   /*plane (tbach) */
                	currentPlane = atoi(string1);
                    if (currentPlane == 0) /* 0 is horizontal, 1 is vertical (tbach) */
                        BPMs[++horizontalBpmCounter].plane = currentPlane;
                    else
                        BPMs[++verticalBpmCounter].plane = currentPlane;
                }

                else if (columnCounter == 1) {   /*bpm name (tbach) */
                    if (currentPlane == 0) {
                        if (horizontalBpmCounter < 0) /* Branch prediction will cry, but well lets have security (tbach) */
                        {
                            fprintf(stderr, "horizontalBpmCounter < 0. Should not happen. Probably malformatted input file?\n");
                            exit(EXIT_FAILURE);
                        }
                        BPMs[horizontalBpmCounter].bpmname = string1;
                        BPMs[horizontalBpmCounter].pickedup = true;
                    } else {
                        BPMs[verticalBpmCounter].bpmname = string1;
                        BPMs[verticalBpmCounter].pickedup = true;
                    }
                }

                else if (columnCounter == 2) {   /*bpm location (tbach) */
                    if (currentPlane == 0)
                    {
                        if (horizontalBpmCounter < 0) /* Branch prediction will cry, but well lets have security (tbach) */
                        {
                            fprintf(stderr, "horizontalBpmCounter < 0. Should not happen. Probably malformatted input file?\n");
                            exit(EXIT_FAILURE);
                        }
                        BPMs[horizontalBpmCounter].bpmpos = atof(string1);
                    }
                    else
                        BPMs[verticalBpmCounter].bpmpos = atof(string1);
                }

                else {    /*bpm data (tbach) */
                    if (currentPlane == 0)
                        BPMs[horizontalBpmCounter].tbtdata[columnCounter - 3] = atof(string1);
                    else
                        BPMs[verticalBpmCounter].tbtdata[columnCounter - 3] = atof(string1);
                    Nturns = columnCounter - 3 + 1;
                    /* If the last line is an empty line, then we can get the number of turns only from here.
                       First 3 are plane, name and location.
                       Plus 1 for index start at 0
                       (tbach) */
                }
                ++columnCounter;
                flag = 1;
                string1[0] = string1[i + 1];
                i = 0;
            }
            if (flag == 0)
                string1[0] = (char)dataFile.get();
        }
        dataFile.close();

        Nbpms = bpmCounter;
        counth = horizontalBpmCounter + 1;
        countv = verticalBpmCounter + 1;

        /* now redefine turns as the minimum of the Nturns read and the DrivingTerms data card */
        /* NB assumes all BPMs have the same number of turns as the last one read is used */
        if (InpData.turns > Nturns) InpData.turns = Nturns;

        /* Some statistics and checks */
        printf("Total number of pick-ups: %d Last turn number: %d, turns to run: %d\n", Nbpms, Nturns, InpData.turns);
        printf("Horizontal pick-ups: %d   Vertical pick-ups: %d\n", counth, -MAXPICK / 2 + countv);
        printf("name of BPM[0]: %s, pos: %f, first turn: %f, second turn: %f, last turn: %f, last turn to run: %f \n",
             BPMs[0].bpmname.c_str(), BPMs[0].bpmpos, BPMs[0].tbtdata[0], BPMs[0].tbtdata[1], BPMs[0].tbtdata[Nturns - 1], BPMs[0].tbtdata[InpData.turns - 1]);
        /* end of data file reading */

        printf("kick: %d \n", InpData.kick);
        /* searching for two working adjacent pick-ups */
        /* after the Q-kickers for finding the kick */
        if (InpData.kick < 0) {
            start = -(InpData.kcase - 1) * MAXPICK / 2 + 2;
            while (BPMs[start].pickedup == false || BPMs[start + 2].pickedup == false) {
                start = start + 2;
            }

            printf("looking for kick in pick-up:%d\n", start + 1);
            /* Find kick here and get kick */
            for (columnCounter = 1; (InpData.kick < 0) && (columnCounter < InpData.turns); ++columnCounter) {
                if (fabs(BPMs[start].tbtdata[columnCounter] - BPMs[start].tbtdata[columnCounter - 1]) > InpData.kper) {
                    InpData.kick = columnCounter;
                }
            }

            if (InpData.kick < 0) {
                fprintf(stderr, "NO KICK FOUND\n");
                exit(EXIT_FAILURE);
            } else
                printf("Found kick in turn:%d\n", InpData.kick + 1);    /*Natural count */
        }

        if (InpData.kick > 0) {
            for (i = 0; i < MAXPICK; i++) {
                if (BPMs[i].pickedup == true) {
                    for (j = InpData.kick; j < InpData.turns; j++)
                        BPMs[i].tbtdata[j - InpData.kick] = BPMs[i].tbtdata[j];
                }
            }
            InpData.turns -= InpData.kick;
        }
        printf("Turns to be processed after kick offset: %d BPMs[0].tbtdata[0]: %f \n", InpData.turns, BPMs[0].tbtdata[0]);

        /* First part of the analysis: Determine  phase of all pick-ups and noise */
        writeSussixInput(sussixInputFilePath, InpData.turns, InpData.istun, InpData.tunex, InpData.tuney);

        if (counth >= (countv - MAXPICK / 2))
            maxcounthv = counth;
        else
            maxcounthv = -MAXPICK / 2 + countv;

        if (maxcounthv > InpData.pickend)
            maxcounthv = InpData.pickend;

        if (maxcounthv >= MAXPICK) {
            fprintf(stderr, "\nNot enough Pick-up mexmory\n");
            exit(EXIT_FAILURE);
        }
        printf("BPMs in loop: %d, pickstart: %d, resulting loop length: %d\n",
             maxcounthv, InpData.pickstart, maxcounthv - InpData.pickstart);

#pragma omp parallel for private(i, horizontalBpmCounter, verticalBpmCounter, j, maxamp, calculatednattunex, calculatednattuney, calculatednatampx, calculatednatampy)
        for (i = InpData.pickstart; i < maxcounthv; ++i) {
            horizontalBpmCounter = i;
            verticalBpmCounter = i + MAXPICK / 2;

            if (verticalBpmCounter >= countv)
                verticalBpmCounter = countv - 1;
            if (horizontalBpmCounter >= counth)
                horizontalBpmCounter = counth - 1;
            if (horizontalBpmCounter < 0 || verticalBpmCounter < 0)
            {
                fprintf(stderr, "horizontal or vertical BpmCounter < 0. Should not happen.\n");
                exit(EXIT_FAILURE);
            }

            for (j = 0; j < MAXTURNS; ++j) {
                doubleToSend[j] = BPMs[horizontalBpmCounter].tbtdata[j];
                doubleToSend[j + MAXTURNS] = BPMs[verticalBpmCounter].tbtdata[j];
                doubleToSend[j + 2 * MAXTURNS] = 0.0;
                doubleToSend[j + 3 * MAXTURNS] = 0.0;
            }


            /* This calls the external Fortran code (tbach)-Different name depending on OS (asherman)*/
            #ifdef _WIN32
                SUSSIX4DRIVENOISE (&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0], (char *)sussixInputFilePath.c_str());
            #else
                sussix4drivenoise_(&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0], (char *)sussixInputFilePath.c_str());
            #endif
            /* Let's look for natural tunes in the istun range if natural tunes input is given*/
            maxamp = 0;
            calculatednattunex = NATTUNE_DEFAULT;
            if (InpData.nattunex > NATTUNE_DEFAULT) {
                for (j = 0; j < 300; ++j) {
                    if ((InpData.nattunex - InpData.istun < allfreqsx[j] && allfreqsx[j] < InpData.nattunex + InpData.istun) && (maxamp < allampsx[j])) {
                        maxamp = allampsx[j];
                        calculatednattunex = allfreqsx[j];
                        calculatednatampx = maxamp;
                    }
                }
            }
            maxamp = 0;
            calculatednattuney = NATTUNE_DEFAULT;
            if (InpData.nattuney > NATTUNE_DEFAULT) {
                for (j = 0; j < 300; ++j) {
                    if ((InpData.nattuney - InpData.istun < allfreqsy[j] && allfreqsy[j] < InpData.nattuney + InpData.istun) && (maxamp < allampsy[j])) {
                        maxamp = allampsy[j];
                        calculatednattuney = allfreqsy[j];
                        calculatednatampy = maxamp;
                    }
                }
            }

            #pragma omp critical
            {
                BPMs[horizontalBpmCounter].pickedup = BPMstatus(1, InpData.turns); /*Always returns true*/
                if (InpData.labelrun == 1)
                    noiseFile << std::scientific << "1 " << horizontalBpmCounter << "  " <<  noise1 << ' ' <<  noiseAve << ' ' << maxpeak << ' ' << maxfreq << ' ' << maxmin << ' ' << nslines << ' ' << BPMs[i].pickedup << ' ' << phase[0] / 360. << std::endl;

                /* PRINT LINEAR FILE */
                if (amplitude[0] > 0 && BPMs[i].pickedup == true && horizontalBpmCounter == i) {
                    linxFile <<  std::scientific << '"' << BPMs[horizontalBpmCounter].bpmname << "\" " << BPMs[horizontalBpmCounter].bpmpos << ' ' << horizontalBpmCounter << ' ' << BPMs[horizontalBpmCounter].pickedup << ' ' << tune[0] << ' ' <<
                            phase[0] / 360. << ' ' << amplitude[0] << ' ' << noise1 << ' ' << maxmin << ' ' << amplitude[2] / amplitude[0] << ' ' << phase[2] / 360. << ' ' << co << ' ' << co2 << ' ' << amplitude[1] / amplitude[0] << ' ' <<
                            phase[1] / 360. << ' ' << amplitude[12] / amplitude[0] << ' ' << phase[12] / 360. << ' ' << amplitude[6] / amplitude[0] << ' ' <<
                            phase[6] / 360. << ' ' << amplitude[14] / amplitude[0]  << ' ' << phase[14] / 360. << ' ' << amplitude[16] / amplitude[0] << ' ' <<
                            phase[16] / 360. << ' ' << amplitude[18] / amplitude[0] << ' ' << phase[18] / 360. << ' ' << calculatednattunex << ' ' << calculatednatampx << std::endl;

                    ++CalcData.tunecountx;
                    CalcData.tunesumx += tune[0];
                    CalcData.tune2sumx += tune[0] * tune[0];
                    if (calculatednattunex > NATTUNE_DEFAULT) { /*  Initialized to -100. Condition true if nat tune found */
                        ++CalcData.nattunexcount;
                        CalcData.nattunexsum += calculatednattunex;
                        CalcData.nattunex2sum += calculatednattunex * calculatednattunex;
                    }

                    /* Horizontal Spectrum output */
                    if (i < 10) {
                        spectrumFilePath = workingDirectoryPath+'/'+BPMs[i].bpmname+".x";
                        spectrumFile.open(spectrumFilePath.c_str());
                        if(cannotOpenFile(spectrumFilePath,'o')){
                            std::cout << "Leaving drive due to error" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        spectrumFile << "* FREQ AMP\n$ %le %le\n";
                        for (j = 0; j < 300; ++j)
                            spectrumFile << std::scientific << allfreqsx[j] << ' ' << allampsx[j] << std::endl;
                        spectrumFile.close();
                    }
                }
                BPMs[verticalBpmCounter].pickedup = BPMstatus(2, InpData.turns); /*Always returns true*/
                if (InpData.labelrun == 1)
                    noiseFile << std::scientific << "2 " << verticalBpmCounter << "  " <<  noise1 << ' ' <<  noiseAve << ' ' << maxpeak << ' ' << maxfreq << ' ' << maxmin << ' ' << nslines << ' ' << BPMs[verticalBpmCounter].pickedup << ' ' << phase[3] / 360. << std::endl;
                if (amplitude[3] > 0 && BPMs[verticalBpmCounter].pickedup == true && verticalBpmCounter == i + MAXPICK / 2) {
                    linyFile <<  std::scientific << '"' << BPMs[verticalBpmCounter].bpmname << "\" " << BPMs[verticalBpmCounter].bpmpos << ' ' << verticalBpmCounter << ' ' << BPMs[verticalBpmCounter].pickedup << ' ' << tune[1] << ' ' <<
                            phase[3] / 360. << ' ' << amplitude[3] << ' ' << noise1 << ' ' << maxmin << ' ' << amplitude[5] / amplitude[3] << ' ' << phase[5] / 360. << ' ' << co << ' ' << co2 << ' ' <<
                            amplitude[13] / amplitude[3] << ' ' << phase[13] / 360. << ' ' << amplitude[15] / amplitude[3] << ' ' << phase[15] / 360. << ' ' <<
                            amplitude[17] / amplitude[3] << ' ' << phase[17] / 360. << ' ' << amplitude[4] / amplitude[3] << ' ' << phase[4] / 360. << ' ' <<
                            amplitude[11] / amplitude[3] << ' ' << phase[11] / 360. << ' ' << calculatednattuney << ' ' << calculatednatampy << std::endl;
                    ++CalcData.tunecounty;
                    CalcData.tunesumy += tune[1];
                    CalcData.tune2sumy += tune[1] * tune[1];
                    if (calculatednattuney > NATTUNE_DEFAULT) { /*  Initialized to -100. Condition true if nat tune found */
                        ++CalcData.nattuneycount;
                        CalcData.nattuneysum += calculatednattuney;
                        CalcData.nattuney2sum += calculatednattuney * calculatednattuney;
                    }
                    if (verticalBpmCounter < MAXPICK / 2 + 10) {
                        spectrumFilePath = workingDirectoryPath+'/'+BPMs[verticalBpmCounter].bpmname+".y";
                        spectrumFile.open(spectrumFilePath.c_str());
                        if(cannotOpenFile(spectrumFilePath,'o')){
                            std::cout << "Leaving drive due to error" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        spectrumFile << "* FREQ AMP\n$ %le %le\n";
                        for (j = 0; j < 300; ++j)
                            spectrumFile << std::scientific << allfreqsy[j] << ' ' << allampsy[j] << std::endl;
                        spectrumFile.close();
                    }
                }
            } /* end of omp critical section */
        } /* end of parallel for */
        linxFile.close();
        linyFile.close();
        if (InpData.labelrun == 1) noiseFile.close();
        /* Sort and move the "@..." lines to the top of the _linx/y files */
        formatLinFile(linxFilePath,
                CalcData.tunecountx, CalcData.tunesumx, CalcData.tune2sumx, CalcData.nattunexcount, CalcData.nattunexsum, CalcData.nattunex2sum, 1);
        formatLinFile(linyFilePath,
                CalcData.tunecounty, CalcData.tunesumy, CalcData.tune2sumy, CalcData.nattuneycount, CalcData.nattuneysum, CalcData.nattuney2sum, 2);
    } /* end of while loop over all files to analyse */
    drivingTermsFile.close();

    if (loopcounter == 0)
        std::cout << "Drivingterms file has bad input, no data ever read\n";

    return EXIT_SUCCESS;
}

void get_inp_data(std::string driveInputFilePath){

    std::string temp_str;
    std::ifstream driveInputFile;
    unsigned int pos;

    driveInputFile.open(driveInputFilePath.c_str());

    while(!driveInputFile.rdstate()){
        std::getline (driveInputFile, temp_str);
        if((pos = temp_str.find("KICK=")) != std::string::npos) InpData.kick = atoi(temp_str.substr(pos+strlen("KICK=")).c_str())-1;
        if((pos = temp_str.find("CASE(1[H], 0[V])=")) != std::string::npos) InpData.kcase = atoi(temp_str.substr(pos+strlen("CASE(1[H], 0[V])=")).c_str());
        if((pos = temp_str.find("KPER(KICK PERCE.)=")) != std::string::npos) InpData.kper = atof(temp_str.substr(pos+strlen("KPER(KICK PERCE.)=")).c_str());
        if((pos = temp_str.find("TUNE X=")) != std::string::npos) InpData.tunex = atof(temp_str.substr(pos+strlen("TUNE X=")).c_str());
        if((pos = temp_str.find("TUNE Y=")) != std::string::npos) InpData.tuney = atof(temp_str.substr(pos+strlen("TUNE Y=")).c_str());
        if((pos = temp_str.find("PICKUP START=")) != std::string::npos) InpData.pickstart = atoi(temp_str.substr(pos+strlen("PICKUP START=")).c_str());
        if((pos = temp_str.find("PICKUP END=")) != std::string::npos) InpData.pickend = atoi(temp_str.substr(pos+strlen("PICKUP END=")).c_str());
        if((pos = temp_str.find("ISTUN=")) != std::string::npos) InpData.istun = atof(temp_str.substr(pos+strlen("ISTUN=")).c_str());
        if((pos = temp_str.find("LABEL RUN (1[yes])=")) != std::string::npos) InpData.labelrun = atoi(temp_str.substr(pos+strlen("LABEL RUN (1[yes])=")).c_str());
        if((pos = temp_str.find("WINDOWa1=")) != std::string::npos) InpData.windowa1 = atof(temp_str.substr(pos+strlen("WINDOWa1=")).c_str());
        if((pos = temp_str.find("WINDOWa2=")) != std::string::npos) InpData.windowa2 = atof(temp_str.substr(pos+strlen("WINDOWa2=")).c_str());
        if((pos = temp_str.find("WINDOWb1=")) != std::string::npos) InpData.windowb1 = atof(temp_str.substr(pos+strlen("WINDOWb1=")).c_str());
        if((pos = temp_str.find("WINDOWb2=")) != std::string::npos) InpData.windowb2 = atof(temp_str.substr(pos+strlen("WINDOWb2=")).c_str());
        if((pos = temp_str.find("NATURAL X=")) != std::string::npos) InpData.nattunex = atof(temp_str.substr(pos+strlen("NATURAL X=")).c_str());
        if((pos = temp_str.find("NATURAL Y=")) != std::string::npos) InpData.nattuney = atof(temp_str.substr(pos+strlen("NATURAL Y=")).c_str());

    }
    driveInputFile.close();
    /* input/option file reading end */
}
void get_and_check_file_paths(std::string *workingDirectoryPath, std::string *drivingTermsFilePath, std::string *driveInputFilePath, std::string *sussixInputFilePath, const char * cmdinput){

    *workingDirectoryPath = cmdinput;

    std::cout << *workingDirectoryPath << std::endl;

    if(cannotOpenFile(*workingDirectoryPath,'i') && OS == "linux"){ //Always fails to open in windows
        std::cout << "Leaving drive due to error" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "\nWorking directory: " << *workingDirectoryPath << std::endl;


    *drivingTermsFilePath = *workingDirectoryPath+"/DrivingTerms";
    *driveInputFilePath = *workingDirectoryPath+"/Drive.inp";
    *sussixInputFilePath = *workingDirectoryPath+"/sussix_v4.inp";

    //check the input files drivingTermsFilePath and Drive.inp

    if(cannotOpenFile(*drivingTermsFilePath,'i')
    || cannotOpenFile(*driveInputFilePath,'i') || cannotOpenFile(*sussixInputFilePath,'o')){
        std::cout << "Leaving drive due to error" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "DrivingTerms: " <<  *drivingTermsFilePath << std::endl;
    std::cout << "Drive.inp: " <<  *driveInputFilePath << std::endl;
    std::cout << "sussix_v4.inp: " <<  *sussixInputFilePath << std::endl;

}

 void check_inp_data(){

    InpData.check_tune(InpData.tunex,'x');
    InpData.check_tune(InpData.tuney,'y');

    if (InpData.kick >= 0)
        printf("Known kick in turn %d\n", InpData.kick + 1);
    if (InpData.kcase == 1)
        printf("Horizontal case\n");
    else if (InpData.kcase == 0)
        printf("Vertical case\n");
    else {
        fprintf(stderr, "No proper kcase in Drive.inp\n");
        exit(EXIT_FAILURE);
    }

    if (InpData.labelrun == 1)
        printf("\n LABELRUN: NOISE FILES WILL BE WRITTEN TO NOISEPATH\n");

    printf("pickstart: %d, pickend: %d\n", InpData.pickstart, InpData.pickend);
    if (InpData.pickstart < 0 || InpData.pickstart > InpData.pickend || InpData.pickstart > MAXPICK) {
        fprintf(stderr, "Bad value for pickstart. Must be >= 0 and < pickend and <= MAXPICK(=%d)\n", MAXPICK);
        exit(EXIT_FAILURE);
    }
}
/* *****************  */
/*    readDrivingTerms*/
/* *****************  */
int readDrivingTerms(std::istream& drivingTermsFile, int *turns, std::string *path)
{
    /* This functions reads from the given FILE one line of this format:
     * <path to datafile> <int start turn> <int end turn>
     * If a line was successfully read, return true
     * If not, return false
     * --tbach */

    int num;
    unsigned int pos;
    *path = "";
    *turns = 0;

    drivingTermsFile >> *path >> num >> *turns;

    if(*path == "" || *turns == 0)
        return 0;

    if(OS == "Windows32") //Path may not exist if drivingtermsfile was created on linux, so this deals with that. (Assumes /afs/cern.ch is 'H' drive)
        if((pos = (*path).find("/afs/cern.ch")) != std::string::npos)
            *path = "H:"+(*path).substr(pos+strlen("/afs/cern.ch"));

    if(OS == "linux"){ //Same story as above, replaces H: with /afs/cern.ch, also changes '\' to '/'
        while((pos = (*path).find("\\")) != std::string::npos)
            (*path)[pos] = '/';
        if((pos = (*path).find("H:")) != std::string::npos)
            *path = "/afs/cern.ch"+(*path).substr(pos+strlen("H:"));
    }

        return 1;
}

/* ***************** */
/*    sussix_inp     */
/* ***************** */
void writeSussixInput(std::string sussixInputFilePath, const int turns, const double istun, const double tunex, const double tuney)
{
    std::ofstream sussixInputFile(sussixInputFilePath.c_str());

    if(cannotOpenFile(sussixInputFilePath,'o')){
        std::cout << "Leaving drive due to error: " << std::endl;
        exit(EXIT_SUCCESS);
    }
    sussixInputFile << "C" << std::endl << "C INPUT FOR SUSSIX_V4 ---17/09/1997---" << std::endl << "C DETAILS ARE IN THE MAIN PROGRAM SUSSIX_V4.F\n";
    sussixInputFile << "C" << std::endl << std::endl;
    sussixInputFile << "ISIX  = 0" << std::endl << "NTOT  = 1" << std::endl << "IANA  = 1" << std::endl << "ICONV = 0" << std::endl;
    sussixInputFile << std::scientific << "TURNS = 1 " << turns << std::endl << "NARM  = 300" << std::endl << "ISTUN = 1 " << istun << ' ' << istun << std::endl;
    sussixInputFile << std::scientific << "TUNES = " << tunex << ' ' << tuney << " .07" << std::endl << "NSUS  = 0" << std::endl << "IDAM  = 2" << std::endl << "NTWIX = 1" << std::endl;
    sussixInputFile << "IR    = 1" << std::endl << "IMETH = 2" << std::endl << "NRC   = 4" << std::endl << "EPS   = 2D-3" << std::endl; /* EPS is the window in the secondary lines, very imp!!! */
    sussixInputFile << "NLINE = 0" << std::endl << "L,M,K = " << std::endl << "IDAMX = 1" << std::endl << "NFIN  = 500" << std::endl;
    sussixInputFile << "ISME  = 0" << std::endl << "IUSME = 200" << std::endl << "INV   = 0" << std::endl << "IINV  = 250" << std::endl;
    sussixInputFile << "ICF   = 0" << std::endl << "IICF  = 350" << std::endl;
    sussixInputFile.close();
}

/************   BPMstatus *************/
/* Analyse fort.300 to detect noise   */
/**************************************/

/*UPDATE: All rejections are removed from this code, as other parts
 * of Beta-beat are in charge of rejecting now.  Calculations are still
 * kept in place. --asherman (07/2013)
 */
#define MINSIGNAL 0.00001
#define SIGMACUT   1.8
#define MINIMUMNOISE 0.0
#define BADPICKUP  8.0
#define MAXSIGNAL 30000
bool BPMstatus(const int plane, const int turns)
{
    double aux = 0, ave = 0, amp = 0,
        maxe = -500000.0, mine = 500000.0;
    int il,counter,counter3 = 0;

    maxpeak = 0;                /*Initialising */
    co = 0.0;
    co2 = 0.0;

    if (plane == 1) {
        for (il = 0; il < turns; il++) {
            co += doubleToSend[il];
            co2 += doubleToSend[il] * doubleToSend[il];
            if (doubleToSend[il] < mine)
                mine = doubleToSend[il];
            if (doubleToSend[il] > maxe)
                maxe = doubleToSend[il];
        }
    }
    if (plane == 2) {
        for (il = MAXTURNS; il < MAXTURNS + turns; il++) {
            co += doubleToSend[il];
            co2 += doubleToSend[il] * doubleToSend[il];
            if (doubleToSend[il] < mine)
                mine = doubleToSend[il];
            if (doubleToSend[il] > maxe)
                maxe = doubleToSend[il];
        }
    }
    co = co / turns;
    co2 = sqrt(co2 / turns - co * co);
    maxmin = maxe - mine;

    /* Compute the spread and average in the intervals [windowa1,windowa2]
       and [windowb1,windowb2] */

    noise1 = 0;

    for (counter = 0; counter < 300; counter++) {
        if (plane == 1) {
            aux = allfreqsx[counter];
            amp = allampsx[counter];
        }
        else if (plane == 2) {
            aux = allfreqsy[counter];
            amp = allampsy[counter];
        }

        if (amp > maxpeak && aux > 0.05) {
            maxpeak = amp;
            maxfreq = aux;
        }

        if ((aux < InpData.windowa2 && aux > InpData.windowa1)
            || (aux < InpData.windowb2 && aux > InpData.windowb1)) {
            if (amp < 0) {      /* Something in sussix went wrong */
                noise1 = 100;
                noiseAve = 100;
                maxpeak = 100;
                maxfreq = 100;
                return true;
            }

            ave = amp + ave;
            noise1 = noise1 + amp * amp;
            ++counter3;
        }

    }
    if (counter3 > 0) {
        if (counter3 > 1)
            noise1 = sqrt((noise1 / counter3 - ave * ave / (counter3*counter3)));
        else
            noise1 = 0;
        noiseAve = ave / counter3;
    } else {
        noise1 = MINIMUMNOISE;
        noiseAve = MINIMUMNOISE;
    }
    nslines = counter3;

    if ((InpData.windowa1 < maxfreq && maxfreq < InpData.windowa2)
     || (InpData.windowb1 < maxfreq && maxfreq < InpData.windowb2))
        printf("NoiseWindow includes largest lines, amp %e freq %e!!!!\n",
               maxpeak, maxfreq);

    return true;
}

bool cannotOpenFile(std::string filePath,char type){
    bool failure;
    if(type == 'o')
        failure = outputFileCheck(filePath);
    else
        failure = inputFileCheck(filePath);
    if(failure){
        std::cout << "Failed to open file " << filePath << std::endl;
        return failure;
        }
    else
        return failure;
}

bool inputFileCheck(std::string filePath){
    std::ifstream file(filePath.c_str());
    bool failure = file.fail();
    file.close();
    return failure;
}

bool outputFileCheck(std::string filePath){
    std::ofstream file(filePath.c_str());
    bool failure = file.fail();
    file.close();
    return failure;
}

void formatLinFile(std::string linFilePath,
        const int tunecount, const double tunesum, const double tune2sum, const int nattunecount, const double nattunesum, const double nattune2sum, int plane_index) {

    std::string tempFilePath, temp_str, fileHolder[MAXPICK];
    std::ifstream linFile;
    std::ofstream tempFile;
    int max=0, min = 1300, bIndex;
    unsigned int pos;
    double temp_double;

    if (tunecount > 0) {
        tempFilePath = linFilePath+"_temp";
        tempFile.open(tempFilePath.c_str());
        if(cannotOpenFile(tempFilePath,'o')){
            std::cout << "Leaving drive due to error" << std::endl;
            exit(EXIT_FAILURE);
        }
        tempFile << std::scientific << "@ Q" << plane_index << " %le " << tunesum / tunecount << "\n@ Q" << plane_index << "RMS %le " << sqrt(tune2sum / tunecount - (tunesum / tunecount) * (tunesum / tunecount)) << std::endl;
        if (nattunecount > 0)
            tempFile << std::scientific << "@ NATQ" << plane_index << " %le " << nattunesum / nattunecount << "\n@ NATQ" << plane_index << "RMS %le " << sqrt(nattune2sum / nattunecount - (nattunesum / nattunecount) * (nattunesum / nattunecount)) << std::endl;

        /*Gets linFile to read*/
        linFile.open(linFilePath.c_str());
        if(cannotOpenFile(linFilePath,'i')){
            std::cout << "Leaving drive due to error" << std::endl;
            exit(EXIT_FAILURE);
        }

        /*Writes first two lines from lineFile into tempFile.*/
        getline(linFile,temp_str);
        tempFile << temp_str << std::endl;
        getline(linFile,temp_str);
        tempFile << temp_str << std::endl;

        /*Begins loop which simultaneously places lines of linFile
         * into fileHolder and sorts them according to bIndex which
         * is an integer in the third column of each line.
         */
        while(!linFile.rdstate()){
            pos = linFile.tellg();
            linFile >> temp_str >> temp_double >> bIndex;
            linFile.seekg(pos);
            getline(linFile, temp_str);
            if(temp_str[0] == '"'){
                fileHolder[bIndex] = temp_str;
                if(bIndex > max) /*Finds highest bIndex in linFile*/
                    max = bIndex;
                if(bIndex < min) /*Finds lowest bIndex in linFile*/
                    min = bIndex;
                }
        }

        /*Writes sorted lines into tempfile*/
       for(bIndex = min;bIndex <= max;bIndex++){
            if(fileHolder[bIndex].size() != 0)
                tempFile << fileHolder[bIndex] << std::endl;
        }

        /*Removes linFile, renames tempFile to linFile*/
        tempFile.close();
        linFile.close();
        remove(linFilePath.c_str());
        rename(tempFilePath.c_str(), linFilePath.c_str());
        std::cout <<  linFilePath << ":\ntune: \n  sum: " << tunesum << ", count: " << tunecount << ", sum2: " << tune2sum << ", \nnatural tune: \n  sum:" << nattunesum << ", count:" << nattunecount << ", sum2:" << nattune2sum << std::endl;
    }
}
