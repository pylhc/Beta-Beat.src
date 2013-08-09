/* Parallelised version of Drive_God_lin.c of 15/05/2011 using Openmp.
   NB this version uses as the number of turns the minimum of the number
   actually read or the number given in the DrivingTerms file.

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
   
   Change 22/07/2013: -Certain local main variables grouped in a struct called 'Data'
   -Removed any rejections in BPMstatus function
   -Changed formatLinFile to be OS compatible, makefile also modified, drive can now 
   be compiled using 'make' on both linux and windows
   -In general, the makefile and code only distinguish between windows and non-windows 
   with the hopes that other OS's will work with the linux version.  If this turns out 
   to not be the case, the code can be modified to allow for another OS.
   -Removed some unused and unnecessary code  
    05/08: Forced input tune to be in [-0.5,0.5] by taking closes member in this interval
    -Updated code to use c++ string manipulations     -asherman
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

int readDrivingTerms(std::istream&, int*, std::string *);
void writeSussixInput(std::string, const int, const double, const double, const double);
int BPMstatus(const int, const int);
void formatLinFile(std::string, const int, const double, const double, const int, const double, const double,int);
bool cannotOpenFile(std::string,char);
bool inputFileCheck(std::string);
bool outputFileCheck(std::string);

    double calculatednattuney, calculatednattunex, calculatednatampy, calculatednatampx, co, co2, maxamp,
           windowa1, windowa2, windowb1, windowb2, noise1, maxfreq, maxmin, maxpeak, nattunex, nattuney, noiseAve;

    double allampsx[300], allampsy[300], allfreqsx[300], allfreqsy[300], amplitude[19], bpmpos[MAXPICK],
               doubleToSend[MAXTURNS4 + 4], matrix[MAXPICK][MAXTURNS], phase[19],
               tune[2];

    int labelrun, nslines=0, Nturns=0;

    int label[MAXPICK], hv[MAXPICK], hvt[MAXPICK];

#pragma omp threadprivate(amplitude, doubleToSend, tune, phase,\
        noise1, noiseAve, maxpeak, maxfreq, maxmin, co, co2,\
        allfreqsx, allampsx, allfreqsy, allampsy)

int main(int argc, char **argv)
{

    struct Data_Struct{
        size_t charCounter;

        double  istun, kper, tunex, tuney, nattunexsum, nattunex2sum, nattuneysum,
                nattuney2sum, tunesumx, tunesumy, tune2sumx, tune2sumy;

        int     counth, countv, kcase, kick, maxcounthv, Nbpms,
                pickstart, pickend, start, turns,
                tunecountx, tunecounty, nattunexcount, nattuneycount;
    };

    char    string1[1000];


    struct Data_Struct Data;


    int i, j, bpmCounter, columnCounter, horizontalBpmCounter, verticalBpmCounter, flag, loopcounter=0;

    char* bpmname[MAXPICK];
    
    std::ofstream linxFile, linyFile, noiseFile, spectrumFile;
    std::ifstream driveInputFile, drivingTermsFile, dataFile;       
    std::string temp_str, dataFilePath, bpmFileName, workingDirectoryPath, sussixInputFilePath, driveInputFilePath, drivingTermsFilePath, noiseFilePath, linxFilePath, linyFilePath, spectrumFilePath;  
   
    unsigned int pos;
    
    #ifdef _WIN32 /*Changes minor formatting difference in windows regarding the output of a number in scientific notation.*/
        _set_output_format(_TWO_DIGIT_EXPONENT);
    #endif  
    
    omp_set_dynamic(0);
    /* Memory allocation */

    for (i = 0; i < MAXPICK; i++)
        bpmname[i] = (char *) calloc(50, sizeof(char));
    
    //To output scientific notation
    std::cout << std::setiosflags (std::ios::scientific);
    
    /*  Path to DrivingTerms and Drive.inp */
    
    workingDirectoryPath = argv[1];
    
    std::cout << workingDirectoryPath << std::endl;
    
    if(cannotOpenFile(workingDirectoryPath,'i') && OS == "linux"){ //Always fails to open in windows
        std::cout << "Leaving drive due to error" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "\nWorking directory: " << workingDirectoryPath << std::endl;
    
        
    drivingTermsFilePath = workingDirectoryPath+"/DrivingTerms";
    driveInputFilePath = workingDirectoryPath+"/Drive.inp";
    sussixInputFilePath = workingDirectoryPath+"/sussix_v4.inp";  

    //check the input files drivingTermsFilePath and Drive.inp
    
    if(cannotOpenFile(drivingTermsFilePath,'i') 
    || cannotOpenFile(driveInputFilePath,'i') || cannotOpenFile(sussixInputFilePath,'o')){
        std::cout << "Leaving drive due to error" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "DrivingTerms: " <<  drivingTermsFilePath << std::endl;
    std::cout << "Drive.inp: " <<  driveInputFilePath << std::endl;
    std::cout << "sussix_v4.inp: " <<  sussixInputFilePath << std::endl;

    /* set all options to defaults, it could happen that they are not included in the file (tbach) */
    Data.kick = Data.kcase = Data.pickstart = Data.pickend = labelrun = 0;
    Data.kper = Data.tunex = Data.tuney = Data.istun = windowa1 = windowa2 = windowb1 = windowb2 = 0.0;
    nattunex = nattuney = NATTUNE_DEFAULT;
    
    /* input/option file reading start */
    driveInputFile.open(driveInputFilePath.c_str());
    
    while(!driveInputFile.rdstate()){
        std::getline (driveInputFile, temp_str);  
        if((pos = temp_str.find("KICK=")) != std::string::npos) Data.kick = atoi(temp_str.substr(pos+strlen("KICK=")).c_str())-1;
        if((pos = temp_str.find("CASE(1[H], 0[V])=")) != std::string::npos) Data.kcase = atoi(temp_str.substr(pos+strlen("CASE(1[H], 0[V])=")).c_str());
        if((pos = temp_str.find("KPER(KICK PERCE.)=")) != std::string::npos) Data.kper = atof(temp_str.substr(pos+strlen("KPER(KICK PERCE.)=")).c_str());
        if((pos = temp_str.find("TUNE X=")) != std::string::npos) Data.tunex = atof(temp_str.substr(pos+strlen("TUNE X=")).c_str());
        if((pos = temp_str.find("TUNE Y=")) != std::string::npos) Data.tuney = atof(temp_str.substr(pos+strlen("TUNE Y=")).c_str());
        if((pos = temp_str.find("PICKUP START=")) != std::string::npos) Data.pickstart = atoi(temp_str.substr(pos+strlen("PICKUP START=")).c_str());
        if((pos = temp_str.find("PICKUP END=")) != std::string::npos) Data.pickend = atoi(temp_str.substr(pos+strlen("PICKUP END=")).c_str());
        if((pos = temp_str.find("ISTUN=")) != std::string::npos) Data.istun = atof(temp_str.substr(pos+strlen("ISTUN=")).c_str());
        if((pos = temp_str.find("LABEL RUN (1[yes])=")) != std::string::npos) labelrun = atoi(temp_str.substr(pos+strlen("LABEL RUN (1[yes])=")).c_str());
        if((pos = temp_str.find("WINDOWa1=")) != std::string::npos) windowa1 = atof(temp_str.substr(pos+strlen("WINDOWa1=")).c_str());
        if((pos = temp_str.find("WINDOWa2=")) != std::string::npos) windowa2 = atof(temp_str.substr(pos+strlen("WINDOWa2=")).c_str());
        if((pos = temp_str.find("WINDOWb1=")) != std::string::npos) windowb1 = atof(temp_str.substr(pos+strlen("WINDOWb1=")).c_str());
        if((pos = temp_str.find("WINDOWb2=")) != std::string::npos) windowb2 = atof(temp_str.substr(pos+strlen("WINDOWb2=")).c_str());
        if((pos = temp_str.find("NATURAL X=")) != std::string::npos) nattunex = atof(temp_str.substr(pos+strlen("NATURAL X=")).c_str());
        if((pos = temp_str.find("NATURAL Y=")) != std::string::npos) nattuney = atof(temp_str.substr(pos+strlen("NATURAL Y=")).c_str());
     
    }
    driveInputFile.close();
    /* input/option file reading end */

    if(Data.tunex > 0.5 || Data.tunex < -0.5){
        i=0;
        while(Data.tunex > 0.5){
            Data.tunex -= 1;
            i--;
        }
        while(Data.tunex < -0.5){
            Data.tunex += 1;
            i++;
        }    
        std::cout << "tune_x input increased by " << i << ". should be less than or equal to 0.5 in magnitude. New value is: " << Data.tunex << std::endl;
    }
    if(Data.tuney > 0.5 || Data.tuney < -0.5){
        i=0;
        while(Data.tuney > 0.5){
            Data.tuney -= 1;
            i--;
        }
        while(Data.tuney < -0.5){
            Data.tuney += 1;
            i++;
        }        
        std::cout << "tune_y input increased by " << i << ". should be less than or equal to 0.5 in magnitude. New value is: " << Data.tuney << std::endl;
    }
    
    if (Data.kick >= 0)
        printf("Known kick in turn %d\n", Data.kick + 1);
    if (Data.kcase == 1)
        printf("Horizontal case\n");
    else if (Data.kcase == 0)
        printf("Vertical case\n");
    else {
        fprintf(stderr, "No proper kcase in Drive.inp\n");
        exit(EXIT_FAILURE);
    }

    if (labelrun == 1)
        printf("\n LABELRUN: NOISE FILES WILL BE WRITTEN TO NOISEPATH\n");
    printf("pickstart: %d, pickend: %d\n", Data.pickstart, Data.pickend);
    if (Data.pickstart < 0 || Data.pickstart > Data.pickend || Data.pickstart > MAXPICK) {
        fprintf(stderr, "Bad value for pickstart. Must be >= 0 and < Data.pickend and <= MAXPICK(=%d)\n", MAXPICK);
        exit(EXIT_FAILURE);
    }


    drivingTermsFile.open(drivingTermsFilePath.c_str());
    if(cannotOpenFile(drivingTermsFilePath,'i')){
        std::cout << "Leaving drive due to error" << std::endl;
        exit(EXIT_FAILURE);
    }
    Data.turns = 0;
    
    /* From drivingTermsFilePath assign dataFilePath, assign turns. */
    while (readDrivingTerms(drivingTermsFile, &(Data.turns), &dataFilePath)) {
        /* set all values to be calculated to default values */     
        Data.tunecountx = Data.tunecounty = Data.nattunexcount = Data.nattuneycount = 0;
        Data.tunesumx = Data.tunesumy = Data.tune2sumx = Data.tune2sumy = Data.nattunexsum = Data.nattunex2sum= Data.nattuneysum = Data.nattuney2sum = 0.0;

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
        if (labelrun == 1) {
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
            label[i] = 0;

        /* start data file reading, constructing a matrix with all the data from the pick-ups */
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
                    hv[bpmCounter] = atoi(string1);
                    if (hv[bpmCounter] == 0) /* 0 is horizontal, 1 is vertical (tbach) */
                        ++horizontalBpmCounter;
                    else
                        ++verticalBpmCounter;
                }

                else if (columnCounter == 1) {   /*bpm name (tbach) */
                    if (hv[bpmCounter] == 0) {
                        if (horizontalBpmCounter < 0) /* Branch prediction will cry, but well lets have security (tbach) */
                        {
                            fprintf(stderr, "horizontalBpmCounter < 0. Should not happen. Probably malformatted input file?\n");
                            exit(EXIT_FAILURE);
                        }
                        hvt[horizontalBpmCounter] = 0;
                        strcpy(bpmname[horizontalBpmCounter], string1);
                        label[horizontalBpmCounter] = 1;
                    } else {
                        hvt[verticalBpmCounter] = 1;
                        strcpy(bpmname[verticalBpmCounter], string1);
                        label[verticalBpmCounter] = 1;
                    }
                }

                else if (columnCounter == 2) {   /*bpm location (tbach) */
                    if (hv[bpmCounter] == 0)
                    {
                        if (horizontalBpmCounter < 0) /* Branch prediction will cry, but well lets have security (tbach) */
                        {
                            fprintf(stderr, "horizontalBpmCounter < 0. Should not happen. Probably malformatted input file?\n");
                            exit(EXIT_FAILURE);
                        }
                        bpmpos[horizontalBpmCounter] = atof(string1);
                    }
                    else
                        bpmpos[verticalBpmCounter] = atof(string1);
                }

                else {    /*bpm data (tbach) */
                    if (hv[bpmCounter] == 0)
                        matrix[horizontalBpmCounter][columnCounter - 3] = atof(string1);
                    else
                        matrix[verticalBpmCounter][columnCounter - 3] = atof(string1);
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

        Data.Nbpms = bpmCounter;
        Data.counth = horizontalBpmCounter + 1;
        Data.countv = verticalBpmCounter + 1;

        /* now redefine turns as the minimum of the Nturns read and the DrivingTerms data card */
        /* NB assumes all BPMs have the same number of turns as the last one read is used */
        if (Data.turns > Nturns) Data.turns = Nturns;

        /* Some statistics and checks */
        printf("Total number of pick-ups: %d Last turn number: %d, turns to run: %d\n", Data.Nbpms, Nturns, Data.turns);
        printf("Horizontal pick-ups: %d   Vertical pick-ups: %d\n", Data.counth, -MAXPICK / 2 + Data.countv);
        printf("name of BPM[0]: %s, pos: %f, first turn: %f, second turn: %f, last turn: %f, last turn to run: %f \n",
             bpmname[0], bpmpos[0], matrix[0][0], matrix[0][1], matrix[0][Nturns - 1], matrix[0][Data.turns - 1]);
        /* end of data file reading */

        printf("kick: %d \n", Data.kick);
        /* searching for two working adjacent pick-ups */
        /* after the Q-kickers for finding the kick */
        if (Data.kick < 0) {
            Data.start = -(Data.kcase - 1) * MAXPICK / 2 + 2;
            while (label[Data.start] == 0 || label[Data.start + 2] == 0) {
                Data.start = Data.start + 2;
            }

            printf("looking for kick in pick-up:%d\n", Data.start + 1);
            /* Find kick here and get kick */
            for (columnCounter = 1; (Data.kick < 0) && (columnCounter < Data.turns); ++columnCounter) {
                if (fabs(matrix[Data.start][columnCounter] - matrix[Data.start][columnCounter - 1]) > Data.kper) {
                    Data.kick = columnCounter;
                }
            }

            if (Data.kick < 0) {
                fprintf(stderr, "NO KICK FOUND\n");
                exit(EXIT_FAILURE);
            } else
                printf("Found kick in turn:%d\n", Data.kick + 1);    /*Natural count */
        }

        if (Data.kick > 0) {
            for (i = 0; i < MAXPICK; i++) {
                if (label[i] == 1) {
                    for (j = Data.kick; j < Data.turns; j++)
                        matrix[i][j - Data.kick] = matrix[i][j];
                }
            }
            Data.turns -= Data.kick;
        }
        printf("Turns to be processed after kick offset: %d matrix[0][0]: %f \n", Data.turns, matrix[0][0]);

        /* First part of the analysis: Determine  phase of all pick-ups and noise */
        writeSussixInput(sussixInputFilePath, Data.turns, Data.istun, Data.tunex, Data.tuney);

        if (Data.counth >= (Data.countv - MAXPICK / 2))
            Data.maxcounthv = Data.counth;
        else
            Data.maxcounthv = -MAXPICK / 2 + Data.countv;

        if (Data.maxcounthv > Data.pickend)
            Data.maxcounthv = Data.pickend;
        /*Shouldn't this check if counth+(countv-MAXPICK/2)>MAXPICK?*/
        if (Data.maxcounthv >= MAXPICK) {
            fprintf(stderr, "\nNot enough Pick-up mexmory\n");
            exit(EXIT_FAILURE);
        }
        printf("BPMs in loop: %d, pickstart: %d, resulting loop length: %d\n",
             Data.maxcounthv, Data.pickstart, Data.maxcounthv - Data.pickstart);

#pragma omp parallel for private(i, horizontalBpmCounter, verticalBpmCounter, j, maxamp, calculatednattunex, calculatednattuney, calculatednatampx, calculatednatampy)
        for (i = Data.pickstart; i < Data.maxcounthv; ++i) {
            horizontalBpmCounter = i;
            verticalBpmCounter = i + MAXPICK / 2;

            if (verticalBpmCounter >= Data.countv)
                verticalBpmCounter = Data.countv - 1;
            if (horizontalBpmCounter >= Data.counth)
                horizontalBpmCounter = Data.counth - 1;
            if (horizontalBpmCounter < 0 || verticalBpmCounter < 0)
            {
                fprintf(stderr, "horizontal or vertical BpmCounter < 0. Should not happen.\n");
                exit(EXIT_FAILURE);
            }
            /*printf("BPM indexes (H,V):%d %d\n", horizontalBpmCounter, verticalBpmCounter); This is not synchronised and can produce random ordered output for multiple threads (tbach) 
            Commented out because it provides no information and makes output even messier (asherman)*/

            for (j = 0; j < MAXTURNS; ++j) {
                doubleToSend[j] = matrix[horizontalBpmCounter][j];
                doubleToSend[j + MAXTURNS] = matrix[verticalBpmCounter][j];
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
            if (nattunex > NATTUNE_DEFAULT) {
                for (j = 0; j < 300; ++j) {
                    if ((nattunex - Data.istun < allfreqsx[j] && allfreqsx[j] < nattunex + Data.istun) && (maxamp < allampsx[j])) {
                        maxamp = allampsx[j];
                        calculatednattunex = allfreqsx[j];
                        calculatednatampx = maxamp;
                    }
                }
            }
            maxamp = 0;
            calculatednattuney = NATTUNE_DEFAULT;
            if (nattuney > NATTUNE_DEFAULT) {
                for (j = 0; j < 300; ++j) {
                    if ((nattuney - Data.istun < allfreqsy[j] && allfreqsy[j] < nattuney + Data.istun) && (maxamp < allampsy[j])) {
                        maxamp = allampsy[j];
                        calculatednattuney = allfreqsy[j];
                        calculatednatampy = maxamp;
                    }
                }
            }

            #pragma omp critical
            {
                label[horizontalBpmCounter] = BPMstatus(1, Data.turns);
                if (labelrun == 1)
                    noiseFile << std::scientific << "1 " << horizontalBpmCounter << "  " <<  noise1 << ' ' <<  noiseAve << ' ' << maxpeak << ' ' << maxfreq << ' ' << maxmin << ' ' << nslines << ' ' << label[i] << ' ' << phase[0] / 360. << std::endl;
                          
                /* PRINT LINEAR FILE */
                if (amplitude[0] > 0 && label[i] == 1 && horizontalBpmCounter == i) {
                    linxFile <<  std::scientific << '"' << bpmname[horizontalBpmCounter] << "\" " << bpmpos[horizontalBpmCounter] << ' ' << horizontalBpmCounter << ' ' << label[horizontalBpmCounter] << ' ' << tune[0] << ' ' <<
                            phase[0] / 360. << ' ' << amplitude[0] << ' ' << noise1 << ' ' << maxmin << ' ' << amplitude[2] / amplitude[0] << ' ' << phase[2] / 360. << ' ' << co << ' ' << co2 << ' ' << amplitude[1] / amplitude[0] << ' ' <<
                            phase[1] / 360. << ' ' << amplitude[12] / amplitude[0] << ' ' << phase[12] / 360. << ' ' << amplitude[6] / amplitude[0] << ' ' << 
                            phase[6] / 360. << ' ' << amplitude[14] / amplitude[0]  << ' ' << phase[14] / 360. << ' ' << amplitude[16] / amplitude[0] << ' ' <<
                            phase[16] / 360. << ' ' << amplitude[18] / amplitude[0] << ' ' << phase[18] / 360. << ' ' << calculatednattunex << ' ' << calculatednatampx << std::endl;
                  
                    ++Data.tunecountx;
                    Data.tunesumx += tune[0];
                    Data.tune2sumx += tune[0] * tune[0];
                    if (calculatednattunex > NATTUNE_DEFAULT) { /*  Initialized to -100. Condition true if nat tune found */
                        ++Data.nattunexcount;
                        Data.nattunexsum += calculatednattunex;
                        Data.nattunex2sum += calculatednattunex * calculatednattunex;
                    }

                    /* Horizontal Spectrum output */
                    if (i < 10) {
                        spectrumFilePath = workingDirectoryPath+'/'+bpmname[i]+".x";
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
                label[verticalBpmCounter] = BPMstatus(2, Data.turns);
                if (labelrun == 1)
                    noiseFile << std::scientific << "2 " << verticalBpmCounter << "  " <<  noise1 << ' ' <<  noiseAve << ' ' << maxpeak << ' ' << maxfreq << ' ' << maxmin << ' ' << nslines << ' ' << label[verticalBpmCounter] << ' ' << phase[3] / 360. << std::endl;
                if (amplitude[3] > 0 && label[verticalBpmCounter] == 1 && verticalBpmCounter == i + MAXPICK / 2) {
                    linyFile <<  std::scientific << '"' << bpmname[verticalBpmCounter] << "\" " << bpmpos[verticalBpmCounter] << ' ' << verticalBpmCounter << ' ' << label[verticalBpmCounter] << ' ' << tune[1] << ' ' <<
                            phase[3] / 360. << ' ' << amplitude[3] << ' ' << noise1 << ' ' << maxmin << ' ' << amplitude[5] / amplitude[3] << ' ' << phase[5] / 360. << ' ' << co << ' ' << co2 << ' ' <<
                            amplitude[13] / amplitude[3] << ' ' << phase[13] / 360. << ' ' << amplitude[15] / amplitude[3] << ' ' << phase[15] / 360. << ' ' <<
                            amplitude[17] / amplitude[3] << ' ' << phase[17] / 360. << ' ' << amplitude[4] / amplitude[3] << ' ' << phase[4] / 360. << ' ' <<
                            amplitude[11] / amplitude[3] << ' ' << phase[11] / 360. << ' ' << calculatednattuney << ' ' << calculatednatampy << std::endl;   
                    ++Data.tunecounty;
                    Data.tunesumy += tune[1];
                    Data.tune2sumy += tune[1] * tune[1];
                    if (calculatednattuney > NATTUNE_DEFAULT) { /*  Initialized to -100. Condition true if nat tune found */
                        ++Data.nattuneycount;
                        Data.nattuneysum += calculatednattuney;
                        Data.nattuney2sum += calculatednattuney * calculatednattuney;
                    }
                    if (verticalBpmCounter < MAXPICK / 2 + 10) {
                        spectrumFilePath = workingDirectoryPath+'/'+bpmname[i]+".y";
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
        if (labelrun == 1) noiseFile.close();
        /* Sort and move the "@..." lines to the top of the _linx/y files */
        formatLinFile(linxFilePath,
                Data.tunecountx, Data.tunesumx, Data.tune2sumx, Data.nattunexcount, Data.nattunexsum, Data.nattunex2sum, 1);
        formatLinFile(linyFilePath,
                Data.tunecounty, Data.tunesumy, Data.tune2sumy, Data.nattuneycount, Data.nattuneysum, Data.nattuney2sum, 2);
    } /* end of while loop over all files to analyse */
    drivingTermsFile.close();

    if (loopcounter == 0)
        std::cout << "Drivingterms file has bad input, no data ever read\n";
    
    return EXIT_SUCCESS;
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
int BPMstatus(const int plane, const int turns)
{
    double aux = 0, ave = 0, amp = 0,
        maxe = -500000.0, mine = 500000.0;
    int il,counter,counter3 = 0;

    maxpeak = 0;                /*Initialising */
    co = 0.0;
    co2 = 0.0;
    /* If peak-to-peak signal smaller than MINSIGNAL reject
     * Update: No longer, see above*/
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

    /*if (maxmin < MINSIGNAL || maxmin > MAXSIGNAL)
        return 0;*/
    
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

        if ((aux < windowa2 && aux > windowa1)
            || (aux < windowb2 && aux > windowb1)) {
            if (amp < 0) {      /* Something in sussix went wrong */
                noise1 = 100;
                noiseAve = 100;
                maxpeak = 100;
                maxfreq = 100;
                return 1;
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

    /* If tune line isn't larger than background reject
     * Update: No longer, see above */

    if ((windowa1 < maxfreq && maxfreq < windowa2)
     || (windowb1 < maxfreq && maxfreq < windowb2))
        printf("NoiseWindow includes largest lines, amp %e freq %e!!!!\n",
               maxpeak, maxfreq);

    /*if (maxpeak <= noiseAve + SIGMACUT * noise1)
        return 0;

    if (noise1 > BADPICKUP)
        return 0;*/

     /*Otherwise pick-up succeeded to first cuts*/ 
    return 1;
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

/* what is happening here? we want to sort the lin file and put 2 lines on top.
     * We could do this in plain c, but it is a bit of work, so we use shell tools.
     * (this breaks OS compatibility, but this is not a priority right now)
     * So we put calculations for tune at a temp file.
     * Then we add the first 2 lines (column headers) to temp file
     * Then we sort all files from the third line to the end, but them in the temp file (sorted)
     * Then we rename the temp file to orig file
     * --tbach*/

/*UPDATE:  Use of shell has been eliminated, uses only c allowing for OS compatibility
 * -- asherman (07/2013) */

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
