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
   */
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

int readDrivingTerms(FILE*, int*, char*, const int);
int getNextInt(FILE*);
void writeSussixInput(const char*, const int, const double, const double, const double);
int BPMstatus(const int, const int);
int containsChar(const char*, const char);
int containsString(const char*, const char*);
int canOpenFile(const char*);
void assertSmaller(const int, const int, const char*);
void setPath(char*, const size_t, const char*, const char*, const char*);
FILE* getFileToWrite(const char*);
FILE* getFileToRead(const char*);
FILE* __getFileWithMode(const char*, const char*, const char*);
void printTuneValue(FILE*, const char*, const int, const double, const double);
void formatLinFile(const char*, const int, const double, const double, const char*, const int, const double, const double, const char*);

char driveInputFilePath[2000], drivingTermsFilePath[2000], noiseFilePath[500]; /*TODO create size dynamically? (tbach)*/

double calculatednattuney, calculatednattunex, co, co2, noiseAve, maxamp,
        maxfreq, maxmin, maxpeak, nattunex, nattuney, noise1,
        windowa1, windowa2, windowb1, windowb2;
double allampsx[300], allampsy[300], allbpmamp[MAXPICK], allbpmphase[MAXPICK],
        allfreqsx[300], allfreqsy[300], amplitude[19], bpmpos[MAXPICK],
        doubleToSend[MAXTURNS4 + 4], matrix[MAXPICK][MAXTURNS], phase[19],
        tune[2];

int labelrun, nslines = 0, Nturns = 0;
int label[MAXPICK], hv[MAXPICK], hvt[MAXPICK];

#pragma omp threadprivate(amplitude, doubleToSend, tune, phase,\
        noise1, noiseAve, maxpeak, maxfreq, maxmin, co, co2,\
        allfreqsx, allampsx, allfreqsy, allampsy)

int main(int argc, char **argv)
{
    double istun, kper, tunex, tuney, nattunexsum, nattunex2sum, nattuneysum,
            nattuney2sum, tunesumx, tunesumy, tune2sumx, tune2sumy;

    int i, bpmCounter, columnCounter, counth, countv, flag,
            horizontalBpmCounter, j, kcase, kick, kk, maxcounthv, Nbpms,
            pickstart, pickend, start, turns, verticalBpmCounter, 
            tunecountx, tunecounty, nattunexcount, nattuneycount;

    size_t charCounter;

    char bpmFileName[300], dataFilePath[500], linxFilePath[2000],
            linyFilePath[2000], spectrumFilePath[400], string1[1000],
            string2[300], sussixInputFilePath[4000], workingDirectoryPath[4000];

    char* lastSlashIndex;
    char* bpmname[MAXPICK];

    FILE *dataFile, *linxFile, *linyFile, *noiseFile, *spectrumFile, *driveInputFile,
            *drivingTermsFile;

    omp_set_dynamic(0);
    /* Memory allocation */
    for (i = 0; i < MAXPICK; i++)
        bpmname[i] = (char *) calloc(50, sizeof(char));

    /*  Path to DrivingTerms and Drive.inp */
    setPath(workingDirectoryPath, sizeof(workingDirectoryPath), argv[1], "", "");
    if (!canOpenFile(workingDirectoryPath)) {
        printf("Directory is not readable: %s\n", workingDirectoryPath);
        exit(EXIT_FAILURE);
    }
    else
        printf("\nWorking directory: %s\n", workingDirectoryPath);

    setPath(drivingTermsFilePath, sizeof(drivingTermsFilePath), workingDirectoryPath, "", "DrivingTerms");
    setPath(driveInputFilePath, sizeof(driveInputFilePath), workingDirectoryPath, "", "Drive.inp");
    setPath(sussixInputFilePath, sizeof(sussixInputFilePath), workingDirectoryPath, "", "sussix_v4.inp");



    /* check the file drivingTermsFilePath */
    if (!canOpenFile(drivingTermsFilePath)) {
        printf("\nNo file %s for reading the name of the Data file\n", drivingTermsFilePath);
        exit(EXIT_FAILURE);
    }

    /* check the input file Drive.inp */
    if (!canOpenFile(driveInputFilePath)) {
        printf("\nNo input file %s\n", driveInputFilePath);
        exit(EXIT_FAILURE);
    }


    /* set all options to defaults, it could happen that they are not included in the file (tbach) */
    kick = kcase = pickstart = pickend = labelrun = 0;
    kper = tunex = tuney = istun = windowa1 = windowa2 = windowb1 = windowb2 = 0.0;
    nattunex = nattuney = NATTUNE_DEFAULT;

    /* input/option file reading start */
    driveInputFile = getFileToRead(driveInputFilePath);
    while (1) {
        /* string1 will be the option name, s the value. expected separator: '=' (tbach) */
        for (charCounter = 0; (charCounter < sizeof(string1)) && ((string1[charCounter] = (char)getc(driveInputFile)) != '=') &&
            (string1[charCounter] != '\n') && (string1[charCounter] != EOF); ++charCounter) ;
        if (charCounter >= sizeof(string1))
        {
            string1[charCounter - 1] = '\0';
            printf("Option name longer than sizeof(ss): %u, read: %s", sizeof(string1), string1);
            exit(EXIT_FAILURE);
        }
        if ((string1[charCounter] == '\n') || (string1[charCounter] == EOF)) /* we expect to have one '=' per line (tbach) */
            break;
        string1[charCounter] = '\0'; /* '=' is replaced by string termination, this is ok (tbach) */


        for (charCounter = 0; ((charCounter < sizeof(string2)) && (string2[charCounter] = (char)getc(driveInputFile)) != EOF) &&
            (string2[charCounter] != '\n');)
        {
            /* if we have ' ' or '\t', do not increase counter and just overwrite them (tbach) */
            if (string2[charCounter] == ' ' || string2[charCounter] == '\t')
                printf("Found trailing(?) whitespace or tab for line with: %s (It is ignored, but should be removed)\n", string1);
            else
                ++charCounter;
        }
        if (charCounter >= sizeof(string2))
        {
            string2[charCounter - 1] = '\0';
            printf("Option value longer than sizeof(string2): %u, read: %s", sizeof(string2), string2);
            exit(EXIT_FAILURE);
        }
        string2[charCounter] = '\0'; /* '\n' or 'EOF' is replaced by string termination, this is ok (tbach) */

        if (containsString(string1, "KICK") && strlen(string1) == 4) kick = atoi(string2) - 1;     /* C arrays start at 0 */
        if (containsString(string1, "CASE")) kcase = atoi(string2);
        if (containsString(string1, "KPER")) kper = atof(string2);
        if (containsString(string1, "TUNE X")) tunex = atof(string2);
        if (containsString(string1, "TUNE Y")) tuney = atof(string2);
        if (containsString(string1, "PICKUP START")) pickstart = atoi(string2);
        if (containsString(string1, "PICKUP END")) pickend = atoi(string2);
        if (containsString(string1, "ISTUN")) istun = atof(string2);
        if (containsString(string1, "LABEL")) labelrun = atoi(string2);
        if (containsString(string1, "WINDOWa1")) windowa1 = atof(string2);
        if (containsString(string1, "WINDOWa2")) windowa2 = atof(string2);
        if (containsString(string1, "WINDOWb1")) windowb1 = atof(string2);
        if (containsString(string1, "WINDOWb2")) windowb2 = atof(string2);
        if (containsString(string1, "NATURAL X")) nattunex = atof(string2);
        if (containsString(string1, "NATURAL Y")) nattuney = atof(string2);
    }
    fclose(driveInputFile);
    /* input/option file reading end */

    if (kick >= 0)
        printf("Known kick in turn %d\n", kick + 1);
    if (kcase == 1)
        printf("Horizontal case\n");
    else if (kcase == 0)
        printf("Vertical case\n");
    else {
        printf("No proper kcase in Drive.inp\n");
        exit(EXIT_FAILURE);
    }

    if (labelrun == 1)
        printf("\n LABELRUN: NOISE FILES WILL BE WRITTEN TO NOISEPATH\n");
    printf("pickstart: %d, pickend: %d\n", pickstart, pickend);
    if (pickstart < 0 || pickstart > pickend || pickstart > MAXPICK)
    {
        printf("Bad value for pickstart. Must be >= 0 and < pickend and <= MAXPICK(=%d)", MAXPICK);
        exit(EXIT_FAILURE);
    }


    drivingTermsFile = getFileToRead(drivingTermsFilePath);
    turns = 0;
    /* From drivingTermsFilePath assign dataFilePath, assign turns. */
    while (readDrivingTerms(drivingTermsFile, &turns, dataFilePath, sizeof(dataFilePath))) {
        /* set all values to be calculated to default values */
        tunecountx = tunecounty = nattunexcount = nattuneycount = 0;
        tunesumx = tunesumy = tune2sumx = tune2sumy = nattunexsum = nattunex2sum= nattuneysum = nattuney2sum = 0.0;

        /* Check the file dataFilePath */
        if (!canOpenFile(dataFilePath)) {
            /* doesn't exist --> try next one */
            printf("\nCan not open data file: %s\n", dataFilePath);
            continue;
        }
        printf("Data file: %s\n", dataFilePath);

        /*constructing name of output files */
        lastSlashIndex = dataFilePath;
        if (containsChar(dataFilePath, '/'))
            lastSlashIndex = strrchr(dataFilePath, '/') + 1; /* search last occurrence, subtract pointer. we search 2 times here, who cares (tbach) */

        /* copy everything from behind the last slash until the end to bpmFileName (tbach) */
        assertSmaller(strlen(lastSlashIndex), sizeof(bpmFileName), "set bpmfile");
        strncpy(bpmFileName, lastSlashIndex, sizeof(bpmFileName)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
        printf("bpmFileName: %s\n", bpmFileName);

        setPath(noiseFilePath, sizeof(noiseFilePath), workingDirectoryPath, bpmFileName, "_noise");
        setPath(linxFilePath, sizeof(linxFilePath), workingDirectoryPath, bpmFileName, "_linx");
        setPath(linyFilePath, sizeof(linyFilePath), workingDirectoryPath, bpmFileName, "_liny");

        assertSmaller(strlen(bpmFileName) + strlen("_bpm"), sizeof(bpmFileName), "modify bpmFileName");
        strncat(bpmFileName, "_bpm", sizeof(bpmFileName) - strlen(bpmFileName) - 1);


        linxFile = getFileToWrite(linxFilePath);
        linyFile = getFileToWrite(linyFilePath);
        fprintf(linxFile,
                "* NAME S    BINDEX SLABEL TUNEX MUX  AMPX NOISE PK2PK AMP01 PHASE01 CO   CORMS AMP_20 PHASE_20 AMP02 PHASE02 AMP_30 PHASE_30 AMP_1_1 PHASE_1_1 AMP2_2 PHASE2_2 AMP0_2 PHASE0_2 NATTUNEX\n");
        fprintf(linxFile,
                "$ %%s  %%le %%le   %%le   %%le  %%le %%le %%le  %%le  %%le  %%le    %%le %%le  %%le   %%le     %%le  %%le    %%le   %%le     %%le    %%le      %%le   %%le     %%le   %%le     %%le\n");
        fprintf(linyFile,
                "* NAME S    BINDEX SLABEL TUNEY MUY  AMPY NOISE PK2PK AMP10 PHASE10 CO   CORMS AMP_1_1 PHASE_1_1 AMP_20 PHASE_20 AMP1_1 PHASE1_1 AMP0_2 PHASE0_2 AMP0_3 PHASE0_3 NATTUNEY\n");
        fprintf(linyFile,
                "$ %%s  %%le %%le   %%le   %%le  %%le %%le %%le  %%le  %%le  %%le    %%le %%le  %%le    %%le      %%le   %%le     %%le   %%le     %%le   %%le     %%le   %%le     %%le\n");

        if (labelrun == 1) noiseFile = getFileToWrite(noiseFilePath);


        flag = 0;
        for (i = 0; i < MAXPICK; i++)
            label[i] = 0;

        /* start data file reading, constructing a matrix with all the data from the pick-ups */
        bpmCounter = 0;
        columnCounter = 0;
        horizontalBpmCounter = -1;
        verticalBpmCounter = MAXPICK / 2 - 1;
        i = 0;
        dataFile = getFileToRead(dataFilePath);
        string1[0] = (char)getc(dataFile);
        while (string1[0] == '#') {       /* then it is a comment line (tbach) */
            while (getc(dataFile) != '\n');       /* read until the end of the line (tbach) */
            string1[0] = (char)getc(dataFile);        /* read the first char of the new line (tbach) */
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
                    string1[i] = (char)getc(dataFile);
                    if (i > 100) {
                        string1[i + 1] = '\0';
                        printf("Found a value which has more than 100 characters, exit parsing."
                            "This is most probably a malformatted file. bpmCounter=%d columnCounter=%d string1=%s\n", bpmCounter, columnCounter, string1);
                        exit(EXIT_FAILURE);
                    }
                }
                string1[i + 1] = string1[i];
                string1[i] = '\0';
                if (LOG_INFO)
                    printf("%s ", string1);
                if (columnCounter >= MAXTURNS) {
                    printf("Found >= %d Turns, this turn size not supported. Reduce amount of turns. bpmCounter:%d", MAXTURNS - 3, bpmCounter); /* 0,1,2 is plane, name and location (tbach) */
                    exit(EXIT_FAILURE);
                }
                if (bpmCounter >= MAXPICK) {
                    printf("Found >= %d BPMs, this size is not supported. Reduce amount of BPMs. columnCounter:%d", MAXPICK, columnCounter);
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
                            printf("horizontalBpmCounter < 0. Should not happen. Probably malformatted input file?");
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
                            printf("horizontalBpmCounter < 0. Should not happen. Probably malformatted input file?");
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
                string1[0] = (char)getc(dataFile);
        }
        fclose(dataFile);

        Nbpms = bpmCounter;
        counth = horizontalBpmCounter + 1;
        countv = verticalBpmCounter + 1;

        /* now redefine turns as the minimum of the Nturns read and the DrivingTerms data card */
        /* NB assumes all BPMs have the same number of turns as the last one read is used */
        if (turns > Nturns) turns = Nturns;

        /* Some statistics and checks */
        printf("Total number of pick-ups: %d Last turn number: %d, turns to run: %d\n", Nbpms, Nturns, turns);
        printf("Horizontal pick-ups: %d   Vertical pick-ups: %d\n", counth, -MAXPICK / 2 + countv);
        printf("name of BPM[0]: %s, pos: %f, first turn: %f, second turn: %f, last turn: %f, last turn to run: %f \n",
             bpmname[0], bpmpos[0], matrix[0][0], matrix[0][1], matrix[0][Nturns - 1], matrix[0][turns - 1]);
        /* end of data file reading */

        printf("kick: %d \n", kick);
        /* searching for two working adjacent pick-ups */
        /* after the Q-kickers for finding the kick */
        if (kick < 0) {
            start = -(kcase - 1) * MAXPICK / 2 + 2;
            while (label[start] == 0 || label[start + 2] == 0) {
                start = start + 2;
            }

            printf("looking for kick in pick-up:%d\n", start + 1);
            /* Find kick here and get kick */
            for (bpmCounter = 1; (kick < 0) && (bpmCounter < turns); ++bpmCounter) {
                if (fabs(matrix[start][bpmCounter] - matrix[start][bpmCounter - 1]) > kper) {
                    kick = bpmCounter;
                    break;
                }
            }

            if (kick < 0) {
                printf("NO KICK FOUND\n");
                exit(EXIT_FAILURE);
            } else
                printf("Found kick in turn:%d\n", kick + 1);    /*Natural count */
        }

        if (kick > 0) {
            for (i = 0; i < MAXPICK; i++) {
                if (label[i] == 1) {
                    for (j = kick; j < turns; j++)
                        matrix[i][j - kick] = matrix[i][j];
                }
            }
            turns -= kick;
        }
        printf("Turns to be processed after kick offset: %d matrix[0][0]: %f \n", turns, matrix[0][0]);

        /* First part of the analysis: Determine  phase of all pick-ups and noise */
        writeSussixInput(sussixInputFilePath, turns, istun, tunex, tuney);

        if (counth >= (countv - MAXPICK / 2))
            maxcounthv = counth;
        else
            maxcounthv = -MAXPICK / 2 + countv;

        if (maxcounthv > pickend)
            maxcounthv = pickend;
        if (maxcounthv >= MAXPICK) {
            printf("\nNot enough Pick-up mexmory\n");
            exit(EXIT_FAILURE);
        }
        printf("BPMs in loop: %d, pickstart: %d, resulting loop length: %d\n",
             maxcounthv, pickstart, maxcounthv - pickstart);

        #pragma omp parallel for private(i, horizontalBpmCounter, verticalBpmCounter, kk, maxamp, calculatednattunex, calculatednattuney)
        for (i = pickstart; i < maxcounthv; ++i) {
            horizontalBpmCounter = i;
            verticalBpmCounter = i + MAXPICK / 2;

            if (verticalBpmCounter >= countv)
                verticalBpmCounter = countv - 1;
            if (horizontalBpmCounter >= counth)
                horizontalBpmCounter = counth - 1;
            if (horizontalBpmCounter < 0 || verticalBpmCounter < 0)
            {
                printf("horizontal or vertical BpmCounter < 0. Should not happen.");
                exit(EXIT_FAILURE);
            }
            printf("BPM indexes (H,V):%d %d\n", horizontalBpmCounter, verticalBpmCounter); /* This is not synchronised and can produce random ordered output for multiple threads (tbach) */

            for (kk = 0; kk < MAXTURNS; ++kk) {
                doubleToSend[kk] = matrix[horizontalBpmCounter][kk];
                doubleToSend[kk + MAXTURNS] = matrix[verticalBpmCounter][kk];
                doubleToSend[kk + 2 * MAXTURNS] = 0.0;
                doubleToSend[kk + 3 * MAXTURNS] = 0.0;
            }

            /* This calls the external Fortran code (tbach) */
            sussix4drivenoise_(&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0], sussixInputFilePath);

            /* Let's look for natural tunes in the istun range if natural tunes input is given*/
            maxamp = 0;
            calculatednattunex = NATTUNE_DEFAULT;
            if (nattunex > NATTUNE_DEFAULT) {
                for (kk = 0; kk < 300; ++kk) {
                    if ((nattunex - istun < allfreqsx[kk] && allfreqsx[kk] < nattunex + istun) && (maxamp < allampsx[kk])) {
                        maxamp = allampsx[kk];
                        calculatednattunex = allfreqsx[kk];
                    }
                }
            }
            maxamp = 0;
            calculatednattuney = NATTUNE_DEFAULT;
            if (nattuney > NATTUNE_DEFAULT) {
                for (kk = 0; kk < 300; ++kk) {
                    if ((nattuney - istun < allfreqsy[kk] && allfreqsy[kk] < nattuney + istun) && (maxamp < allampsy[kk])) {
                        maxamp = allampsy[kk];
                        calculatednattuney = allfreqsy[kk];
                    }
                }
            }

            #pragma omp critical
            {
                allbpmphase[horizontalBpmCounter] = phase[0];
                allbpmphase[verticalBpmCounter] = phase[3];
                allbpmamp[horizontalBpmCounter] = amplitude[0];
                allbpmamp[verticalBpmCounter] = amplitude[3];
                label[horizontalBpmCounter] = BPMstatus(1, turns);
                if (labelrun == 1)
                    fprintf(noiseFile, "1 %d  %e %e %e %e %e %d %d %f\n",
                            horizontalBpmCounter, noise1, noiseAve, maxpeak, maxfreq, maxmin, nslines, label[i], phase[0] / 360.);

                /* PRINT LINEAR FILE */
                if (amplitude[0] > 0 && label[i] == 1 && horizontalBpmCounter == i) {
                    fprintf(linxFile, "\"%s\" %e %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                            bpmname[horizontalBpmCounter], bpmpos[horizontalBpmCounter], horizontalBpmCounter, label[horizontalBpmCounter], tune[0],
                            phase[0] / 360., amplitude[0], noise1, maxmin, amplitude[2] / amplitude[0], phase[2] / 360.,
                            co, co2, amplitude[1] / amplitude[0],
                            phase[1] / 360., amplitude[12] / amplitude[0], phase[12] / 360., amplitude[6] / amplitude[0],
                            phase[6] / 360., amplitude[14] / amplitude[0], phase[14] / 360., amplitude[16] / amplitude[0],
                            phase[16] / 360., amplitude[18] / amplitude[0], phase[18] / 360.,  calculatednattunex);
                    ++tunecountx;
                    tunesumx += tune[0];
                    tune2sumx += tune[0] * tune[0];
                    if (calculatednattunex > NATTUNE_DEFAULT) { /*  Initialized to -100. Condition true if nat tune found */
                        ++nattunexcount;
                        nattunexsum += calculatednattunex;
                        nattunex2sum += calculatednattunex * calculatednattunex;
                    }

                    /* Horizontal Spectrum output */
                    if (i < 10) {
                        setPath(spectrumFilePath, sizeof(spectrumFilePath), workingDirectoryPath, bpmname[i], ".x");
                        spectrumFile = getFileToWrite(spectrumFilePath);
                        fprintf(spectrumFile, "%s %s %s\n", "*", "FREQ", "AMP");
                        fprintf(spectrumFile, "%s %s %s\n", "$", "%le", "%le");
                        for (kk = 0; kk < 300; ++kk)
                            fprintf(spectrumFile, "%e %e\n", allfreqsx[kk], allampsx[kk]);
                        fclose(spectrumFile);
                    }
                }
                label[verticalBpmCounter] = BPMstatus(2, turns);
                if (labelrun == 1)
                    fprintf(noiseFile, "2 %d  %e %e %e %e %e %d %d %f\n",
                            verticalBpmCounter, noise1, noiseAve, maxpeak, maxfreq, maxmin, nslines, label[verticalBpmCounter], phase[3] / 360.);
                if (amplitude[3] > 0 && label[verticalBpmCounter] == 1 && verticalBpmCounter == i + MAXPICK / 2) {
                    fprintf(linyFile, "\"%s\" %e %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                            bpmname[verticalBpmCounter], bpmpos[verticalBpmCounter], verticalBpmCounter, label[verticalBpmCounter], tune[1], phase[3] / 360., amplitude[3], noise1,
                            maxmin, amplitude[5] / amplitude[3], phase[5] / 360., co, co2,
                            amplitude[13] / amplitude[3], phase[13] / 360., amplitude[15] / amplitude[3], phase[15] / 360.,
                            amplitude[17] / amplitude[3], phase[17] / 360., amplitude[4] / amplitude[3], phase[4] / 360.,
                            amplitude[11] / amplitude[3], phase[11] / 360., calculatednattuney);
                    ++tunecounty;
                    tunesumy += tune[1];
                    tune2sumy += tune[1] * tune[1];
                    if (calculatednattuney > NATTUNE_DEFAULT) { /*  Initialized to -100. Condition true if nat tune found */
                        ++nattuneycount;
                        nattuneysum += calculatednattuney;
                        nattuney2sum += calculatednattuney * calculatednattuney;
                    }
                    if (verticalBpmCounter < MAXPICK / 2 + 10) {
                        setPath(spectrumFilePath, sizeof(spectrumFilePath), workingDirectoryPath, bpmname[verticalBpmCounter], ".y");
                        spectrumFile = getFileToWrite(spectrumFilePath);
                        fprintf(spectrumFile, "%s %s %s\n", "*", "FREQ", "AMP");
                        fprintf(spectrumFile, "%s %s %s\n", "$", "%le", "%le");
                        for (kk = 0; kk < 300; ++kk)
                            fprintf(spectrumFile, "%e %e \n", allfreqsy[kk], allampsy[kk]);
                        fclose(spectrumFile);
                    }
                }
            } /* end of omp critical section */
        } /* end of parallel for */
        fclose(linxFile);
        fclose(linyFile);
        if (labelrun == 1) fclose(noiseFile);

        /* Sort and move the "@..." lines to the top of the _linx/y files */
        formatLinFile(linxFilePath,
                tunecountx, tunesumx, tune2sumx, "@ Q1 %%le %e\n@ Q1RMS %%le %e\n",
                nattunexcount, nattunexsum, nattunex2sum, "@ NATQ1 %%le %e\n@ NATQ1RMS %%le %e\n");
        formatLinFile(linyFilePath,
                tunecounty, tunesumy, tune2sumy, "@ Q2 %%le %e\n@ Q2RMS %%le %e\n",
                nattuneycount, nattuneysum, nattuney2sum, "@ NATQ2 %%le %e\n@ NATQ2RMS %%le %e\n");
    } /* end of while loop over all files to analyse */
    fclose(drivingTermsFile);

    return EXIT_SUCCESS;
}

/* *****************  */
/*    readDrivingTerms*/
/* *****************  */
int readDrivingTerms(FILE* drivingTermsFile, int* turns, char* path, const int sizeOfPath)
{
    /* This functions reads from the given FILE one line of this format:
     * <path to datafile> <int start turn> <int end turn>
     * If a line was successfully read, return true
     * If not, return false
     * --tbach */

    int charCounter = 0;

    /* this block reads the filepath */
    while (isspace((int)(path[0] = (char)getc(drivingTermsFile)))) ;
    for (charCounter = 1; ((charCounter < sizeOfPath && (path[charCounter] = (char)getc(drivingTermsFile)) != EOF) &&
         (path[charCounter] != '\n') && (path[charCounter] != '%') &&
         (path[charCounter] != ' ')); charCounter++) ;
    if (charCounter >= sizeOfPath)
    {
        printf("Error: path longer than sizeOfPath: %d", sizeOfPath);
        exit(EXIT_FAILURE);
    }
    if (path[charCounter] == EOF) { /* we do not expect an EOF here (tbach)*/
        path[charCounter] = '\0';
        return 0;
    }
    path[charCounter] = '\0';

    /* this block reads over the first number (tbach) */
    if (getNextInt(drivingTermsFile) == EOF)
        return 0;

    /* this block reads the second number (tbach) */
    *turns = getNextInt(drivingTermsFile);
    if (*turns == EOF)
        return 0;
    return 1;
}

int getNextInt(FILE* file)
{
    /* This reads the next int from the given file.
     * - skips ' '
     * - reads everything until EOF, '\n' or ' ' into some <var>
     * - converts the found char in <var> to int
     * - returns the int or EOF if something went wrong
     * --tbach */
    int result = 0;
    char string1[1000];
    size_t charCounter = 0;
    while (isspace((int)(string1[0] = (char)getc(file)))) ;  /* skip all spaces (tbach) */
    if (string1[0] == EOF)
        return EOF;
    for (charCounter = 1; ((charCounter < sizeof(string1)) && (string1[charCounter] = (char)getc(file)) != EOF) &&
         (string1[charCounter] != '\n') && (string1[charCounter] != ' '); charCounter++) ;
    if (charCounter >= sizeof(string1))
    {
        string1[charCounter - 1] = '\0';
        printf("Error: input longer than sizeof(string1): %u, string1: %s", sizeof(string1), string1);
        exit(EXIT_FAILURE);
    }
    string1[charCounter] = '\0';
    result = atoi(string1);
    return result;
}

/* ***************** */
/*    sussix_inp     */
/* ***************** */
void writeSussixInput(const char* sussixInputFilePath, const int turns, const double istun, const double tunex, const double tuney)
{
    FILE* sussixInputFile = getFileToWrite(sussixInputFilePath);
    fprintf(sussixInputFile, "C\n");
    fprintf(sussixInputFile, "C INPUT FOR SUSSIX_V4 ---17/09/1997---\n");
    fprintf(sussixInputFile, "C DETAILS ARE IN THE MAIN PROGRAM SUSSIX_V4.F\n");
    fprintf(sussixInputFile, "C\n");
    fprintf(sussixInputFile, "\n");
    fprintf(sussixInputFile, "ISIX  = 0\n");
    fprintf(sussixInputFile, "NTOT  = 1\n");
    fprintf(sussixInputFile, "IANA  = 1\n");
    fprintf(sussixInputFile, "ICONV = 0\n");
    fprintf(sussixInputFile, "TURNS = 1 %d\n", turns);
    fprintf(sussixInputFile, "NARM  = 300\n");
    fprintf(sussixInputFile, "ISTUN = 1 %e %e\n", istun, istun);
    fprintf(sussixInputFile, "TUNES = %e %e .07\n", tunex, tuney);
    fprintf(sussixInputFile, "NSUS  = 0\n");
    fprintf(sussixInputFile, "IDAM  = 2\n");
    fprintf(sussixInputFile, "NTWIX = 1\n");
    fprintf(sussixInputFile, "IR    = 1\n");
    fprintf(sussixInputFile, "IMETH = 2\n");
    fprintf(sussixInputFile, "NRC   = 4\n");
    fprintf(sussixInputFile, "EPS   = 2D-3\n"); /* EPS is the window in the secondary lines, very imp!!! */
    fprintf(sussixInputFile, "NLINE = 0\n");
    fprintf(sussixInputFile, "L,M,K = \n");
    fprintf(sussixInputFile, "IDAMX = 1\n");
    fprintf(sussixInputFile, "NFIN  = 500\n");
    fprintf(sussixInputFile, "ISME  = 0\n");
    fprintf(sussixInputFile, "IUSME = 200\n");
    fprintf(sussixInputFile, "INV   = 0\n");
    fprintf(sussixInputFile, "IINV  = 250\n");
    fprintf(sussixInputFile, "ICF   = 0\n");
    fprintf(sussixInputFile, "IICF  = 350\n");
    fclose(sussixInputFile);
}

/************   BPMstatus *************/
/* Analyse fort.300 to detect noise   */
/**************************************/
#define MINSIGNAL 0.00001
#define SIGMACUT   1.8
#define MINIMUMNOISE 0.0
#define BADPICKUP  8.0
#define MAXSIGNAL 30000
int BPMstatus(const int plane, const int turns)
{
    double aux = 0, ave = 0, amp = 0,
        maxe = -500000.0, mine = 500000.0;
    int counter, il, counter3 = 0;

    maxpeak = 0;                /*Initialising */
    co = 0.0;
    co2 = 0.0;
    /* If peak-to-peak signal smaller than MINSIGNAL reject */
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
    if (maxmin < MINSIGNAL || maxmin > MAXSIGNAL)
        return 0;

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
                return 0;
            }

            ave = amp + ave;
            noise1 = noise1 + amp * amp;
            ++counter3;
        }

    }
    if (counter3 > 0) {
        if (counter3 > 1)
            noise1 = sqrt((noise1 / counter3 - ave * ave / (counter3 * counter3)));
        else
            noise1 = 0;
        noiseAve = ave / counter3;
    } else {
        noise1 = MINIMUMNOISE;
        noiseAve = MINIMUMNOISE;
    }
    nslines = counter3;

    /* If tune line isn't larger than background reject */

    if ((windowa1 < maxfreq && maxfreq < windowa2)
     || (windowb1 < maxfreq && maxfreq < windowb2))
        printf("NoiseWindow includes largest lines, amp %e freq %e!!!!\n",
               maxpeak, maxfreq);

    if (maxpeak <= noiseAve + SIGMACUT * noise1)
        return 0;

    if (noise1 > BADPICKUP)
        return 0;

    /* Otherwise pick-up succeeded to first cuts */
    return 1;
}

int containsChar(const char* string, const char character)
{
    return (strchr(string, character) != NULL);
}

int containsString(const char* string1, const char* string2)
{
    return (strstr(string1, string2) != NULL);
}

int canOpenFile(const char* const filename)
{
    FILE* file;
    if ((file = fopen(filename, "r")) == NULL)
        return 0;
    fclose(file);
    return 1;
}

FILE* getFileToWrite(const char* const filename)
{
    return __getFileWithMode(filename, "w", "Cannot open to write: %s\n");
}

FILE* getFileToRead(const char* const filename)
{
    return __getFileWithMode(filename, "r", "Cannot open to read: %s\n");
}

FILE* __getFileWithMode(const char* const filename, const char* const mode, const char* const errormessage)
{
    FILE* file = fopen(filename, mode);
    if (file == NULL)
    {
        printf(errormessage, filename);
        exit(EXIT_FAILURE);
    }
    return file;
}

void assertSmaller(const int a, const int b, const char* message)
{
    if (a < b)
        return;
    printf("Value1: %d is not < than Value2: %d. Message: %s\n", a, b, message);
    exit(EXIT_FAILURE);
}

void setPath(char* path, const size_t sizeOfPath, const char* workingDirectoryPath, const char* bpmfile, const char* fileEnding)
{
    /* sets path to <workingDirectoryPath>/<bpmfilefile><Ending> (tbach) */
    assertSmaller(strlen(workingDirectoryPath) + 1 + strlen(bpmfile) + strlen(fileEnding), sizeOfPath, fileEnding);
    sprintf(path, "%s/%s%s", workingDirectoryPath, bpmfile, fileEnding);
    printf("%s: %s\n", fileEnding, path);
}

char cmd[6000];
char tempFilePath[2000];
FILE* tempFile;
void formatLinFile(const char* linFilePath,
        const int tunecount, const double tunesum, const double tune2sum, const char* tuneheader,
        const int nattunecount, const double nattunesum, const double nattune2sum, const char* nattuneheader) {
    /* what is happening here? we want to sort the lin file and put 2 line son top.
     * We could do this in plain c, but it is a bit of work, so we use shell tools.
     * (this breaks OS compatibility, but this is not a priority right now)
     * So we put calculations for tune at a temp file.
     * Then we add the first 2 lines (column headers) to temp file
     * Then we sort all files from the third line to the end, but them in the temp file (sorted)
     * Then we rename the temp file to orig file
     * --tbach */
    if (tunecount > 0) {
        sprintf(tempFilePath, "%s_temp", linFilePath);
        tempFile = getFileToWrite(tempFilePath);

        printTuneValue(tempFile, tuneheader, tunecount, tunesum, tune2sum);
        if (nattunecount > 0)
            printTuneValue(tempFile, nattuneheader, nattunecount, nattunesum, nattune2sum);

        fclose(tempFile);

        /* TODO fix for cross platform (tbach) */
        sprintf(cmd, "head -2 %s >> %s_temp", linFilePath,
                linFilePath);
        system(cmd);
        sprintf(cmd, "grep \\\" %s | sort -n --key=3 >> %s_temp",
                linFilePath, linFilePath);
        system(cmd);

        remove(linFilePath);
        rename(tempFilePath, linFilePath);

        printf("%s:\ntune: \n  sum: %e, count: %d, sum2: %e \nnatural tune: \n  sum: %e, count: %d, sum2: %e\n",
                linFilePath, tunesum, tunecount, tune2sum, nattunesum,
                nattunecount, nattune2sum);
    }
}

void printTuneValue(FILE* linFile, const char* header, const int count, const double tunesum, const double tune2sum)
{
    /* calculates the average and the RMS(?) (tbach) */
    fprintf(linFile, header,
            tunesum / count,
            sqrt(tune2sum / count - (tunesum / count) * (tunesum / count)));
}
