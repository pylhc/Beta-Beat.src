/* Parallelised version of Drive_God_lin.c of 15/05/2011 using Openmp.
   NB this version uses as the number of turns the minimum of the number
   actually read or the number given in the DrivingTerms file.
   
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
#define MAXTURNS 10000   /* Critial parameter in sussixfordiveNoO.f, same parameter name !*/
#define MAXTURNS4 40000 /*Always four times  MAXTURNS*/
#define MAXRUNS 100
#define M_PI 3.14159265358979310862446895044

#if !defined(LOG_INFO)
#define LOG_INFO 0              /*set to 1 to enable log prints (tbach) */
#endif

char *get_word(FILE *);
void get_name(), sussix_inp(int, char *), betabeat(), savefort300(int);
int BPMstatus(int), norm(int, int);
int omp_get_thread_num(void);
int containsChar(const char*, const char);
int containsString(const char*, const char*);
int canOpenFile(const char*);
void assertSmaller(const int, const int, const char*);
FILE* getFileToWrite(const char*);
FILE* getFileToRead(const char*);
FILE* __getFileWithMode( const char*, const char*, const char*);

char drivingTermsFilePath[2000], dataFilePath[500], noiseFilePath[500],
    driveInputFilePath[2000], ss[1000], *bpmname[MAXPICK];

double matrix[MAXPICK][MAXTURNS], bpmpos[MAXPICK], tunex, tuney, tune[2],
    amplitude[19], phase[19], istun, noise1, noiseAve,
    maxpeak, maxfreq, maxmin, co, co2, windowa1, windowa2, windowb1,
    windowb2, allfreqsx[300], allampsx[300], allfreqsy[300], allampsy[300],
    doubleToSend[MAXTURNS4 + 4], allbpmamp[MAXPICK], allbpmphase[MAXPICK];

int turns, normalisation, scum = 0, in = 0, iir = 0,
    labelrun = 0, label[MAXPICK], NoVerticalData = 0, hv[MAXPICK],
    hvt[MAXPICK], nslines, Nturns = 0;

size_t charCounter = 0;

#pragma omp threadprivate(doubleToSend,tune,tunex,tuney,amplitude,phase,\
        istun, noise1, noiseAve,\
        maxpeak, maxfreq, maxmin, co, co2,\
        allfreqsx,allampsx,allfreqsy,allampsy)

int main(int argc, char **argv)
{
    double tunesum[2], tune2sum[2], kper;

    int i = 0, kick = 0, ij, counth = 0, countv = 0, maxcounthv,
        start, flag = 0, columnCounter = 0, j, horizontalBpmCounter = 0, verticalBpmCounter = 0,
        count0[2], pickstart = 0, pickend = 0, kcase = 0, runBetaBeat = 0, bpmCounter = 0,
        lastSlashIndex, kk, Nbpms;

    char s[300], bpmfile[300], linxFilePath[2000], workingDirectoryPath[4000], 
        linyFilePath[2000], cmd[6000], sussixInputFilePath[4000], bsfile[400];

    FILE *dataFile, *linxFile, *linyFile, *noiseFile, *bpmFile, *driveInputFile;
    omp_set_dynamic(0);
    /* Memory allocation */
    for (i = 0; i < MAXPICK; i++)
        bpmname[i] = (char *) calloc(50, sizeof(char));

/*  Path to DrivingTerms and Drive.inp */    
    assertSmaller(strlen(argv[1]), sizeof(workingDirectoryPath), "copy argv[1] to workingDirectoryPath");
    strncpy(workingDirectoryPath, argv[1], sizeof(workingDirectoryPath)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
    if (!canOpenFile(workingDirectoryPath)) {
        printf("Directory is not readable: %s\n", workingDirectoryPath);
        exit(EXIT_FAILURE);
    }
    else
        printf("\n Path to Drive input folder: %s\n", workingDirectoryPath);

    /*TODO create size dynamically? (tbach)*/
    assertSmaller(strlen(workingDirectoryPath) + strlen("/sussix_v4.inp"), sizeof(sussixInputFilePath), "modify sussixInputFilePath");
    strncpy(sussixInputFilePath, workingDirectoryPath, sizeof(sussixInputFilePath)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
    strncat(sussixInputFilePath, "/sussix_v4.inp", sizeof(sussixInputFilePath) - strlen(sussixInputFilePath) - 1);

    assertSmaller(strlen(workingDirectoryPath) + strlen("/DrivingTerms"), sizeof(drivingTermsFilePath), "modify drivingTermsFilePath");
    strncpy(drivingTermsFilePath, workingDirectoryPath, sizeof(drivingTermsFilePath)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
    strncat(drivingTermsFilePath, "/DrivingTerms", sizeof(drivingTermsFilePath) - strlen("/DrivingTerms") - 1);

    assertSmaller(strlen(workingDirectoryPath) + strlen("/Drive.inp"), sizeof(driveInputFilePath), "modify driveInputFilePath");
    strncpy(driveInputFilePath, workingDirectoryPath, sizeof(driveInputFilePath)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
    strncat(driveInputFilePath, "/Drive.inp", sizeof(driveInputFilePath) - strlen(driveInputFilePath) - 1);

    printf("\n Drive.inp: %s\n", driveInputFilePath);



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


    count0[0] = 0;
    count0[1] = 0;
    tunesum[0] = 0;
    tunesum[1] = 0;
    tune2sum[0] = 0;
    tune2sum[1] = 0;


    /* input/option file reading start */
    driveInputFile = getFileToRead(driveInputFilePath);
    while (1) {
        /* ss will be the option name, s the value. expected separator: '=' (tbach) */
        for (charCounter = 0; (charCounter < sizeof(ss)) && ((ss[charCounter] = getc(driveInputFile)) != '=') && 
            (ss[charCounter] != '\n') && (ss[charCounter] != EOF); ++charCounter) ;
        if (charCounter >= sizeof(ss))
        {
            ss[charCounter - 1] = '\0';
            printf("Option name longer than sizeof(ss): %u, read: %s", sizeof(ss), ss);
            exit(EXIT_FAILURE);
        }
        if ((ss[charCounter] == '\n') || (ss[charCounter] == EOF)) /* we expect to have one '=' per line (tbach) */
            break;
        ss[charCounter] = '\0'; /* '=' is replaced by string termination, this is ok (tbach) */
        
        
        for (charCounter = 0; ((charCounter < sizeof(s)) && (s[charCounter] = getc(driveInputFile)) != EOF) && 
            (s[charCounter] != '\n') && (s[charCounter] != ' '); ++charCounter) ;
        if (charCounter >= sizeof(s))
        {
            s[charCounter - 1] = '\0';
            printf("Option value longer than sizeof(s): %u, read: %s", sizeof(s), s);
            exit(EXIT_FAILURE);
        }
        s[charCounter] = '\0';
        
        
        if (containsString(ss, "KICK") && strlen(ss) == 4) kick = atoi(s) - 1;     /* C arrays start at 0 */
        if (containsString(ss, "CASE")) kcase = atoi(s);
        if (containsString(ss, "KPER")) kper = atof(s);
        if (containsString(ss, "TUNE X")) tunex = atof(s);
        if (containsString(ss, "TUNE Y")) tuney = atof(s);
        if (containsString(ss, "PICKUP START")) pickstart = atoi(s);
        if (containsString(ss, "PICKUP END")) pickend = atoi(s);
        if (containsString(ss, "NORMALISATION")) normalisation = atoi(s);
        if (containsString(ss, "ISTUN")) istun = atof(s);
        if (containsString(ss, "BETABEATING")) runBetaBeat = atoi(s);
        if (containsString(ss, "IR")) iir = atoi(s);
        if (containsString(ss, "LABEL")) labelrun = atoi(s);
        if (containsString(ss, "WINDOWa1")) windowa1 = atof(s);
        if (containsString(ss, "WINDOWa2")) windowa2 = atof(s);
        if (containsString(ss, "WINDOWb1")) windowb1 = atof(s);
        if (containsString(ss, "WINDOWb2")) windowb2 = atof(s);
    }
    fclose(driveInputFile);
    /* input/option file reading end */
    
    if (kick >= 0)
        printf("Known kick in %d turn\n", kick + 1);
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
    printf("Normalisation: %d\n", normalisation);
    printf("pickstart: %d, pickend: %d\n", pickstart, pickend);
    if (pickstart < 0 || pickstart > pickend || pickstart > MAXPICK)
    {
        printf("Bad value for pickstart. Must be >= 0 and < pickend and <= MAXPICK(=%d)", MAXPICK);
        exit(EXIT_FAILURE);
    }


    while (scum == 0) {
        /* From drivingTermsFilePath assign dataFilePath, asign turns. */
        get_name();

        /* Check the file dataFilePath */
        if (!canOpenFile(dataFilePath)) {
            /* doesn't exist --> stop */
            printf("\nNo Data file %s\n", dataFilePath);
            break;
        }
        printf("Data File: %s\n", dataFilePath);
        /*constructing name of BPM files */
        lastSlashIndex = 0;
        if (containsChar(dataFilePath, '/'))
            lastSlashIndex = strrchr(dataFilePath, '/') - dataFilePath; /* search last occurence, substract pointer. we search 2 times here, who cares (tbach) */

        /* copy everything from behind the last slash until the end to bpmfile (tbach) */
        assertSmaller(strlen(dataFilePath + lastSlashIndex), sizeof(bpmfile), "modify bpmfile");
        strncpy(bpmfile, dataFilePath + lastSlashIndex, sizeof(bpmfile)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */

        assertSmaller(strlen(workingDirectoryPath) + 1 + strlen(bpmfile) + strlen("_noise"), sizeof(noiseFilePath), "modify noiseFilePath");
        strncpy(noiseFilePath, workingDirectoryPath, sizeof(noiseFilePath)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
        strncat(noiseFilePath, "/", sizeof(noiseFilePath) - strlen(noiseFilePath) - 1);
        strncat(noiseFilePath, bpmfile, sizeof(noiseFilePath) - strlen(noiseFilePath) - 1);
        strncat(noiseFilePath, "_noise", sizeof(noiseFilePath) - strlen(noiseFilePath) - 1);

        assertSmaller(strlen(workingDirectoryPath) + 1 + strlen(bpmfile) + strlen("_linx"), sizeof(linxFilePath), "modify linxFilePath");
        strncpy(linxFilePath, workingDirectoryPath, sizeof(linxFilePath)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
        strncat(linxFilePath, "/", sizeof(linxFilePath) - strlen(linxFilePath) - 1);
        strncat(linxFilePath, bpmfile, sizeof(linxFilePath) - strlen(linxFilePath) - 1);
        strncat(linxFilePath, "_linx", sizeof(linxFilePath) - strlen(linxFilePath) - 1);

        assertSmaller(strlen(workingDirectoryPath) + 1 + strlen(bpmfile) + strlen("_liny"), sizeof(linxFilePath), "modify linyFilePath");
        strncpy(linyFilePath, workingDirectoryPath, sizeof(linyFilePath)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
        strncat(linyFilePath, "/", sizeof(linyFilePath) - strlen(linyFilePath) - 1);
        strncat(linyFilePath, bpmfile, sizeof(linyFilePath) - strlen(linyFilePath) - 1);
        strncat(linyFilePath, "_liny", sizeof(linyFilePath) - strlen(linyFilePath) - 1);

        assertSmaller(strlen(bpmfile) + strlen("_bpm"), sizeof(bpmfile), "modify bpmfile");
        strncat(bpmfile, "_bpm", sizeof(bpmfile) - strlen(bpmfile) - 1);
        
        
        linxFile = getFileToWrite(linxFilePath);
        linyFile = getFileToWrite(linyFilePath);
        fprintf(linxFile,
                "* NAME  S   BINDEX SLABEL  TUNEX   MUX  AMPX   NOISE   PK2PK  AMP01 PHASE01 CO CORMS AMP_20  PHASE_20  AMP02 PHASE02 AMP_30  PHASE_30  AMP_1_1  PHASE_1_1  AMP2_2  PHASE2_2  AMP0_2 PHASE0_2\n");
        fprintf(linxFile,
                "$ %%s  %%le   %%le    %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le    %%le     %%le   %%le     %%le   %%le     %%le     %%le     %%le     %%le     %%le     %%le     %%le ");
        fprintf(linyFile,
                "* NAME  S   BINDEX SLABEL  TUNEY  MUY  AMPY   NOISE   PK2PK AMP10 PHASE10 CO CORMS AMP_1_1 PHASE_1_1  AMP0_1  PHASE0_1  AMP1_1 PHASE1_1  AMP0_2  PHASE0_2  AMP0_3  PHASE0_3\n");
        fprintf(linyFile,
                "$ %%s  %%le  %%le    %%le  %%le  %%le  %%le  %%le  %%le %%le  %%le %%le  %%le %%le %%le %%le %%le %%le %%le %%le %%le %%le %%le ");
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
        s[0] = getc(dataFile);
        while (s[0] == '#') {       /* then it is a comment line (tbach) */
            while (getc(dataFile) != '\n');       /* read until the end of the line (tbach) */
            s[0] = getc(dataFile);        /* read the first char of the new line (tbach) */
        }
        /* after this, we have skipped all the comment lines, and s[0] is the first character of a new line which is not a "#" (tbach) */
        if (LOG_INFO)
            printf("BPM file content:\n");
        while (s[0] != EOF) {
            if (s[0] == '\n') {
                ++bpmCounter;
                if (LOG_INFO)
                    printf("\n");
                columnCounter = 0;
            }
            if (isspace(s[0]) && flag == 1)
                flag = 0;
            if (!isspace(s[0]) && flag == 0) {
                while (!isspace(s[i]) && s[i] != EOF) {
                    ++i;
                    s[i] = getc(dataFile);
                    if (i > 100) {
                        s[i + 1] = '\0';
                        printf("Found a value which has more than 100 characters, exit parsing."
                            "This is most probably a malformatted file. bpmCounter=%d columnCounter=%d s=%s\n", bpmCounter, columnCounter, s);
                        exit(EXIT_FAILURE);
                    }
                }
                s[i + 1] = s[i];
                s[i] = '\0';
                if (LOG_INFO)
                    printf("%s ", s);
                if (columnCounter >= MAXTURNS) {
                    printf("Found >= %d Turns, this turn size not supported. Reduce amount of turns. bpmCounter:%d", MAXTURNS - 3, bpmCounter); /* 0,1,2 is plane, name and location (tbach) */
                    exit(EXIT_FAILURE);
                }
                if (bpmCounter >= MAXPICK) {
                    printf("Found >= %d BPMs, this size is not supported. Reduce amount of BPMs. columnCounter:%d", MAXPICK, columnCounter);
                    exit(EXIT_FAILURE);
                }
                if (columnCounter == 0) {   /*plane (tbach) */
                    hv[bpmCounter] = atoi(s);
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
                        strcpy(bpmname[horizontalBpmCounter], s);
                        label[horizontalBpmCounter] = 1;
                    } else {
                        hvt[verticalBpmCounter] = 1;
                        strcpy(bpmname[verticalBpmCounter], s);
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
                        bpmpos[horizontalBpmCounter] = atof(s);
                    }
                    else
                        bpmpos[verticalBpmCounter] = atof(s);
                }

                else {    /*bpm data (tbach) */
                    if (hv[bpmCounter] == 0)
                        matrix[horizontalBpmCounter][columnCounter - 3] = atof(s);
                    else
                        matrix[verticalBpmCounter][columnCounter - 3] = atof(s);
                    Nturns = columnCounter - 3 + 1;
                    /* If the last line is an empty line, then we can get the number of turns only from here.
                       First 3 are plane, name and location.
                       Plus 1 for index start at 0
                       (tbach) */
                }
                ++columnCounter;
                flag = 1;
                s[0] = s[i + 1];
                i = 0;
            }
            if (flag == 0)
                s[0] = getc(dataFile);
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
        printf("name of pick-up 0: %s, pos: %f, first turn: %f, second turn: %f, last turn overall: %f, last turn to run: %f \n",
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
            bpmCounter = 1;
            while (kick < 0 && bpmCounter < turns) {
                if (fabs(matrix[start][bpmCounter] - matrix[start][bpmCounter - 1]) > kper) {
                    kick = bpmCounter;
                    break;
                }
                ++bpmCounter;
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
        printf("Turns to be procesed after kick offset: %d matrix[0][0]: %f \n", turns, matrix[0][0]);

        /* First part of the analysis: Determine  phase of all pick-ups and noise */
        sussix_inp(1, sussixInputFilePath);

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
        #pragma omp parallel for private(i,columnCounter,ij,kk)
        for (i = pickstart; i < maxcounthv; ++i) {
            columnCounter = i;
            ij = i + MAXPICK / 2;

            if (ij >= countv)
                ij = countv - 1;
            if (columnCounter >= counth)
                columnCounter = counth - 1;
            if (columnCounter < 0)
            {
                printf("Columncounter < 0. Should not happen.");
                exit(EXIT_FAILURE);
            }
            printf("BPM indexes (H,V):%d %d\n", columnCounter, ij); /* This is not synchronized and can produce random ordered output for multiple threads (tbach) */

            for (kk = 0; kk < MAXTURNS; ++kk) {
                doubleToSend[kk] = matrix[columnCounter][kk];
                doubleToSend[kk + MAXTURNS] = matrix[ij][kk];   /* BUG to solve TUNES TOO CLOSE  */
                doubleToSend[kk + 2 * MAXTURNS] = 0.0;
                doubleToSend[kk + 3 * MAXTURNS] = 0.0;
            }

            sussix4drivenoise_(&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0], sussixInputFilePath);
            /* This calls the external fortan code (tbach) */

            #pragma omp critical
            {
                allbpmphase[columnCounter] = phase[0];
                allbpmphase[ij] = phase[3];
                allbpmamp[columnCounter] = amplitude[0];
                allbpmamp[ij] = amplitude[3];
                label[columnCounter] = BPMstatus(1);
                if (labelrun == 1)
                    fprintf(noiseFile, "1 %d  %e %e %e %e %e %d %d %f\n",
                            columnCounter, noise1, noiseAve, maxpeak, maxfreq, maxmin, nslines, label[i], phase[0] / 360.);

                /* PRINT LINEAR FILE */
                if (amplitude[0] > 0 && label[i] == 1 && columnCounter == i) {
                    fprintf(linxFile, "\n\"%s\" %e %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ",
                            bpmname[columnCounter], bpmpos[columnCounter],
                            columnCounter, label[columnCounter], tune[0],
                            phase[0] / 360., amplitude[0], noise1, maxmin,
                            amplitude[2] / amplitude[0], phase[2] / 360.,
                            co, co2, amplitude[1] / amplitude[0],
                            phase[1] / 360., amplitude[12] / amplitude[0],
                            phase[12] / 360., amplitude[6] / amplitude[0],
                            phase[6] / 360., amplitude[14] / amplitude[0],
                            phase[14] / 360., amplitude[16] / amplitude[0],
                            phase[16] / 360., amplitude[18] / amplitude[0],
                            phase[18] / 360.);
                    ++count0[0];
                    tunesum[0] += tune[0];
                    tune2sum[0] += tune[0] * tune[0];
                    /* Horizontal Spectrum output */
                    if (i < 10) {
                        sprintf(bsfile, "%s/%s.x", workingDirectoryPath, bpmname[i]);
                        bpmFile = getFileToWrite(bsfile);
                        fprintf(bpmFile, "%s %s %s\n", "*", "FREQ", "AMP");
                        fprintf(bpmFile, "%s %s %s\n", "$", "%le", "%le");
                        for (kk = 0; kk < 300; ++kk)
                            fprintf(bpmFile, "%e %e\n", allfreqsx[kk], allampsx[kk]);
                        fclose(bpmFile);
                    }
                }
                label[ij] = BPMstatus(2);
                if (labelrun == 1)
                    fprintf(noiseFile, "2 %d  %e %e %e %e %e %d %d %f\n",
                            ij, noise1, noiseAve, maxpeak, maxfreq, maxmin, nslines, label[ij], phase[3] / 360.);
                if (amplitude[3] > 0 && label[ij] == 1 && ij == i + MAXPICK / 2) {
                    fprintf(linyFile, "\n\"%s\" %e %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
                            bpmname[ij], bpmpos[ij], ij, label[ij],
                            tune[1], phase[3] / 360., amplitude[3], noise1,
                            maxmin, amplitude[5] / amplitude[3],
                            phase[5] / 360., co, co2,
                            amplitude[13] / amplitude[3], phase[13] / 360.,
                            amplitude[15] / amplitude[3], phase[15] / 360.,
                            amplitude[17] / amplitude[3], phase[17] / 360.,
                            amplitude[4] / amplitude[3], phase[4] / 360.,
                            amplitude[11] / amplitude[3], phase[11] / 360.);
                    ++count0[1];
                    tunesum[1] += tune[1];
                    tune2sum[1] += tune[1] * tune[1];
                    if (ij < MAXPICK / 2 + 10) {
                        sprintf(bsfile, "%s/%s.y", workingDirectoryPath, bpmname[ij]);
                        bpmFile = getFileToWrite(bsfile);
                        fprintf(bpmFile, "%s %s %s\n", "*", "FREQ", "AMP");
                        fprintf(bpmFile, "%s %s %s\n", "$", "%le", "%le");
                        for (kk = 0; kk < 300; ++kk)
                            fprintf(bpmFile, "%e %e \n", allfreqsy[kk], allampsy[kk]);
                        fclose(bpmFile);
                    }
                }
            }                   /* end of omp critical section */
        }
        fprintf(linxFile, "\n");
        fprintf(linyFile, "\n");
        fclose(linxFile);
        fclose(linyFile);
        /* What follows is some dirty operations to move the @ Q1 ... to the top of _linx _liny */
        if (count0[0] > 0) {
            sprintf(cmd, "head -2 %s > %s/tx", linxFilePath, workingDirectoryPath);
            system(cmd);
            sprintf(cmd, "grep \\\" %s | sort -n --key=3 >> %s/tx", linxFilePath, workingDirectoryPath);
            system(cmd);
            sprintf(cmd, "rm %s", linxFilePath);
            system(cmd);
            
            linxFile = getFileToWrite(linxFilePath);
            fprintf(linxFile, "@ Q1 %%le %e\n@ Q1RMS %%le %e\n",
                    tunesum[0] / count0[0], sqrt(tune2sum[0] / count0[0] - (tunesum[0] / count0[0]) * (tunesum[0] / count0[0])));
            fclose(linxFile);
            sprintf(cmd, "cat  %s/tx >> %s ; rm %s/tx", workingDirectoryPath, linxFilePath, workingDirectoryPath);
            system(cmd);
            printf("linx %e %d %e\n", tunesum[0], count0[0], tune2sum[0]);
        }
        if (count0[1] > 0) {
            sprintf(cmd, "head -2 %s > %s/ty", linyFilePath, workingDirectoryPath);
            system(cmd);
            sprintf(cmd, "grep \\\" %s | sort -n --key=3 >> %s/ty", linyFilePath, workingDirectoryPath);
            system(cmd);
            sprintf(cmd, "rm %s", linyFilePath);
            system(cmd);

            linyFile = getFileToWrite(linyFilePath);
            fprintf(linyFile, "@ Q2 %%le %e\n@ Q2RMS %%le %e\n",
                    tunesum[1] / count0[1], sqrt(tune2sum[1] / count0[1] - (tunesum[1] / count0[1]) * (tunesum[1] / count0[1])));
            fclose(linyFile);
            sprintf(cmd, "cat  %s/ty >> %s ; rm %s/ty", workingDirectoryPath, linyFilePath, workingDirectoryPath);
            system(cmd);
            printf("liny %e %d %e\n", tunesum[1], count0[1], tune2sum[1]);
        }

/*     1     2     3       4   5      6      7   8       9      10    11   12*/
/*lin:HBPM LABEL HTUNE HPHASE HAMP -20AMP -20PH 01AMP -30AMP -40AMP 20AMP 20PH*/
/*    LABEL VTUNE VAMP VPHASE 10AMPv 0-2AMPv 0-3AMPv  -1-1AMPv   02AMPh POSH  POSV PPH PPV 01PH 10PHv*/
/*     13   14    15    16     17     18      19       20         21    22    23   24  25  26   27*/

        if (labelrun == 1) fclose(noiseFile);
        exit(0);
    }
    if (runBetaBeat == 1) {
        printf("\nBETAbeating process running\n");
        betabeat();
    }
    
    return 0;
}

/* ***************** */
/*    GET_NAME       */
/* ***************** */
void get_name()
{
    char s[1000];
    size_t charCounter = 0;
    FILE *drivingTermsFile = getFileToRead(drivingTermsFilePath);

    /* this block reads the filepath */
    for (charCounter = 0; isspace(dataFilePath[0] = getc(drivingTermsFile)); charCounter++) ;
    for (charCounter = 1; ((charCounter < sizeof(dataFilePath) && (dataFilePath[charCounter] = getc(drivingTermsFile)) != EOF) &&
         (dataFilePath[charCounter] != '\n') && (dataFilePath[charCounter] != '%') &&
         (dataFilePath[charCounter] != ' ')); charCounter++) ;
    if (charCounter >= sizeof(dataFilePath))
    {
        printf("Error: path longer than sizeof(dataFilePath): %u", sizeof(dataFilePath));
        exit(EXIT_FAILURE);
    }
    if (dataFilePath[charCounter] == EOF) {
        scum = 1;
        dataFilePath[charCounter] = '\0';
        fclose(drivingTermsFile);
        return;
    }
    dataFilePath[charCounter] = '\0';
    
    /* this block reads over the first number (tbach) */
    for (charCounter = 0; isspace(s[0] = getc(drivingTermsFile)); charCounter++);
    for (charCounter = 1; ((s[charCounter] = getc(drivingTermsFile)) != EOF) &
         (s[charCounter] != '\n') & (s[charCounter] != ' '); charCounter++) ;
    if (s[charCounter] == EOF)
        scum = 1;
    s[charCounter] = '\0';

    /* this block reads the second number (tbach) */
    for (charCounter = 0; isspace(s[0] = getc(drivingTermsFile)); charCounter++) ;
    for (charCounter = 1; ((charCounter < sizeof(s)) && (s[charCounter] = getc(drivingTermsFile)) != EOF) &&
         (s[charCounter] != '\n') && (s[charCounter] != ' '); charCounter++) ;
    if (charCounter >= sizeof(s))
    {
        printf("Error: input longer than sizeof(s): %u", sizeof(s));
        exit(EXIT_FAILURE);
    }
    if (s[charCounter] == EOF)
        scum = 1;
    s[charCounter] = '\0';
    turns = atoi(s);
    
    fclose(drivingTermsFile);
}

/* ***************** */
/*    sussix_inp     */
/* ***************** */
void sussix_inp(int ir, char *sussixInputFilePath)
{
    FILE *sussixInputFile = getFileToWrite(sussixInputFilePath);
    fprintf(sussixInputFile, "C\nC INPUT FOR SUSSIX_V4 ---17/09/1997---\n");
    fprintf(sussixInputFile, "C DETAILS ARE IN THE MAIN PROGRAM SUSSIX_V4.F\nC\n\n");
    fprintf(sussixInputFile, "ISIX  = 0\nNTOT  = 1\nIANA  = 1\nICONV = 0\n");
    fprintf(sussixInputFile, "TURNS = 1 %d\n", turns);
    fprintf(sussixInputFile, "NARM  = 160\nISTUN = 1 %e %e\n", istun, istun);
    fprintf(sussixInputFile, "TUNES = %e %e .07\n", tunex, tuney);
    fprintf(sussixInputFile, "NSUS  = 0\nIDAM  = %d\n", 2);
    fprintf(sussixInputFile, "NTWIX = 1\nIR    = %d\nIMETH = 2\nNRC   = 4\nEPS   = 2D-3\n", ir); /* EPS is the window in the secondary lines, very imp!!! */
    fprintf(sussixInputFile, "NLINE = 0\nL,M,K = \nIDAMX = 1\nNFIN  = 500\nISME  = 0\n");
    fprintf(sussixInputFile, "IUSME = 200\nINV   = 0\nIINV  = 250\nICF   = 0\nIICF  = 350\n");
    fclose(sussixInputFile);
}

void get_nameb();
char *namefileb = "BPMfiles", datafileb[100];
int iparameterb = 0, scumb = 0;

/* ***************** */
/*    BetaBeat       */
/* ***************** */
void betabeat()
{
    char s[100];
    int i = 0, j, ii = 0, jj = 0, horizontalBpmCounter = 0, verticalBpmCounter = 0,
        flag = 0, contx[MAXPICK], conty[MAXPICK], xpick[MAXPICK], ypick[MAXPICK];
    double ampx[MAXPICK], ampy[MAXPICK], sumampx[MAXPICK],
        sum2ampx[MAXPICK], sumampy[MAXPICK], sum2ampy[MAXPICK], ax, ay;
    FILE *df, *bf;
    for (j = 0; j < MAXPICK; j++) {
        contx[j] = 0;
        conty[j] = 0;
        xpick[j] = 0;
        ypick[j] = 0;
        sumampx[j] = 0;
        sum2ampx[j] = 0;
        sumampy[j] = 0;
        sum2ampy[j] = 0;
        ampx[j] = 0;
        ampy[j] = 0;
    }

    if (!canOpenFile(namefileb))
        printf("no file %s", namefileb);

    while (scumb == 0) {
        get_nameb();
        if (!canOpenFile(datafileb)) {
            printf("No file %s", datafileb);
            break;
        } else
            printf("opening %s\n", datafileb);

        df = getFileToRead(datafileb);
        jj = 0;
        i = 0;
        ii = 0;
        horizontalBpmCounter = 0;
        verticalBpmCounter = 0;
        ax = 0;
        ay = 0;

        while (isspace(s[0] = getc(df)));

        while (s[0] != EOF) {
            if (s[0] == '\n') {
                jj++;
                ii = 0;
            }
            if ((s[0] == '\n' || s[0] == ' ' || s[0] == '\t') && flag == 1)
                flag = 0;
            if ((s[0] != '\n' && s[0] != ' ' && s[0] != '\t') && flag == 0) {
                while (!isspace(s[i]) && (s[i] != EOF)) {
                    i++;
                    if (i >= 100) {
                        printf("Found a value which has >= 100 characters, exit parsing."
                            "This is most probably a malformatted file: %s.", s);
                        exit(EXIT_FAILURE);
                    }
                    s[i] = getc(df);
                }
                s[i + 1] = s[i];
                s[i] = '\0';

                if (ii == 0)
                    xpick[horizontalBpmCounter] = atoi(s);
                if (ii == 2 && i >= 1 && s[i - 1] != 'n') {
                    ampx[xpick[horizontalBpmCounter]] = atof(s);
                    ax = ax + ampx[xpick[horizontalBpmCounter]];
                    printf("ASSIGN %e\t", ampx[horizontalBpmCounter]);
                    horizontalBpmCounter++;
                }
                if (ii == 5)
                    ypick[verticalBpmCounter] = atoi(s);
                if (ii == 7 && i >= 1 && s[i - 1] != 'n') {
                    ampy[ypick[verticalBpmCounter]] = atof(s);
                    ay = ay + ampy[ypick[verticalBpmCounter]];
                    verticalBpmCounter++;
                }

                ii++;
                flag = 1;
                s[0] = s[i + 1];
                i = 0;
            }
            if (flag == 0)
                s[0] = getc(df);
        }
        fclose(df);
        ay = ay / verticalBpmCounter;
        ax = ax / horizontalBpmCounter;

        for (i = 0; i < horizontalBpmCounter; i++) {
            sumampx[xpick[i]] = sumampx[xpick[i]] + ampx[xpick[i]] / ax;
            sum2ampx[xpick[i]] = sum2ampx[xpick[i]] + ampx[xpick[i]] * ampx[xpick[i]] / (ax * ax);
            contx[xpick[i]]++;
        }
        for (i = 0; i < verticalBpmCounter; i++) {
            sumampy[ypick[i]] = sumampy[ypick[i]] + ampy[ypick[i]] / ay;
            sum2ampy[ypick[i]] = sum2ampy[ypick[i]] + ampy[ypick[i]] * ampy[ypick[i]] / (ay * ay);
            conty[ypick[i]]++;
        }
        iparameterb++;
    }
    bf = getFileToWrite("betabeating"); /* file without any path does not make sense? */
    for (i = 0; i <= horizontalBpmCounter; i++) {
        fprintf(bf, "%d  %e  %e  %d  %e  %e\n",
                xpick[i], sumampx[xpick[i]] / contx[xpick[i]],
                sqrt(sum2ampx[xpick[i]] / contx[xpick[i]] - sumampx[xpick[i]] * sumampx[xpick[i]] / (contx[xpick[i]] * contx[xpick[i]])),
                ypick[i], sumampy[ypick[i]] / conty[ypick[i]],
                sqrt(sum2ampy[ypick[i]] / conty[ypick[i]] - sumampy[ypick[i]] * sumampy[ypick[i]] / (conty[ypick[i]] * conty[ypick[i]])));
    }
    fclose(bf);
}

/* ***************** */
/*    GET_NAMEB      */
/* ***************** */
void get_nameb()
{
    int i = 0;
    FILE *file = getFileToRead(namefileb);
    while (i < iparameterb) {
        if ((datafileb[0] = getc(file)) == '\n')
            i++;
    }
    for (i = 0; isspace(datafileb[0] = getc(file)); i++) ;
    for (i = 1; ((datafileb[i] = getc(file)) != EOF) &
         (datafileb[i] != '\n') && (datafileb[i] != '%') &
         (datafileb[i] != ' '); i++) ;
    if (datafileb[i] == EOF) {
        scumb = 1;
        datafileb[i] = '\0';
        return;
    }
    datafileb[i] = '\0';
    fclose(file);
}

/************   GET_WORD *************/
/* get next string in between blanks */
/*************************************/

char *get_word(FILE * ff)
{
    int i = 1;
    while (isspace(ss[0] = getc(ff))) ;
    while (!isspace(ss[i] = getc(ff)))
        ++i;
    ss[i] = '\0';
    return &ss[0];
}

/************   BPMstatus *************/
/* Analyse fort.300 to detect noise   */
/**************************************/
#define MINSIGNAL 0.00001
#define SIGMACUT   1.8
#define MINIMUMNOISE 0.0
#define BADPICKUP  8.0
#define MAXSIGNAL 30000
int BPMstatus(int plane)
{
    double aux, ave = 0, amp,
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
        if (plane == 2) {
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

    if ((maxfreq < windowa2 && maxfreq > windowa1)
        || (maxfreq < windowb2 && maxfreq > windowb1))
        printf("NoiseWindow includes largest lines, amp %e freq %e!!!!\n",
               maxpeak, maxfreq);

    if (maxpeak <= noiseAve + SIGMACUT * noise1)
        return 0;

    if (noise1 > BADPICKUP)
        return 0;

    /* Otherwise pick-up succeded to first cuts */
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
    FILE *file;
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
