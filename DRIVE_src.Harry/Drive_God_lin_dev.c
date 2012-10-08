/* Parallelised version of Drive_God_lin.c of 15/05/2011 using Openmp.
   NB this version requires the DrivingTerms input file to define the
   correct number of turns required - must not be more than in the data
   else bpm values of zero will be used.
   Last Change 29/03/2012 at line 43: remove window variables from the
   threadprivate pragma. These are global constants which were undefined
   other than in the primary thread so noise1, co and co2 were being
   calculated as zero in all secondary threads.
   Change 29/09/2011 at lines 715 and 724 to find lines with any
   bpm name to sort in order by looking for a " rather than a name string. 
   Has matching sussix4drivexxNoO.f      H.Renshall & E.Maclean

   Version <x> 20121001 (tbach):
   - fixed comment line reading (by rtomas)
   - changed variable names to more meaningful and readable names
   - fixed errors from static code analysis
   - removed outcommented code
   - removed unused variables and functions
   - formatted the source
   - changed error messages to more helpful content for the user
   Hint: use ICC (13.0 tested) to compile. GCC has problems with the critical pragma sections.
   */
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <ctype.h>
#define MAXPICK 1100
#define MAXTURNS 4000           /* Critial parameter in sussixfordiveNoO.f ! */
#define MAXTURNS4 16000         /*Always four times  MAXTURNS */
#define MAXRUNS 100
#define M_PI 3.14159265358979310862446895044

#if !defined(LOG_INFO)
#define LOG_INFO 0              /*set to 1 to enable log prints (tbach) */
#endif

char *get_word(FILE *);
void get_name(), sussix_inp(int, char *), betabeat(), savefort300(int);
int BPMstatus(int, int), norm(int, int);
int omp_get_thread_num(void);
int canOpenFile(const char*);
void assertSmaller(const int, const int, const char*);
FILE* getFileToWrite(const char*);
FILE* getFileToRead(const char*);
FILE* __getFileWithMode( const char*, const char*, const char*);

char namefile[2000], datafile[200], labelfile[200], noisefile[200],
    datafileaux[200], inputfile[2000], ss[1000], labelfilepath[200],
    *bpmname[MAXPICK];

double matrix[MAXPICK][MAXTURNS], bpmpos[MAXPICK], tunex, tuney, tune[2],
    amplitude[14], phase[14], parameter[MAXRUNS], istun, noise1, noiseAve,
    maxpeak, maxfreq, maxmin, co, co2, windowa1, windowa2, windowb1,
    windowb2, allfreqsx[300], allampsx[300], allfreqsy[300], allampsy[300],
    doubleToSend[MAXTURNS4 + 4], allbpmamp[MAXPICK], allbpmphase[MAXPICK];

int turns, normalisation, iparameter = 0, scum = 0, in = 0, iir = 0,
    labelrun = 0, label[MAXPICK], NoVerticalData = 0, hv[MAXPICK],
    hvt[MAXPICK], nslines;

#pragma omp threadprivate(doubleToSend,tune,tunex,tuney,amplitude,phase,\
        parameter, istun, noise1, noiseAve,\
        maxpeak, maxfreq, maxmin, co, co2,\
        allfreqsx,allampsx,allfreqsy,allampsy)

main(int argc, char **argv)
{
    double amax, amp, backtune[2], backphase[3],
        backamp[14], offrms, tunesum[2], tune2sum[2],
        zsum, ysum, zysum, z2sum, amp10bpm[MAXPICK], amp30bpm[MAXPICK],
        kper, amp01bpm[MAXPICK], amp30ave[MAXRUNS], amp01ave[MAXRUNS],
        amp00ave[MAXRUNS], amp1sum[2], amp1sum2[2], amp2sum[2],
        amp2sum2[2], phase01ave[MAXRUNS], phase02sum[2], amp00bpm[MAXPICK],
        amp00sum[2], amp00sum2[2], raizbeta[262], amp4sum[2], amp4sum2[2],
        amp5sum[2], amp6sum[2], amp40bpm[MAXPICK], amp50bpm[MAXPICK],
        amp60bpm[MAXPICK], amp40ave[MAXRUNS], amp50ave[MAXRUNS], amp5sum2,
        amp60ave[MAXRUNS], amp20bpm[MAXPICK], amp20sum[2], amp20sum2[2],
        amp20ave[MAXRUNS], meantune, fact, upperpeak, lowerpeak,
        peaktopeak[MAXPICK], averamp[2], averamp2[2];

    int i = 0, kick = 0, ij, finii, finij, mini, counth = 0, countv = 0, maxcounthv,
        pickupaxis1[MAXPICK][2], pickupaxis2[MAXPICK][2], start, flag = 0,
        columnCounter = 0, j, horizontalBpmCounter = 0, verticalBpmCounter = 0,
        count0[2], count1[2], count2[2], count8[2], count9[2], pickstart = 0,
        pickend = 0, kcase = 0, ibetabeat = 0, opticlabel = 1, bpmCounter = 0,
        lastSlashIndex, count4[2], count5[2], count6[2], iim, iformat = 0,
        ipick, ipickm, backlai, backlaij, kk, Nbpms, Nturns = 0, avercount[2],
        counthl = 0, countvl = 0;

    char s[300], bpmfile[300], linfile[300], linfilex[2000],
        linfiley[2000], *opticfile = "monitors-hpi40.txt", *pointer,
        *ps, cmd[6000], path[4000], bsfile[400];

    FILE *df, *di, *fb, *optic, *lf, *nf, *kf, *linfx, *linfy, *fBS;
    omp_set_dynamic(0);
    /* Memory allocation */
    for (i = 0; i < MAXPICK; i++)
        bpmname[i] = (char *) calloc(50, sizeof(char));

/*  Path to DrivingTerms and Drive.inp */
    printf("\n Path to Drive input: %s\n", argv[1]);

    /*TODO create size dynamically? (tbach)*/
    assertSmaller(strlen(argv[1]) + strlen("/sussix_v4.inp"), sizeof(path), "modify path");
    strncpy(path, argv[1], sizeof(path)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
    strncat(path, "/sussix_v4.inp", sizeof(path) - strlen(path) - 1);

    assertSmaller(strlen(argv[1]) + strlen("/DrivingTerms"), sizeof(namefile), "modify namefile");
    strncpy(namefile, argv[1], sizeof(namefile)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
    strncat(namefile, "/DrivingTerms", sizeof(namefile) - strlen("/DrivingTerms") - 1);

    assertSmaller(strlen(argv[1]) + strlen("/Drive.inp"), sizeof(inputfile), "modify inputfile");
    strncpy(inputfile, argv[1], sizeof(inputfile)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
    strncat(inputfile, "/Drive.inp", sizeof(inputfile) - strlen(inputfile) - 1);

    printf("\n Drive.inp: %s\n", inputfile);



/* check the file namefile */
    if (!canOpenFile(namefile)) {
        printf("\nNo file %s for reading the name of the Data file\n", namefile);
        exit(EXIT_FAILURE);
    }

/* check the input file Drive.inp */
    if (!canOpenFile(inputfile)) {
        printf("\nNo input file %s\n", inputfile);
        exit(EXIT_FAILURE);
    }


    count0[0] = 0;
    count0[1] = 0;
    tunesum[0] = 0;
    tune2sum[0] = 0;
    tunesum[1] = 0;
    tune2sum[1] = 0;
    /* check the OPTIC  file and take Beta(s) in the monitors if it exists */
    if (!canOpenFile(opticfile)) {
        /* doesn't exist --> forget */
        printf("\nNO OPTIC FILE: %s\n", opticfile);
        opticlabel = 0;

    } else {
        optic = getFileToRead(opticfile);
        for (i = 0; i < 262; i++) {
            raizbeta[i] = sqrt(atof(get_word(optic)));
            get_word(optic);
            get_word(optic);
            get_word(optic);

        }
        printf("\nOPTIC FILE: %s USED\n", opticfile);
        fclose(optic);
    }


    i = 0;
    di = getFileToRead(inputfile);
    while (i <= 17) {
        *s = '\0';
        if ((s[0] = getc(di)) == '=') {
            i++;
            for (bpmCounter = 0; ((s[bpmCounter] = getc(di)) != EOF) & (s[bpmCounter] != '\n') & (s[bpmCounter] != ' ');
                 ++bpmCounter);
            s[bpmCounter] = '\0';
            switch (i) {
            case 1:
                kick = atof(s) - 1;     /* C arrays start at 0 */
                break;
            case 2:
                kcase = atof(s);
                break;
            case 3:
                kper = atof(s);
                break;
            case 4:
                tunex = atof(s);
                break;
            case 5:
                tuney = atof(s);
                break;
            case 6:
                pickstart = atoi(s);
                break;
            case 7:
                pickend = atoi(s);
                break;
            case 8:
                normalisation = atoi(s);
                break;
            case 9:
                istun = atof(s);
                break;
            case 10:
                ibetabeat = atoi(s);
                break;
            case 11:
                iir = atoi(s);
                break;
            case 12:
                labelrun = atoi(s);
                break;
            case 13:
                iformat = atoi(s);
                break;
            case 14:
                strncpy(labelfilepath, s, sizeof(labelfile));
                labelfilepath[sizeof(labelfile) - 1] = '\0';
                break;
            case 15:
                windowa1 = atof(s);
                break;
            case 16:
                windowa2 = atof(s);
                break;
            case 17:
                windowb1 = atof(s);
                break;
            case 18:
                windowb2 = atof(s);
                break;
            }
        }
    }
    fclose(di);
    if (kick >= 0)
        printf("Known kick in %d turn\n", kick + 1);
    if (kcase == 1)
        printf("Horizontal case\n");
    else if (kcase == 0)
        printf("Vertical case\n");
    else {
        printf("No proper case in Drive.inp\n");
        exit(EXIT_FAILURE);
    }

    if (labelrun == 1)
        printf("\n LABELRUN: NOISE FILES WILL BE WRITTEN TO NOISEPATH\n");
    printf("Normalisation: %d\n", normalisation);
    printf("pickstart: %d, pickend: %d\n", pickstart, pickend);
    if (pickstart < 0 || pickstart > pickend || pickstart > MAXPICK)
    {
        printf("Bad value for pickstart. Must be >= 0, < pickend, <= %d", MAXPICK);
        exit(EXIT_FAILURE);
    }


    while (scum == 0) {
        /* From namefile assign datafile, increase iparameter,
           asign parameter[iparameter] and turns.               */
        get_name();

        /* Check the file datafile */
        if (!canOpenFile(datafile)) {
            /* doesn't exist --> stop */
            printf("\nNo Data file %s\n", datafile);
            break;
        }
        printf("Data File: %s\n", datafile);
        /*constructing name of BPM files and labelfile */
        lastSlashIndex = 0;
        if (strrchr(datafile, '/') != NULL)
            lastSlashIndex = strrchr(datafile, '/') - datafile; /* search last occurence, substract pointer. we search 2 times here, who cares (tbach) */

        /* copy everything from behind the last slash until the end to bpmfile (tbach) */
        assertSmaller(strlen(datafile + lastSlashIndex), sizeof(bpmfile), "modify bpmfile");
        strncpy(bpmfile, datafile + lastSlashIndex, sizeof(bpmfile)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */

        /* If labelfilepath does not end with '/' put it */
        if (strrchr(labelfilepath, '/') == NULL || /* if none found, try to add a slash (tbach) */
            strrchr(labelfilepath, '/') - labelfilepath != strlen(labelfilepath) - 1) /* search last occurence, substract pointer, have index (tbach) */
        {
            assertSmaller(strlen(labelfilepath) + 1, sizeof(labelfilepath), "modify labelfilepath");
            strncat(labelfilepath, "/", sizeof(labelfilepath) - strlen(labelfilepath) - 1);
        }

        assertSmaller(strlen(labelfilepath) + strlen(bpmfile) + strlen("_label"), sizeof(labelfile), "modify labelfile");
        strncpy(labelfile, labelfilepath, sizeof(labelfile)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
        strncat(labelfile, bpmfile, sizeof(labelfile) - strlen(labelfile) - 1);
        strncat(labelfile, "_label", sizeof(labelfile) - strlen(labelfile) - 1);

        assertSmaller(strlen(labelfilepath) + strlen(bpmfile) + strlen("_noise"), sizeof(noisefile), "modify noisefile");
        strncpy(noisefile, labelfilepath, sizeof(noisefile)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
        strncat(noisefile, bpmfile, sizeof(noisefile) - strlen(noisefile) - 1);
        strncat(noisefile, "_noise", sizeof(noisefile) - strlen(noisefile) - 1);

        assertSmaller(strlen(argv[1]) + 1 + strlen(bpmfile) + strlen("_linx"), sizeof(linfilex), "modify linfilex");
        strncpy(linfilex, argv[1], sizeof(linfilex)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
        strncat(linfilex, "/", sizeof(linfilex) - strlen(linfilex) - 1);
        strncat(linfilex, bpmfile, sizeof(linfilex) - strlen(linfilex) - 1);
        strncat(linfilex, "_linx", sizeof(linfilex) - strlen(linfilex) - 1);

        assertSmaller(strlen(argv[1]) + 1 + strlen(bpmfile) + strlen("_liny"), sizeof(linfilex), "modify linfiley");
        strncpy(linfiley, argv[1], sizeof(linfiley)); /* This could produce a not null terminated result, which is prevented by the assert (tbach) */
        strncat(linfiley, "/", sizeof(linfiley) - strlen(linfiley) - 1);
        strncat(linfiley, bpmfile, sizeof(linfiley) - strlen(linfiley) - 1);
        strncat(linfiley, "_liny", sizeof(linfiley) - strlen(linfiley) - 1);

        assertSmaller(strlen(bpmfile) + strlen("_bpm"), sizeof(bpmfile), "modify bpmfile");
        strncat(bpmfile, "_bpm", sizeof(bpmfile) - strlen(bpmfile) - 1);
        
        
        linfx = getFileToWrite(linfilex);
        linfy = getFileToWrite(linfiley);
        fprintf(linfx,
                "* NAME  S   BINDEX SLABEL  TUNEX   MUX  AMPX   NOISE   PK2PK  AMP01 PHASE01 CO CORMS AMP_20  PHASE_20  AMP02 PHASE02 AMP_30  PHASE_30\n");
        fprintf(linfx,
                "$ %%s  %%le   %%le    %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le    %%le     %%le   %%le     %%le   %%le     %%le ");
        fprintf(linfy,
                "* NAME  S   BINDEX SLABEL  TUNEY  MUY  AMPY   NOISE   PK2PK AMP10 PHASE10 CO CORMS AMP_1_1 PHASE_1_1\n");
        fprintf(linfy,
                "$ %%s  %%le  %%le    %%le  %%le  %%le  %%le  %%le  %%le %%le  %%le %%le  %%le %%le %%le ");
        if (labelrun == 1) {
            lf = getFileToWrite(labelfile);
            nf = getFileToWrite(noisefile);
        }

        /* Constructing a matrix with all the data from the pick-ups for iformat i=0,1 */
        if (iformat < 2) {
            bpmCounter = 0;
            i = 0;
            df = getFileToRead(datafile);
            while ((s[0] = getc(df)) != '\n') ;
            s[0] = getc(df);
            while (s[0] != EOF) {
                if (s[0] == '\n') {
                    ++bpmCounter;
                    columnCounter = 0;
                }
                if ((s[0] == '\n' || s[0] == ' ' || s[0] == '\t') && flag == 1)
                    flag = 0;
                if ((s[0] != '\n' && s[0] != ' ' && s[0] != '\t') && flag == 0) {
                    while (!isspace(s[i]) && s[i] != EOF) {
                        ++i;
                        s[i] = getc(df);
                    }
                    s[i + 1] = s[i];
                    s[i] = '\0';
                    iim = columnCounter / 2;

                    if (opticlabel == 0)
                        fact = 1.0;
                    else
                        fact = raizbeta[columnCounter];

                    if (iformat == 0) {
                        if (iim == (columnCounter / 2.0))
                            matrix[iim + MAXPICK / 2][bpmCounter] = atof(s) / fact;
                        else
                            matrix[iim][bpmCounter] = atof(s) / fact;
                    } else
                        matrix[columnCounter][bpmCounter] = atof(s) / fact;
                    ++columnCounter;
                    flag = 1;
                    s[0] = s[i + 1];
                    i = 0;
                }
                if (flag == 0)
                    s[0] = getc(df);
            }

            fclose(df);

            /* Label spare pick-ups from datafile & labelfile */
            for (i = 0; i < MAXPICK; i++)
                label[i] = 0;
            df = getFileToRead(datafile);
            for (i = 0; i < MAXPICK; i++) {
                while (isspace(s[0] = getc(df))) ;
                if (s[0] == 'B') {
                    iim = i / 2;
                    if (iformat == 0)
                        if (iim == (i / 2.0))
                            label[iim + MAXPICK / 2] = 1;
                        else
                            label[iim] = 1;
                    if (iformat == 1)
                        label[i] = 1;
                }
                while (!isspace(s[0] = getc(df))) ;
            }
            fclose(df);

        }

        /* Closes if iformat<2 */
        /* Constructing a matrix with all the data from the pick-ups for iformat i=2, RHIC */
        flag = 0;
        for (i = 0; i < MAXPICK; i++)
            label[i] = 0;
        if (iformat == 2) {
            bpmCounter = 0;
            columnCounter = 0;
            horizontalBpmCounter = -1;
            verticalBpmCounter = MAXPICK / 2 - 1;
            i = 0;
            df = getFileToRead(datafile);
            s[0] = getc(df);
            while (s[0] == '#') {       /* then it is a comment line (tbach) */
                while (getc(df) != '\n');       /* read until the end of the line (tbach) */
                s[0] = getc(df);        /* read the first char of the new line (tbach) */
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
                        s[i] = getc(df);
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
                        if (opticlabel == 0)
                            fact = 1.0;
                        else
                            fact = raizbeta[columnCounter];
                        if (hv[bpmCounter] == 0)
                            matrix[horizontalBpmCounter][columnCounter - 3] = atof(s) / fact;
                        else
                            matrix[verticalBpmCounter][columnCounter - 3] = atof(s) / fact;
                        Nturns = columnCounter - 3;     /* If the last line is an empty line, then we can get the number of turns only from here. First 3 are plane, name and location. (tbach) */
                    }
                    ++columnCounter;
                    flag = 1;
                    s[0] = s[i + 1];
                    i = 0;
                }
                if (flag == 0)
                    s[0] = getc(df);
            }
            fclose(df);
            Nbpms = bpmCounter;
            counth = horizontalBpmCounter + 1;
            countv = verticalBpmCounter + 1;

            /* Some statistics and checks */
            printf("Total number of pick-ups: %d Last turn number: %d\n", Nbpms, Nturns);
            printf("Horizontal pick-ups: %d   Vertical pick-ups: %d\n ", counth, -MAXPICK / 2 + countv);
            printf("name of pick-up 0: %s pos: %f first turn:%f, second turn: %f\n",
                 bpmname[0], bpmpos[0], matrix[0][0], matrix[0][1]);
        }/* CLoses if iformat==2 */

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
                if (abs(matrix[start][bpmCounter] - matrix[start][bpmCounter - 1]) > kper) {
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
        printf("Turns to be procesed: %d %f\n", turns, matrix[0][0]);


        /* Measuring the peak-to-peak amplitude in the first 10 turns after kick */
        averamp[0] = 0.0;
        averamp[1] = 0.0;
        averamp2[0] = 0.0;
        averamp2[1] = 0.0;
        avercount[0] = 0;
        avercount[1] = 0;
        for (i = 0; i < MAXPICK; i++) {
            upperpeak = -10000;
            lowerpeak = 10000;
            if (label[i] == 1) {        /* This is not really needed but ok... */
                for (j = 0; j < 10; j++) {
                    if (matrix[i][j] == 99)
                        printf("99 found in %d %d %s\n", i, j, bpmname[i]);
                    if (matrix[i][j] > upperpeak)
                        upperpeak = matrix[i][j];
                    if (matrix[i][j] < lowerpeak)
                        lowerpeak = matrix[i][j];
                }
                peaktopeak[i] = upperpeak - lowerpeak;
            }
        }

        /* First part of the analysis: Determine  phase of all pick-ups and noise */
        sussix_inp(1, path);

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

            sussix4drivenoise_(&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0], path);
            /* This calls the external fortan code (tbach) */

            #pragma omp critical
            {
                allbpmphase[columnCounter] = phase[0];
                allbpmphase[ij] = phase[3];
                allbpmamp[columnCounter] = amplitude[0];
                allbpmamp[ij] = amplitude[3];
                label[columnCounter] = BPMstatus(1, columnCounter);
                if (labelrun == 1)
                    fprintf(nf, "1 %d  %e %e %e %e %e %d %d %f\n",
                            columnCounter, noise1, noiseAve, maxpeak,
                            maxfreq, maxmin, nslines, label[i],
                            phase[0] / 360.);

                /* PRINT LINEAR FILE */
                if (amplitude[0] > 0 && label[i] == 1 && columnCounter == i) {
                    fprintf(linfx, "\n\"%s\" %e %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ",
                            bpmname[columnCounter], bpmpos[columnCounter],
                            columnCounter, label[columnCounter], tune[0],
                            phase[0] / 360., amplitude[0], noise1, maxmin,
                            amplitude[2] / amplitude[0], phase[2] / 360.,
                            co, co2, amplitude[1] / amplitude[0],
                            phase[1] / 360., amplitude[12] / amplitude[0],
                            phase[12] / 360., amplitude[6] / amplitude[0],
                            phase[6] / 360.);
                    ++count0[0];
                    tunesum[0] += tune[0];
                    tune2sum[0] += tune[0] * tune[0];
                    /* Horizontal Spectrum output */
                    if (i < 10) {
                        sprintf(bsfile, "%s/%s.x", argv[1], bpmname[i]);
                        fBS = getFileToWrite(bsfile);
                        fprintf(fBS, "%s %s %s\n", "*", "FREQ", "AMP");
                        fprintf(fBS, "%s %s %s\n", "$", "%le", "%le");
                        for (kk = 0; kk < 300; ++kk)
                            fprintf(fBS, "%e %e\n", allfreqsx[kk], allampsx[kk]);
                        fclose(fBS);
                    }
                }
                label[ij] = BPMstatus(2, ij);
                if (labelrun == 1)
                    fprintf(nf, "2 %d  %e %e %e %e %e %d %d %f\n", ij,
                            noise1, noiseAve, maxpeak, maxfreq, maxmin,
                            nslines, label[ij], phase[3] / 360.);
                if (amplitude[3] > 0 && label[ij] == 1 && ij == i + MAXPICK / 2) {
                    fprintf(linfy, "\n\"%s\" %e %d %d %e %e %e %e %e %e %e %e %e %e %e",
                            bpmname[ij], bpmpos[ij], ij, label[ij],
                            tune[1], phase[3] / 360., amplitude[3], noise1,
                            maxmin, amplitude[5] / amplitude[3],
                            phase[5] / 360., co, co2,
                            amplitude[13] / amplitude[3],
                            phase[13] / 360.);
                    ++count0[1];
                    tunesum[1] += tune[1];
                    tune2sum[1] += tune[1] * tune[1];
                    if (ij < MAXPICK / 2 + 10) {
                        sprintf(bsfile, "%s/%s.y", argv[1], bpmname[ij]);
                        fBS = getFileToWrite(bsfile);
                        fprintf(fBS, "%s %s %s\n", "*", "FREQ", "AMP");
                        fprintf(fBS, "%s %s %s\n", "$", "%le", "%le");
                        for (kk = 0; kk < 300; ++kk)
                            fprintf(fBS, "%e %e \n", allfreqsy[kk], allampsy[kk]);
                        fclose(fBS);
                    }
                }
            }                   /* end of omp critical section */
        }
        fprintf(linfx, "\n");
        fprintf(linfy, "\n");
        fclose(linfx);
        fclose(linfy);
        /* What follows is some dirty operations to move the @ Q1 ... to the top of _linx _liny */
        if (count0[0] > 0) {
            /* HRR June 2011 parallel version change from mv to sort on bpm number:  sprintf(cmd, "mv %s %s/tx", linfilex,argv[1]); */
            sprintf(cmd, "head -2 %s > %s/tx", linfilex, argv[1]);
            system(cmd);
            /* HRR Sep 2011 change grep -i bpm to grep \" as only bpm lines contain a " */
            sprintf(cmd, "grep \\\" %s | sort -n --key=3 >> %s/tx",
                    linfilex, argv[1]);
            system(cmd);
            sprintf(cmd, "rm %s", linfilex);
            system(cmd);
        }
        if (count0[1] > 0) {
            /* HRR June 2011 parallel version change from mv to sort on bpm number:  sprintf(cmd, "mv %s %s/ty", linfiley,argv[1]); */
            sprintf(cmd, "head -2 %s > %s/ty", linfiley, argv[1]);
            system(cmd);
            sprintf(cmd, "grep \\\" %s | sort -n --key=3 >> %s/ty",
                    linfiley, argv[1]);
            system(cmd);
            sprintf(cmd, "rm %s", linfiley);
            system(cmd);
        }
        if (count0[0] > 0) {
            linfx = getFileToWrite(linfilex);
            fprintf(linfx, "@ Q1 %%le %e\n@ Q1RMS %%le %e\n",
                    tunesum[0] / count0[0],
                    sqrt(tune2sum[0] / count0[0] - (tunesum[0] / count0[0]) * (tunesum[0] / count0[0])));
            fclose(linfx);
            sprintf(cmd, "cat  %s/tx >> %s ; rm %s/tx", argv[1], linfilex, argv[1]);
            system(cmd);
            printf("linx %e %d %e\n", tunesum[0], count0[0], tune2sum[0]);
        }
        if (count0[1] > 0) {
            linfy = getFileToWrite(linfiley);
            fprintf(linfy, "@ Q2 %%le %e\n@ Q2RMS %%le %e\n",
                    tunesum[1] / count0[1],
                    sqrt(tune2sum[1] / count0[1] - (tunesum[1] / count0[1]) * (tunesum[1] / count0[1])));
            fclose(linfy);
            sprintf(cmd, "cat  %s/ty >> %s ; rm %s/ty", argv[1], linfiley, argv[1]);
            system(cmd);
            printf("liny %e %d %e\n", tunesum[1], count0[1], tune2sum[1]);
        }

/*     1     2     3       4   5      6      7   8       9      10    11   12*/
/*lin:HBPM LABEL HTUNE HPHASE HAMP -20AMP -20PH 01AMP -30AMP -40AMP 20AMP 20PH*/
/*    LABEL VTUNE VAMP VPHASE 10AMPv 0-2AMPv 0-3AMPv  -1-1AMPv   02AMPh POSH  POSV PPH PPV 01PH 10PHv*/
/*     13   14    15    16     17     18      19       20         21    22    23   24  25  26   27*/

        if (labelrun == 1) {
            fclose(nf);
            fclose(lf);
        }
        exit(0);                /*FIXME this does not make any sense? (that we have _dead_ code after this) (tbach) */
    }
    if (ibetabeat == 1) {
        printf("\nBETAbeating process running\n");
        betabeat();
    }
}

/* ***************** */
/*    GET_NAME       */
/* ***************** */
void get_name()
{
    char s[1000];
    int i = 0;
    FILE *file = getFileToRead(namefile);
    while (i < iparameter) { /* FIXME This does not make any sense, iparameter is always 0 at this point (tbach) */
        if ((datafile[0] = getc(file)) == '\n')
            i++;
    }
    for (i = 0; isspace(datafile[0] = getc(file)); i++) ;
    for (i = 1; ((datafile[i] = getc(file)) != EOF) &&
         (datafile[i] != '\n') && (datafile[i] != '%') &&
         (datafile[i] != ' '); i++) ;

    if (datafile[i] == EOF) {
        scum = 1;
        datafile[i] = '\0';
        iparameter++;
        return;
    }
    datafile[i] = '\0';
    for (i = 0; isspace(s[0] = getc(file)); i++);
    for (i = 1; ((s[i] = getc(file)) != EOF) &
         (s[i] != '\n') & (s[i] != ' '); i++) ;
    if (s[i] == EOF)
        scum = 1;
    s[i] = '\0';
    parameter[iparameter] = atof(s);

    for (i = 0; isspace(s[0] = getc(file)); i++) ;
    for (i = 1; ((s[i] = getc(file)) != EOF) &
         (s[i] != '\n') & (s[i] != ' '); i++) ;
    if (s[i] == EOF)
        scum = 1;
    s[i] = '\0';
    turns = atof(s);
    fclose(file);
}

/* ***************** */
/*    sussix_inp     */
/* ***************** */
void sussix_inp(int ir, char *path)
{
    FILE *file = getFileToWrite(path);
    fprintf(file, "C\nC INPUT FOR SUSSIX_V4 ---17/09/1997---\n");
    fprintf(file, "C DETAILS ARE IN THE MAIN PROGRAM SUSSIX_V4.F\nC\n\n");
    fprintf(file, "ISIX  = 0\nNTOT  = 1\nIANA  = 1\nICONV = 0\n");
    fprintf(file, "TURNS = 1 %d\n", turns);
    fprintf(file, "NARM  = 160\nISTUN = 1 %e %e\n", istun, istun);
    fprintf(file, "TUNES = %e %e .07\n", tunex, tuney);
    fprintf(file, "NSUS  = 0\nIDAM  = %d\n", 2);
    fprintf(file, "NTWIX = 1\nIR    = %d\nIMETH = 2\nNRC   = 4\nEPS   = 2D-3\n", ir);     /* EPS is the window in the secondary lines, very imp!!! */
    fprintf(file, "NLINE = 0\nL,M,K = \nIDAMX = 1\nNFIN  = 500\nISME  = 0\n");
    fprintf(file, "IUSME = 200\nINV   = 0\nIINV  = 250\nICF   = 0\nIICF  = 350\n");
    fclose(file);
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
                    s[i] = getc(df);
                }
                s[i + 1] = s[i];
                s[i] = '\0';

                if (ii == 0)
                    xpick[horizontalBpmCounter] = atof(s);
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
    bf = getFileToWrite("betabeating");
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
    char s[80];
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
int BPMstatus(int plane, int pickupN)
{
    double norma = 0.0, weight, aux, freq = 0.0, freq2 = 0.0, ave = 0, amp,
        maxe = -500000.0, mine = 500000.0;
    int counter, il, counter2 = 0, counter3 = 0;
    FILE *fort, *f90;
    char line[500], ssw[80], str[80];

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
    return __getFileWithMode(filename, "w", "Cannot open to write: %s");
}
FILE* getFileToRead(const char* const filename)
{
    return __getFileWithMode(filename, "r", "Cannot open to read: %s");
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
    printf("Value1: %d is not < than Value2: %d. Message: %s", a, b, message);
    exit(EXIT_FAILURE);
}
