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



//======================================================================================================================
// Declarations and definitions
//======================================================================================================================
class InputData{
	void check_inp_data();
	void check_tune(double& tune, char plane);
	void readInputDataFromFileDriveInp(std::string driveInputFilePath);

public:
    void initAndReadInputValues(std::string& driveInpFilePath);

    double  istun, kper, tunex, tuney, windowa1, windowa2, windowb1, windowb2, nattunex, nattuney;
    int turns, pickStart, pickEnd, kcase, kick, labelrun;
}inpData;


class IoHelper{
	static bool inputFileCheck(std::string filePath);
	static bool outputFileCheck(std::string filePath);
public:
	static void getAndCheckFilePaths(std::string *drivingTermsFilePath, std::string *driveInputFilePath, const char * cmdinput);
	static bool cannotOpenFile(std::string filePath,char type);
	static void writeSussixInputFile(const int, const double, const double, const double);
	static void formatLinFile(std::string linFilePath, const int tunecount, const double tunesum, const double tune2sum,
			const int nattunecount, const double nattunesum, const double nattune2sum, int plane_index);

	static std::string workingDirectoryPath;
	static std::string sussixInputFilePath;
};
std::string IoHelper::workingDirectoryPath;
std::string IoHelper::sussixInputFilePath;


class OutFilesHandler{
public:
	OutFilesHandler(std::string &dataFilePath);
	~OutFilesHandler();

	void writeLine();

	std::ofstream linxFile, linyFile, noiseFile, rejectedBpmsFileX, rejectedBpmsFileY;
	std::string linxFilePath, linyFilePath;
};


class TuneCalcData{
public:
    double tuneSumX, tuneSumY, tune2sumX, tune2sumY, nattuneSumX, nattune2sumX, nattuneSumY, nattune2sumY;
    int    tuneCountX, tuneCountY, nattuneCountX, nattuneCountY;
    void init_calc_values(){ //initializes all values to be calculated to zero
        tuneSumX = tuneSumY = tune2sumX = tune2sumY = nattuneSumX = nattune2sumX = nattuneSumY = nattune2sumY = 0.0;
        tuneCountX = tuneCountY = nattuneCountX = nattuneCountY = 0;
    }
    void addTuneX(double tune, double natTune) {
    	addTune(tuneCountX, tuneSumX, tune2sumX, nattuneCountX, nattuneSumX, nattune2sumX, tune, natTune);
    }
    void addTuneY(double tune, double natTune) {
    	addTune(tuneCountY, tuneSumY, tune2sumY, nattuneCountY, nattuneSumY, nattune2sumY, tune, natTune);
    }
private:
    void addTune(int& tuneCount, double& tuneSum, double& tune2Sum, int& natTuneCount, double& natTuneSum, double& natTune2Sum,
    		double& tune, double& natTune) {
    	++tuneCount;
    	tuneSum += tune;
    	tune2Sum += tune * tune;
		if (natTune > NATTUNE_DEFAULT) { /*  Initialized to -100. Condition true if nat tune found */
			++natTuneCount;
			natTuneSum += natTune;
			natTune2Sum += natTune * natTune;
		}
    }
}tuneCalcData;


void harmonicAnalysisForSingleTbtDataFile(std::string&);
void readTbtDataFile(std::string&, int&, int&, int&);
void findKick();
bool readLineInDrivingTerms(std::istream&, int*, std::string *);
bool BPMstatus(const int, const int);

// Functions inside of the omp parallel for loop
inline void callExternalFortranFunctionForHarmonicAnalysis();
inline void calculateNaturalTune();
inline void writeLineIntoLinXFile(std::ofstream& linXFile, const int& horizontalBpmIndex);
inline void writeLineIntoLinYFile(std::ofstream& linYFile, const int& verticalBpmIndex);
inline void createSpectrumFileForCurrentHorizontalBpm(const int& horizontalBpmIndex);
inline void createSpectrumFileForCurrentVerticalBpm(const int& verticalBpmIndex);
inline void writeLineIntoRejecetedBpmsFileX(std::ofstream& rejectedBpmsFileX, int& horizontalBpmIndex);
inline void writeLineIntoRejecetedBpmsFileY(std::ofstream& rejectedBpmsFileY, int& verticalBpmIndex);

double calculatednattuney, calculatednattunex, calculatednatampy, calculatednatampx, co, co2,
	   noise1, maxfreq, maxmin, maxpeak, noiseAve;

double allampsx[300], allampsy[300], allfreqsx[300], allfreqsy[300], amplitude[19],
		   doubleToSend[MAXTURNS4 + 4], phase[19], tune[2];

int nslines=0;

struct BPM{ /*Structure for each BPM- has name, position, plane, if it's pickedup, and it's tbtdata*/
	std::string bpmName;
	double bpmPos;
	int plane; /* 0 is horizontal, 1 is vertical*/
	bool pickedUp;
	double tbtData[MAXTURNS];
}BPMs[MAXPICK];

#pragma omp threadprivate(amplitude, doubleToSend, tune, phase,\
        noise1, noiseAve, maxpeak, maxfreq, maxmin, co, co2,\
        allfreqsx, allampsx, allfreqsy, allampsy)


//======================================================================================================================
// main function
//======================================================================================================================
int main(int argc, char **argv)
{

    #ifdef _WIN32 /*Changes minor formatting difference in windows regarding the output of a number in scientific notation.*/
        _set_output_format(_TWO_DIGIT_EXPONENT);
    #endif
    //To output scientific notation
    std::cout << std::setiosflags (std::ios::scientific);

    omp_set_dynamic(0);


    std::string drivingTermsFilePath, driveInpFilePath;

    IoHelper::getAndCheckFilePaths(&drivingTermsFilePath, &driveInpFilePath, argv[1]);

    inpData.initAndReadInputValues(driveInpFilePath);


    std::ifstream drivingTermsFile;
    drivingTermsFile.open(drivingTermsFilePath.c_str());
    std::string dataFilePath;
	int loopcounter = 0;
    while (readLineInDrivingTerms(drivingTermsFile, &(inpData.turns), &dataFilePath)) {
        /* Check the file dataFilePath */
        if(IoHelper::cannotOpenFile(dataFilePath,'i'))
            continue;
        loopcounter++;
        harmonicAnalysisForSingleTbtDataFile(dataFilePath);
    } /* end of while loop over all files to analyse */
    drivingTermsFile.close();

    if (loopcounter == 0)
        std::cerr << "Drivingterms file has bad input, no data ever read\n";

    return EXIT_SUCCESS;
}


void harmonicAnalysisForSingleTbtDataFile(std::string &dataFilePath){
    /* set all values to be calculated to default values */
    tuneCalcData.init_calc_values();
	OutFilesHandler filesHandler(dataFilePath);
	int maxOfHorBpmsAndVerBpms = 0; /** max(numOfHorBpms, numOfVerBpms) */
	int lastBpmIndexHor = 0;
	int lastBpmIndexVer = 0;
	readTbtDataFile(dataFilePath, maxOfHorBpmsAndVerBpms, lastBpmIndexHor, lastBpmIndexVer);

	findKick();

	IoHelper::writeSussixInputFile(inpData.turns, inpData.istun, inpData.tunex, inpData.tuney);

#pragma omp parallel for private(calculatednattunex, calculatednattuney, calculatednatampx, calculatednatampy)
	for (int i = inpData.pickStart; i < maxOfHorBpmsAndVerBpms; ++i) {
		int horizontalBpmIndex = i;
		int verticalBpmIndex = i + MAXPICK/2;

		if (verticalBpmIndex > lastBpmIndexVer)
			verticalBpmIndex = lastBpmIndexVer;
		if (horizontalBpmIndex > lastBpmIndexHor)
			horizontalBpmIndex = lastBpmIndexHor;
		if (horizontalBpmIndex < 0 || verticalBpmIndex < MAXPICK/2)
		{
			fprintf(stderr, "horizontal < 0 or vertical BpmIndex < MAXPICK/2. Should not happen.\n");
			exit(EXIT_FAILURE);
		}

		for (int j = 0; j < MAXTURNS; ++j) {
			doubleToSend[j] = BPMs[horizontalBpmIndex].tbtData[j];
			doubleToSend[j + MAXTURNS] = BPMs[verticalBpmIndex].tbtData[j];
			doubleToSend[j + 2 * MAXTURNS] = 0.0;
			doubleToSend[j + 3 * MAXTURNS] = 0.0;
		}

		callExternalFortranFunctionForHarmonicAnalysis();

		calculateNaturalTune();

		#pragma omp critical
		{
			BPMs[horizontalBpmIndex].pickedUp = BPMstatus(1, inpData.turns); /*Always returns true*/
			if (inpData.labelrun == 1)
				filesHandler.noiseFile << std::scientific << "1 " << horizontalBpmIndex << "  " <<  noise1 << ' ' <<  noiseAve << ' ' << maxpeak << ' ' << maxfreq << ' ' << maxmin << ' ' << nslines << ' ' << BPMs[i].pickedUp << ' ' << phase[0] / 360. << std::endl;

			/* PRINT LINEAR FILE */
			if (amplitude[0] > 0 && BPMs[horizontalBpmIndex].pickedUp == true && horizontalBpmIndex == i) {
				writeLineIntoLinXFile(filesHandler.linxFile, horizontalBpmIndex);
				tuneCalcData.addTuneX(tune[0], calculatednattunex);
				createSpectrumFileForCurrentHorizontalBpm(horizontalBpmIndex);
			}else if(true == (BPMs[horizontalBpmIndex].pickedUp == true && horizontalBpmIndex == i)){
				writeLineIntoRejecetedBpmsFileX(filesHandler.rejectedBpmsFileX, horizontalBpmIndex);

			}

			BPMs[verticalBpmIndex].pickedUp = BPMstatus(2, inpData.turns); /*Always returns true*/
			if (inpData.labelrun == 1)
				filesHandler.noiseFile << std::scientific << "2 " << verticalBpmIndex << "  " <<  noise1 << ' ' <<  noiseAve << ' ' << maxpeak << ' ' << maxfreq << ' ' << maxmin << ' ' << nslines << ' ' << BPMs[verticalBpmIndex].pickedUp << ' ' << phase[3] / 360. << std::endl;

			if (amplitude[3] > 0 && BPMs[verticalBpmIndex].pickedUp == true && verticalBpmIndex == i + MAXPICK/2) {
				writeLineIntoLinYFile(filesHandler.linyFile, verticalBpmIndex);
				tuneCalcData.addTuneY(tune[1], calculatednattuney);
				createSpectrumFileForCurrentVerticalBpm(verticalBpmIndex);
			}else if(true == (BPMs[verticalBpmIndex].pickedUp == true && verticalBpmIndex == i + MAXPICK/2)){
				writeLineIntoRejecetedBpmsFileY(filesHandler.rejectedBpmsFileY, verticalBpmIndex);
			}
		} /* end of omp critical section */
	} /* end of parallel for */

	/* Sort and move the "@..." lines to the top of the _linx/y files */
	IoHelper::formatLinFile(filesHandler.linxFilePath,
			tuneCalcData.tuneCountX, tuneCalcData.tuneSumX, tuneCalcData.tune2sumX, tuneCalcData.nattuneCountX, tuneCalcData.nattuneSumX, tuneCalcData.nattune2sumX, 1);
	IoHelper::formatLinFile(filesHandler.linyFilePath,
			tuneCalcData.tuneCountY, tuneCalcData.tuneSumY, tuneCalcData.tune2sumY, tuneCalcData.nattuneCountY, tuneCalcData.nattuneSumY, tuneCalcData.nattune2sumY, 2);
}


/**
 * Reads the turn-by-turn file and stores each BPM into struct array BPMs.
 * @param maxOfHorBpmsAndVerBpms: 1 + max(numOfReadHorBpms, numOfReadVerBpms)
 * @param lastUsedIndexHor: The last used index of horizontal BPMs in BPMs struct array will be stored here
 * @param lastUsedIndexVer: The last used index of vertical BPMs in BPMs struct array will be stored here
 */
void readTbtDataFile(std::string& dataFilePath, int& maxOfHorBpmsAndVerBpms, int& lastUsedIndexHor, int& lastUsedIndexVer) {
	for (int i = 0; i < MAXPICK; i++)
			BPMs[i].pickedUp = false;

	char tempStr[1000];
	std::ifstream dataFile;
	dataFile.open(dataFilePath.c_str());
	tempStr[0] = (char) (dataFile.get());
	while (tempStr[0] == '#') {
		/* then it is a comment line (tbach) */
		while (dataFile.get() != '\n')
			; /* read until the end of the line (tbach) */

		tempStr[0] = (char) (dataFile.get()); /* read the first char of the new line (tbach) */
	}
	/* after this, we have skipped all the comment lines, and s[0] is the first character of a new line which is not a "#" (tbach) */
	if (LOG_INFO)
		printf("BPM file content:\n");

	int i = 0;
	int flag = 0;
	lastUsedIndexHor = -1;
	lastUsedIndexVer = MAXPICK / 2 - 1;
	int columnCounter = 0;
	int numTurns = 0;
	int minNumTurns = -1;
	int bpmCounter = 0;
	int currentPlane = -1;
	while (tempStr[0] != EOF) {
		if (tempStr[0] == '\n') {
			++bpmCounter;
			if (LOG_INFO)
				printf("\n");
			columnCounter = 0;
			if(numTurns != 0)
				if(-1 == minNumTurns)
					minNumTurns = numTurns;
				if(minNumTurns > numTurns){
					// We have a BPM with less turns than previous BPMs (vimaier)
					fprintf(stderr, "##Probably unequal lengths of turns detected.\n"
							"##  lastNumOfTurns: %d\n##  newNumOfTurns: %d\n##  bpmCounter: %d\n",
							minNumTurns, numTurns, bpmCounter-1);
					std::string bpmName;
					if (currentPlane == 0) /* 0 is horizontal, 1 is vertical (tbach) */
						bpmName = BPMs[lastUsedIndexHor].bpmName;
					else
						bpmName = BPMs[lastUsedIndexVer].bpmName;
					fprintf(stderr, "##  Check BPM %s\n", bpmName.c_str());
					minNumTurns = numTurns;
				}
		}
		if (isspace((int) tempStr[0]) && flag == 1)
			flag = 0;
		if (!isspace((int) tempStr[0]) && flag == 0) {
			while (!isspace((int) tempStr[i]) && tempStr[i] != EOF) {
				++i;
				tempStr[i] = (char) dataFile.get();
				if (i > 100) {
					tempStr[i + 1] = '\0';
					fprintf(stderr,"Found a value which has more than 100 characters, exit parsing.\n"
						"This is most probably a malformatted file. bpmCounter=%d columnCounter=%d string1=%s\n", bpmCounter, columnCounter, tempStr);
					exit(EXIT_FAILURE);
				}
			}
			tempStr[i + 1] = tempStr[i];
			tempStr[i] = '\0';
			if (LOG_INFO)
				printf("%s ", tempStr);
			if (columnCounter >= MAXTURNS) {
				fprintf(stderr, "Found >= %d Turns, this turn size is not supported. Reduce amount of turns. bpmCounter:%d\n", MAXTURNS - 3, bpmCounter); /* 0,1,2 is plane, name and location (tbach) */
				exit(EXIT_FAILURE);
			}
			if (bpmCounter >= MAXPICK) {
				fprintf(stderr, "Found >= %d BPMs, this size is not supported. Reduce amount of BPMs. columnCounter:%d\n", MAXPICK, columnCounter);
				exit(EXIT_FAILURE);
			}
			if (columnCounter == 0) { /*plane (tbach) */
				currentPlane = atoi(tempStr);
				if (currentPlane == 0) /* 0 is horizontal, 1 is vertical (tbach) */
					BPMs[++lastUsedIndexHor].plane = currentPlane;
				else
					BPMs[++lastUsedIndexVer].plane = currentPlane;
			}

			else if (columnCounter == 1) { /*bpm name (tbach) */
				if (currentPlane == 0) {
					if (lastUsedIndexHor < 0) /* Branch prediction will cry, but well lets have security (tbach) */
					{
						fprintf(stderr, "horizontalBpmCounter < 0. Should not happen. Probably malformatted input file?\n");
						exit(EXIT_FAILURE);
					}
					BPMs[lastUsedIndexHor].bpmName = tempStr;
					BPMs[lastUsedIndexHor].pickedUp = true;
				} else {
					BPMs[lastUsedIndexVer].bpmName = tempStr;
					BPMs[lastUsedIndexVer].pickedUp = true;
				}
			}

			else if (columnCounter == 2) { /*bpm location (tbach) */
				if (currentPlane == 0)
				{
					if (lastUsedIndexHor < 0) /* Branch prediction will cry, but well lets have security (tbach) */
					{
						fprintf(stderr, "horizontalBpmCounter < 0. Should not happen. Probably malformatted input file?\n");
						exit(EXIT_FAILURE);
					}
					BPMs[lastUsedIndexHor].bpmPos = atof(tempStr);
				}
				else
				BPMs[lastUsedIndexVer].bpmPos = atof(tempStr);
			}

			else { /*bpm data (tbach) */
				numTurns = columnCounter - 3 + 1;
				/* If the last line is an empty line, then we can get the number of turns only from here.
				 First 3 are plane, name and location.
				 Plus 1 for index start at 0
				 (tbach) */
				if (currentPlane == 0)
					BPMs[lastUsedIndexHor].tbtData[numTurns-1] = atof(tempStr);
				else
					BPMs[lastUsedIndexVer].tbtData[numTurns-1] = atof(tempStr);
			}
			++columnCounter;
			flag = 1;
			tempStr[0] = tempStr[i + 1];
			i = 0;
		}
		if (flag == 0)
			tempStr[0] = (char) dataFile.get();
	}
	dataFile.close();

	int counth = lastUsedIndexHor + 1;
	int countv = lastUsedIndexVer + 1;

	/* now redefine turns as the minimum of the Nturns read and the DrivingTerms data card */
	/* NB assumes all BPMs have the same number of turns as the last one read is used */
	numTurns = minNumTurns;
	if (inpData.turns > numTurns) inpData.turns = numTurns;

	/* Some statistics and checks */
	printf("Total number of pick-ups: %d Last turn number: %d, turns to run: %d\n", bpmCounter, numTurns, inpData.turns);
	printf("Horizontal pick-ups: %d   Vertical pick-ups: %d\n", counth, -MAXPICK / 2 + countv);
	printf("name of BPM[0]: %s, pos: %f, first turn: %f, second turn: %f, last turn: %f, last turn to run: %f \n",
		 BPMs[0].bpmName.c_str(), BPMs[0].bpmPos, BPMs[0].tbtData[0], BPMs[0].tbtData[1], BPMs[0].tbtData[numTurns - 1], BPMs[0].tbtData[inpData.turns - 1]);

	if (counth >= (countv - MAXPICK / 2))
		maxOfHorBpmsAndVerBpms = counth;
	else
		maxOfHorBpmsAndVerBpms = -MAXPICK / 2 + countv;

	if (maxOfHorBpmsAndVerBpms > inpData.pickEnd)
		maxOfHorBpmsAndVerBpms = inpData.pickEnd;

	if (maxOfHorBpmsAndVerBpms >= MAXPICK) {
		fprintf(stderr, "\nNot enough Pick-up memory\n");
		exit(EXIT_FAILURE);
	}
	printf("BPMs in loop: %d, pickstart: %d, resulting loop length: %d\n",
		 maxOfHorBpmsAndVerBpms, inpData.pickStart, maxOfHorBpmsAndVerBpms - inpData.pickStart);
}

void findKick() {
	printf("kick: %d \n", inpData.kick);
	/* searching for two working adjacent pick-ups */
	/* after the Q-kickers for finding the kick */
	if (inpData.kick < 0) {
		int start = -(inpData.kcase - 1) * MAXPICK / 2 + 2;
		while (BPMs[start].pickedUp == false
				|| BPMs[start + 2].pickedUp == false) {
			start = start + 2;
		}

		printf("looking for kick in pick-up:%d\n", start + 1);
		/* Find kick here and get kick */
		for (int columnCounter = 1;
				(inpData.kick < 0) && (columnCounter < inpData.turns);
				++columnCounter) {
			if (fabs(
					BPMs[start].tbtData[columnCounter]
							- BPMs[start].tbtData[columnCounter - 1])
					> inpData.kper) {
				inpData.kick = columnCounter;
			}
		}

		if (inpData.kick < 0) {
			fprintf(stderr, "NO KICK FOUND\n");
			exit(EXIT_FAILURE);
		} else
		printf("Found kick in turn:%d\n", inpData.kick + 1); /*Natural count */
	}

	if (inpData.kick > 0) {
		// Shift array by inpData.kick
		for (int i = 0; i < MAXPICK; i++) {
			if (BPMs[i].pickedUp == true) {
				for (int j = inpData.kick; j < inpData.turns; j++)
					BPMs[i].tbtData[j - inpData.kick] = BPMs[i].tbtData[j];
			}
		}
		inpData.turns -= inpData.kick;
	}
	printf(
			"Turns to be processed after kick offset: %d BPMs[0].tbtdata[0]: %f \n",
			inpData.turns, BPMs[0].tbtData[0]);
}


/* *****************  */
/*    readLineInDrivingTerms*/
/* *****************  */
bool readLineInDrivingTerms(std::istream& drivingTermsFile, int *turns, std::string *path)
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
        return false;

    if(OS == "Windows32") //Path may not exist if drivingtermsfile was created on linux, so this deals with that. (Assumes /afs/cern.ch is 'H' drive)
        if((pos = (*path).find("/afs/cern.ch")) != std::string::npos)
            *path = "H:"+(*path).substr(pos+strlen("/afs/cern.ch"));

    if(OS == "linux"){ //Same story as above, replaces H: with /afs/cern.ch, also changes '\' to '/'
        while((pos = (*path).find("\\")) != std::string::npos)
            (*path)[pos] = '/';
        if((pos = (*path).find("H:")) != std::string::npos)
            *path = "/afs/cern.ch"+(*path).substr(pos+strlen("H:"));
    }

	return true;
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

        if ((aux < inpData.windowa2 && aux > inpData.windowa1)
            || (aux < inpData.windowb2 && aux > inpData.windowb1)) {
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

    if ((inpData.windowa1 < maxfreq && maxfreq < inpData.windowa2)
     || (inpData.windowb1 < maxfreq && maxfreq < inpData.windowb2))
        printf("NoiseWindow includes largest lines, amp %e freq %e!!!!\n",
               maxpeak, maxfreq);

    return true;
}

//======================================================================================================================
// IoHelper definitions
//======================================================================================================================
void IoHelper::getAndCheckFilePaths(std::string *drivingTermsFilePath, std::string *driveInputFilePath, const char * cmdinput){

    workingDirectoryPath = cmdinput;

    std::cout << "Working directory path: " << workingDirectoryPath << std::endl;

    if(OS == "linux" && cannotOpenFile(workingDirectoryPath,'i')){ //Always fails to open in windows
        exit(EXIT_FAILURE);
    }

    *drivingTermsFilePath = workingDirectoryPath+"/DrivingTerms";
    *driveInputFilePath = workingDirectoryPath+"/Drive.inp";
    sussixInputFilePath = workingDirectoryPath+"/sussix_v4.inp";

    //check the input files drivingTermsFilePath and Drive.inp
    if(cannotOpenFile(*drivingTermsFilePath,'i')
    || cannotOpenFile(*driveInputFilePath,'i') || cannotOpenFile(sussixInputFilePath,'o')){
        exit(EXIT_FAILURE);
    }
    std::cout << "DrivingTerms: " <<  *drivingTermsFilePath << std::endl;
    std::cout << "Drive.inp: " <<  *driveInputFilePath << std::endl;
    std::cout << "sussix_v4.inp: " <<  sussixInputFilePath << std::endl;
}


bool IoHelper::cannotOpenFile(std::string filePath,char type){
    bool failure;
    if(type == 'o')
        failure = outputFileCheck(filePath);
    else
        failure = inputFileCheck(filePath);
    if(failure){
        std::cerr << "Failed to open file " << filePath << std::endl;
        return failure;
        }
    else
        return failure;
}

bool IoHelper::inputFileCheck(std::string filePath){
    std::ifstream file(filePath.c_str());
    bool failure = file.fail();
    file.close();
    return failure;
}

bool IoHelper::outputFileCheck(std::string filePath){
    std::ofstream file(filePath.c_str());
    bool failure = file.fail();
    file.close();
    return failure;
}

void IoHelper::writeSussixInputFile(const int turns, const double istun, const double tunex, const double tuney)
{
    std::ofstream sussixInputFile(sussixInputFilePath.c_str());
    if(IoHelper::cannotOpenFile(sussixInputFilePath,'o')){
        exit(EXIT_FAILURE);
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

/**
 * Restructures x_lin* file. Brings the descriptors(lines starting with '@') to the top of the file.
 */
void IoHelper::formatLinFile(std::string linFilePath,
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
            std::cerr << "Leaving drive due to error" << std::endl;
            exit(EXIT_FAILURE);
        }
        tempFile << std::scientific << "@ Q" << plane_index << " %le " << tunesum / tunecount << "\n@ Q" << plane_index << "RMS %le " << sqrt(tune2sum / tunecount - (tunesum / tunecount) * (tunesum / tunecount)) << std::endl;
        if (nattunecount > 0)
            tempFile << std::scientific << "@ NATQ" << plane_index << " %le " << nattunesum / nattunecount << "\n@ NATQ" << plane_index << "RMS %le " << sqrt(nattune2sum / nattunecount - (nattunesum / nattunecount) * (nattunesum / nattunecount)) << std::endl;

        /*Gets linFile to read*/
        linFile.open(linFilePath.c_str());
        if(cannotOpenFile(linFilePath,'i')){
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


//======================================================================================================================
// Inputdata definitions
//======================================================================================================================
void InputData::initAndReadInputValues(std::string& driveInpFilePath){ //initializes all input data to 0 except for natural tunes which have default value defined
        istun = kper = tunex = tuney = windowa1 = windowa2 = windowb1 = windowb2 = 0.0;
        turns = pickStart = pickEnd = kcase = kick = labelrun = 0;
        nattunex = nattuney = NATTUNE_DEFAULT;

        //Reads input data from Drive.inp into global InpData class
        this->readInputDataFromFileDriveInp(driveInpFilePath);

        //Checks tune_x/y, kcase, kick, labelrun, pickend, and pickend
        this->check_inp_data();
    }


void InputData::check_inp_data(){

   this->check_tune(this->tunex,'x');
   this->check_tune(this->tuney,'y');

   if (this->kick >= 0)
	   printf("Known kick in turn %d\n", this->kick + 1);
   if (this->kcase == 1)
	   printf("Horizontal case\n");
   else if (this->kcase == 0)
	   printf("Vertical case\n");
   else {
	   fprintf(stderr, "No proper kcase in Drive.inp\n");
	   exit(EXIT_FAILURE);
   }

   if (this->labelrun == 1)
	   printf("\n LABELRUN: NOISE FILES WILL BE WRITTEN TO NOISEPATH\n");

   printf("pickStart: %d, pickEnd: %d\n", this->pickStart, this->pickEnd);
   if (this->pickStart < 0 || this->pickStart > this->pickEnd || this->pickStart > MAXPICK) {
	   fprintf(stderr, "Bad value for pickStart. Must be >= 0 and < pickEnd and <= MAXPICK(=%d)\n", MAXPICK);
	   exit(EXIT_FAILURE);
   }
}

void InputData::check_tune(double& tuneValue, char plane){ //checks that magnitude of tune is less than 0.5, changes it by integer value until it is
	int change;
	if(tuneValue > 0.5 || tuneValue < -0.5){
		change=0;
		while(tuneValue > 0.5){
			tuneValue -= 1;
			change--;
		}
		while(tuneValue < -0.5){
			tuneValue += 1;
			change++;
		}
		std::cout << "tune_" << plane <<" input increased by " << change << ". Should be less than or equal to 0.5 in magnitude. New value is: " << tuneValue << std::endl;
	}
}

void InputData::readInputDataFromFileDriveInp(std::string driveInputFilePath){
	std::string temp_str;
	std::ifstream driveInputFile;
	unsigned int pos;
	driveInputFile.open(driveInputFilePath.c_str());
	while(!driveInputFile.rdstate()){
		std::getline (driveInputFile, temp_str);
		if((pos = temp_str.find("KICK=")) != std::string::npos)
			this->kick = atoi(temp_str.substr(pos+strlen("KICK=")).c_str())-1;
		if((pos = temp_str.find("CASE(1[H], 0[V])=")) != std::string::npos)
			this->kcase = atoi(temp_str.substr(pos+strlen("CASE(1[H], 0[V])=")).c_str());
		if((pos = temp_str.find("KPER(KICK PERCE.)=")) != std::string::npos)
			this->kper = atof(temp_str.substr(pos+strlen("KPER(KICK PERCE.)=")).c_str());
		if((pos = temp_str.find("TUNE X=")) != std::string::npos)
			this->tunex = atof(temp_str.substr(pos+strlen("TUNE X=")).c_str());
		if((pos = temp_str.find("TUNE Y=")) != std::string::npos)
			this->tuney = atof(temp_str.substr(pos+strlen("TUNE Y=")).c_str());
		if((pos = temp_str.find("PICKUP START=")) != std::string::npos)
			this->pickStart = atoi(temp_str.substr(pos+strlen("PICKUP START=")).c_str());
		if((pos = temp_str.find("PICKUP END=")) != std::string::npos)
			this->pickEnd = atoi(temp_str.substr(pos+strlen("PICKUP END=")).c_str());
		if((pos = temp_str.find("ISTUN=")) != std::string::npos)
			this->istun = atof(temp_str.substr(pos+strlen("ISTUN=")).c_str());
		if((pos = temp_str.find("LABEL RUN (1[yes])=")) != std::string::npos)
			this->labelrun = atoi(temp_str.substr(pos+strlen("LABEL RUN (1[yes])=")).c_str());
		if((pos = temp_str.find("WINDOWa1=")) != std::string::npos)
			this->windowa1 = atof(temp_str.substr(pos+strlen("WINDOWa1=")).c_str());
		if((pos = temp_str.find("WINDOWa2=")) != std::string::npos)
			this->windowa2 = atof(temp_str.substr(pos+strlen("WINDOWa2=")).c_str());
		if((pos = temp_str.find("WINDOWb1=")) != std::string::npos)
			this->windowb1 = atof(temp_str.substr(pos+strlen("WINDOWb1=")).c_str());
		if((pos = temp_str.find("WINDOWb2=")) != std::string::npos)
			this->windowb2 = atof(temp_str.substr(pos+strlen("WINDOWb2=")).c_str());
		if((pos = temp_str.find("NATURAL X=")) != std::string::npos)
			this->nattunex = atof(temp_str.substr(pos+strlen("NATURAL X=")).c_str());
		if((pos = temp_str.find("NATURAL Y=")) != std::string::npos)
			this->nattuney = atof(temp_str.substr(pos+strlen("NATURAL Y=")).c_str());
	}
	driveInputFile.close();
}

//======================================================================================================================
// OutFilesHandler definitions
//======================================================================================================================
OutFilesHandler::OutFilesHandler(std::string &dataFilePath){
	std::cout << "Data file: " << dataFilePath << std::endl;

	std::string bpmFileName = dataFilePath.substr(dataFilePath.rfind('/',dataFilePath.size()-2)+1);

	std::cout << "bpmFileName: " << bpmFileName << std::endl;

	std::string noiseFilePath = IoHelper::workingDirectoryPath+'/'+bpmFileName+"_noise";
	linxFilePath = IoHelper::workingDirectoryPath+'/'+bpmFileName+"_linx";
	linyFilePath = IoHelper::workingDirectoryPath+'/'+bpmFileName+"_liny";

	rejectedBpmsFileX.open((IoHelper::workingDirectoryPath+'/'+bpmFileName+"_rejectedBpms_x").c_str());
	rejectedBpmsFileX << "* NAME S    BINDEX SLABEL TUNEX MUX  AMPX NOISE PK2PK AMP01 PHASE01 CO   CORMS AMP_20 PHASE_20 AMP02 PHASE02 AMP_30 PHASE_30 AMP_1_1 PHASE_1_1 AMP2_2 PHASE2_2 AMP0_2 PHASE0_2 NATTUNEX NATAMPX\n";
	rejectedBpmsFileX << std::scientific << "$ %s  %le %le   %le   %le  %le %le %le  %le  %le  %le    %le %le  %le   %le     %le  %le    %le   %le     %le    %le      %le   %le     %le   %le     %le     %le\n";
	rejectedBpmsFileY.open((IoHelper::workingDirectoryPath+'/'+bpmFileName+"_rejectedBpms_y").c_str());
	rejectedBpmsFileY << "* NAME S    BINDEX SLABEL TUNEY MUY  AMPY NOISE PK2PK AMP10 PHASE10 CO   CORMS AMP_1_1 PHASE_1_1 AMP_20 PHASE_20 AMP1_1 PHASE1_1 AMP0_2 PHASE0_2 AMP0_3 PHASE0_3 NATTUNEY NATAMPY\n";
	rejectedBpmsFileY << std::scientific << "$ %s  %le %le   %le   %le  %le %le %le  %le  %le  %le    %le %le  %le    %le      %le   %le     %le   %le     %le   %le     %le   %le     %le       %le\n";

	if (inpData.labelrun == 1) {
		noiseFilePath = IoHelper::workingDirectoryPath+'/'+bpmFileName+"_bpm_noise";
		noiseFile.open(noiseFilePath.c_str());
		if(IoHelper::cannotOpenFile(noiseFilePath,'o')){
			exit(EXIT_FAILURE);
		}
	}
	if(IoHelper::cannotOpenFile(linxFilePath,'o') || IoHelper::cannotOpenFile(linyFilePath,'o')){
		exit(EXIT_FAILURE);
	}
	linxFile.open(linxFilePath.c_str());
	linxFile << "* NAME S    BINDEX SLABEL TUNEX MUX  AMPX NOISE PK2PK AMP01 PHASE01 CO   CORMS AMP_20 PHASE_20 AMP02 PHASE02 AMP_30 PHASE_30 AMP_1_1 PHASE_1_1 AMP2_2 PHASE2_2 AMP0_2 PHASE0_2 NATTUNEX NATAMPX\n";
	linxFile << std::scientific << "$ %s  %le %le   %le   %le  %le %le %le  %le  %le  %le    %le %le  %le   %le     %le  %le    %le   %le     %le    %le      %le   %le     %le   %le     %le     %le\n";
	linyFile.open(linyFilePath.c_str());
	linyFile << "* NAME S    BINDEX SLABEL TUNEY MUY  AMPY NOISE PK2PK AMP10 PHASE10 CO   CORMS AMP_1_1 PHASE_1_1 AMP_20 PHASE_20 AMP1_1 PHASE1_1 AMP0_2 PHASE0_2 AMP0_3 PHASE0_3 NATTUNEY NATAMPY\n";
	linyFile << std::scientific << "$ %s  %le %le   %le   %le  %le %le %le  %le  %le  %le    %le %le  %le    %le      %le   %le     %le   %le     %le   %le     %le   %le     %le       %le\n";
}

OutFilesHandler::~OutFilesHandler() {
	rejectedBpmsFileX.close();
	rejectedBpmsFileY.close();
	linxFile.close();
	linyFile.close();
	if (inpData.labelrun == 1)
		noiseFile.close();
}


//======================================================================================================================
// Helper functions in omp parallel for loop
//======================================================================================================================
inline void callExternalFortranFunctionForHarmonicAnalysis() {
	/* This calls the external Fortran code (tbach)-Different name depending on OS (asherman)*/
	#ifdef _WIN32
		SUSSIX4DRIVENOISE (&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0], (char*)IoHelper::sussixInputFilePath.c_str());
	#else
		sussix4drivenoise_(&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0], (char*)IoHelper::sussixInputFilePath.c_str());
	#endif

}

inline void calculateNaturalTune() {
	/* Let's look for natural tunes in the istun range if natural tunes input is given*/
	double maxamp = 0.0;
	calculatednattunex = NATTUNE_DEFAULT;
	if (inpData.nattunex > NATTUNE_DEFAULT) {
		for (int j = 0; j < 300; ++j) {
			if ((inpData.nattunex - inpData.istun < allfreqsx[j] && allfreqsx[j] < inpData.nattunex + inpData.istun) && (maxamp < allampsx[j])) {
				maxamp = allampsx[j];
				calculatednattunex = allfreqsx[j];
				calculatednatampx = maxamp;
			}
		}
	}
	maxamp = 0;
	calculatednattuney = NATTUNE_DEFAULT;
	if (inpData.nattuney > NATTUNE_DEFAULT) {
		for (int j = 0; j < 300; ++j) {
			if ((inpData.nattuney - inpData.istun < allfreqsy[j] && allfreqsy[j] < inpData.nattuney + inpData.istun) && (maxamp < allampsy[j])) {
				maxamp = allampsy[j];
				calculatednattuney = allfreqsy[j];
				calculatednatampy = maxamp;
			}
		}
	}
}


inline void writeLineIntoLinXFile(std::ofstream& linXFile, const int& horizontalBpmIndex){
	linXFile <<  std::scientific << '"' << BPMs[horizontalBpmIndex].bpmName << "\" " << BPMs[horizontalBpmIndex].bpmPos << ' ' << horizontalBpmIndex << ' ' << BPMs[horizontalBpmIndex].pickedUp << ' ' << tune[0] << ' ' <<
			phase[0] / 360. << ' ' << amplitude[0] << ' ' << noise1 << ' ' << maxmin << ' ' << amplitude[2] / amplitude[0] << ' ' << phase[2] / 360. << ' ' << co << ' ' << co2 << ' ' << amplitude[1] / amplitude[0] << ' ' <<
			phase[1] / 360. << ' ' << amplitude[12] / amplitude[0] << ' ' << phase[12] / 360. << ' ' << amplitude[6] / amplitude[0] << ' ' <<
			phase[6] / 360. << ' ' << amplitude[14] / amplitude[0]  << ' ' << phase[14] / 360. << ' ' << amplitude[16] / amplitude[0] << ' ' <<
			phase[16] / 360. << ' ' << amplitude[18] / amplitude[0] << ' ' << phase[18] / 360. << ' ' << calculatednattunex << ' ' << calculatednatampx << std::endl;

}
inline void writeLineIntoLinYFile(std::ofstream& linYFile, const int& verticalBpmIndex){
	linYFile <<  std::scientific << '"' << BPMs[verticalBpmIndex].bpmName << "\" " << BPMs[verticalBpmIndex].bpmPos << ' ' << verticalBpmIndex << ' ' << BPMs[verticalBpmIndex].pickedUp << ' ' << tune[1] << ' ' <<
							phase[3] / 360. << ' ' << amplitude[3] << ' ' << noise1 << ' ' << maxmin << ' ' << amplitude[5] / amplitude[3] << ' ' << phase[5] / 360. << ' ' << co << ' ' << co2 << ' ' <<
							amplitude[13] / amplitude[3] << ' ' << phase[13] / 360. << ' ' << amplitude[15] / amplitude[3] << ' ' << phase[15] / 360. << ' ' <<
							amplitude[17] / amplitude[3] << ' ' << phase[17] / 360. << ' ' << amplitude[4] / amplitude[3] << ' ' << phase[4] / 360. << ' ' <<
							amplitude[11] / amplitude[3] << ' ' << phase[11] / 360. << ' ' << calculatednattuney << ' ' << calculatednatampy << std::endl;

}

inline void writeSpectrumToFile(std::string& spectrumFilePath, const double* allFrequencies, const double* allAmplitudes) {
	if(IoHelper::cannotOpenFile(spectrumFilePath,'o')){
		exit(EXIT_FAILURE);
	}
	std::ofstream spectrumFile(spectrumFilePath.c_str());
	spectrumFile << "* FREQ AMP\n$ %le %le\n";
	for (int j = 0; j < 300; ++j)
		spectrumFile << std::scientific << allFrequencies[j] << ' ' << allAmplitudes[j] << std::endl;
	spectrumFile.close();
}
inline void createSpectrumFileForCurrentHorizontalBpm(const int& horizontalBpmIndex) {
	/* Horizontal Spectrum output */
	if (horizontalBpmIndex >= 10) {
		return;
	}
	std::string spectrumFilePath = IoHelper::workingDirectoryPath+'/'+BPMs[horizontalBpmIndex].bpmName+".x";
	writeSpectrumToFile(spectrumFilePath, allfreqsx, allampsx);
}
inline void createSpectrumFileForCurrentVerticalBpm(const int& verticalBpmIndex) {
	/* Horizontal Spectrum output */
	if (verticalBpmIndex >= MAXPICK/2 + 10) {
		return;
	}
	std::string spectrumFilePath = IoHelper::workingDirectoryPath+'/'+BPMs[verticalBpmIndex].bpmName+".y";
	writeSpectrumToFile(spectrumFilePath, allfreqsy, allampsy);
}


inline void writeLineIntoRejecetedBpmsFileX(std::ofstream& rejectedBpmsFileX, int& horizontalBpmIndex){
	printf("Hor. BPM %s not in lin file because following condition failed: ", BPMs[horizontalBpmIndex].bpmName.c_str());
	printf("false == amplitude[0] > 0\n");
	printf("false == %12f > 0\n", amplitude[0]);
	rejectedBpmsFileX <<  std::scientific << '"' << BPMs[horizontalBpmIndex].bpmName << "\" " << BPMs[horizontalBpmIndex].bpmPos << ' ' << horizontalBpmIndex << ' ' << BPMs[horizontalBpmIndex].pickedUp << ' ' << tune[0] << ' ' <<
								phase[0] / 360. << ' ' << amplitude[0] << ' ' << noise1 << ' ' << maxmin << ' ' << amplitude[2] / amplitude[0] << ' ' << phase[2] / 360. << ' ' << co << ' ' << co2 << ' ' << amplitude[1] / amplitude[0] << ' ' <<
								phase[1] / 360. << ' ' << amplitude[12] / amplitude[0] << ' ' << phase[12] / 360. << ' ' << amplitude[6] / amplitude[0] << ' ' <<
								phase[6] / 360. << ' ' << amplitude[14] / amplitude[0]  << ' ' << phase[14] / 360. << ' ' << amplitude[16] / amplitude[0] << ' ' <<
								phase[16] / 360. << ' ' << amplitude[18] / amplitude[0] << ' ' << phase[18] / 360. << ' ' << calculatednattunex << ' ' << calculatednatampx << std::endl;
}


inline void writeLineIntoRejecetedBpmsFileY(std::ofstream& rejectedBpmsFileY, int& verticalBpmIndex){
	printf("Ver. BPM %s not in lin file because following condition failed: ", BPMs[verticalBpmIndex].bpmName.c_str());
	printf("false == amplitude[3] > 0\n");
	printf("false == %12f > 0\n", amplitude[3]);
	rejectedBpmsFileY <<  std::scientific << '"' << BPMs[verticalBpmIndex].bpmName << "\" " << BPMs[verticalBpmIndex].bpmPos << ' ' << verticalBpmIndex << ' ' << BPMs[verticalBpmIndex].pickedUp << ' ' << tune[1] << ' ' <<
								phase[3] / 360. << ' ' << amplitude[3] << ' ' << noise1 << ' ' << maxmin << ' ' << amplitude[5] / amplitude[3] << ' ' << phase[5] / 360. << ' ' << co << ' ' << co2 << ' ' <<
								amplitude[13] / amplitude[3] << ' ' << phase[13] / 360. << ' ' << amplitude[15] / amplitude[3] << ' ' << phase[15] / 360. << ' ' <<
								amplitude[17] / amplitude[3] << ' ' << phase[17] / 360. << ' ' << amplitude[4] / amplitude[3] << ' ' << phase[4] / 360. << ' ' <<
								amplitude[11] / amplitude[3] << ' ' << phase[11] / 360. << ' ' << calculatednattuney << ' ' << calculatednatampy << std::endl;

}


