#---------------------------------
# Reads SDDS ASCII files
# Removes bad BPM (flat signal, SVD mode over threshold)
# Also removes the noise floor of the data to enhance physical signal
# Ram Calaga, Dec 30, 2003
# @ Glenn => made it possible to filter on space
#
# version 2 20120905 (tbach):
# - major refactoring (introduced real class objects, removed all globals, removed a lot of not necessary switches and conversions, and some things more)
# - order of bad BPMs can be different to previous versions
# - formatting of bad BPMs is different to previous versions (The column order is the same as a normal SDDS ASCII file now, added a reason why the BPM is bad)
# - output formating is different to previous versions (only for whitespaces to make it human readable)
# - Improved speed a lot (biggest change was from O(n*n*n) to O(n*n)). Previous runs had around 20 seconds, now 3 seconds)
# - Data output is in all test cases the same, but can be different if different bad BPM contribute to the same mode in the SVD
#
# version 3 20120913 (tbach):
# - writing for bad BPM changed
# - flat BPM detection made before matrix creation
#
# version 4 20120921 (tbach):
# - only one svd call for every plane needed now (reduces previous 4 calls to now 2 calls)
# - the rearrangeent of the SVD calculation changes some part of the data.
# Previous: Matrix M, do SVD from M, find bad BPMs from SVD, remove this bad BPMs from M, redo SVD from M, cut singular values, create M from cutted SVD
# Now:      Matrix M, do SVD from M, find bad BPMs from SVD, remove this bad BPMs from SVD results,        cut singular values, create M from cutted SVD
# From tests, the data was different, but only with a small derivation which looked ok
#
# version 5 20121015 (tbach):
# - fixed the startturn and maxturn argument handling, did not work in previous versions
# - beautified the print output
# - introduced newfile argument. This can be used to specify a new filename if overriding is not wanted
#
# usage: python svd_clean -f filename
# -> outputs filename (overrides)

#---------------------------------

#------------------------No Changes below this----------------------------------------------
import sys, os, time

from numpy import dot as matrixmultiply

from optparse import OptionParser
import numpy

#-------------Base Class for parsing sdds file and make nturnsxnbpms matrix with tbt--------



# option parser
parser = OptionParser(usage="python %prog -f filename [other options]", version="%prog 5")
parser.add_option("-t", "--turn",
                help="Turn number to start. Default is first turn: 1",
                default="1",dest="startTurn", type="int")
parser.add_option("-m", "--maxturns",
                help="Maximum number of turns to be analysed. Default is a number that is lower than the maximum which can be handled by drive: 9.500",
                default="9500", dest="maxturns", type="int")

parser.add_option("-v", "--sing_val",
                help="Keep this amount of singular values in decreasing order, rest will be cutted (set to 0). Default is a large number: 100.000",
                default="100000", dest="sing", type="int")

parser.add_option("-p", "--pk-2-pk",
                help="Peak to peak amplitude cut. This removes BPMs where abs(max(turn values) - min(turn values)) <= threshold. Default is 10e-8",
                default="0.00000001", dest="peak", type="float")
parser.add_option("-s", "--sumsquare",
                help="Threshold for single BPM dominating a mode. Should be > 0.9 for LHC",
                default="0.925", dest="sum", type="float")

parser.add_option("-f", "--file",
                help="File to clean",
                default="./", dest="file")
parser.add_option("-n", "--newfile",
                help="File name for new file. Default is override current file",
                default="", dest="newfile")


(options, args) = parser.parse_args()


#----INPUT VARIABLES
startTurnHuman = options.startTurn      # startTurn number to start
startTurn = startTurnHuman - 1          # startTurn number to start for internal usage. indices start with 0, not with 1
pk_pk_cut = options.peak                # peak to peak amplitude 
sumsquare = options.sum                 # sqroot(sumsquare of the 4 bpm values) should be > 0.9
sing_val = options.sing                 # keep svdcut number of singular values in decreasing order
maxTurnsHuman = options.maxturns        # maximum number of turns which should be parsed
maxTurns = maxTurnsHuman - 1            # maximum number of turns which should be parsed. indices start with 0, not with 1
newFile = options.newfile               # new file for SVDClean output
if newFile == "":                       # set to input file if no other output is given
    newFile = options.file
    print "no newfile given, existing file will be replaced with output"
#----------------------------------------
#-----Input variables check
if len(options.file) <= 2:
    print "Missing file name"
    sys.exit(1)
#----------------------------------------

# internal options
printTimes = False  # If True, the execution >>Time for each part of the script and overall time are printed
printDebug = False  # If True, internal debug information will be printed

# internal constants
planeX = "0"
planeY = "1"

class SddsFile(object):
    """This class represents a SDDS file.
    It is used to read in an existing SDDS file and to create and write the new one"""
    def __init__(self, pathToSddsFile):
        self.pathToSddsFile = pathToSddsFile
        self.badBpmFile = BadBpmFile(pathToSddsFile)
        self.parsed = False
        self.header = ""
        self.numberOfTurns = 0
        self.dictionaryPlaneToBpms = { planeX:Bpms(), planeY:Bpms() }
            
    def init(self):
        """Parse the current SDDS file and sets all member variables"""
        if (self.parsed):
            return
        timeStart = time.time()
        
        lastNumberOfTurns = 0
        detectedNumberOfTurns = 0
        flatBpmCounter = 0
        
        fileSdds = open(self.pathToSddsFile, "r")
        print "Extracting data from file..."
        for fileSddsItem in fileSdds: # Iterator over all lines (tbach)
            if fileSddsItem.startswith("#"): # then it is a comment line (tbach)
                self.header += fileSddsItem
                continue
            
            # from here, we have a data line (tbach)
            listSplittedLineValues = fileSddsItem.split()
            numberOfColumns = len(listSplittedLineValues)
            if (numberOfColumns < 3): # this is not a valid data line (tbach)
                continue
            plane = listSplittedLineValues[0]
            bpmName = listSplittedLineValues[1]
            location = float(listSplittedLineValues[2])
            
            # first 3 entries metadata, rest should be turn data (tbach)
            detectedNumberOfTurns = numberOfColumns - 3
            
            if lastNumberOfTurns > 0:
                if lastNumberOfTurns != detectedNumberOfTurns:
                    print "(plane", plane, ") BPM has a different number of turns then previous. BPM:", bpmName, "previous turns:", lastNumberOfTurns, "current turns:", detectedNumberOfTurns
                    sys.exit(1)
            else: # this will happen only once, when lastNumberOfTurns <= 0 (or multiple times if we do not have any turn data) (tbach)
                if (startTurnHuman > detectedNumberOfTurns):
                    print "startTurn > detectedNumberOfTurns. startTurn:", startTurnHuman, "detectedNumberOfTurns:", detectedNumberOfTurns
                    sys.exit(1)
            lastNumberOfTurns = detectedNumberOfTurns
            
            bpmsNameLocationPlane = (bpmName, location, plane)
            # be careful with the bpmsNameLocationPlane format, it is important for other parts of the program
            # sadly, it is a lot slower if we create an object for every BPM here instead of using the list
            # -- tbach
            
            self.numberOfTurns = min(maxTurns, detectedNumberOfTurns) - startTurn
            ndarrayLineData = numpy.array(listSplittedLineValues[3 + startTurn:3 + startTurn + self.numberOfTurns], dtype=numpy.float64)
            # this block handles BPMs with the same values for all turns (tbach)
            peakToPeakDifference = numpy.abs(numpy.max(ndarrayLineData) - numpy.min(ndarrayLineData))
            if peakToPeakDifference <= pk_pk_cut: # then do not use this BPM (tbach)
                reasonForBadBpm = "Flat BPM, the difference between all values is smaller than " + str(pk_pk_cut)
                badBpm = BadBpm(bpmsNameLocationPlane, ndarrayLineData, reasonForBadBpm)
                self.badBpmFile.addBadBpm(badBpm)
                flatBpmCounter += 1
                continue
            
            self.dictionaryPlaneToBpms[plane].bpmData.append(ndarrayLineData)
            self.dictionaryPlaneToBpms[plane].bpmsNameLocationPlane.append(bpmsNameLocationPlane)
        
        fileSdds.close()
        self.parsed = True
        if printTimes: print ">>Time for init (read in file):", time.time() - timeStart, "s"
        if flatBpmCounter > 0 : print "Flat BPMs detected. Number of BPMs removed:", flatBpmCounter
        print "Startturn:", startTurnHuman, "Maxturns:", maxTurnsHuman
        print "Number of turns:", self.numberOfTurns
        print "Horizontal BPMs:", self.dictionaryPlaneToBpms[planeX].getNumberOfBpms(), 
        print "Vertical BPMs:", self.dictionaryPlaneToBpms[planeY].getNumberOfBpms()
        print "Reading done"
                
    
    def removeBpmsWithIndicesForPlaneWithReason(self, badBpmIndices, plane, reason):
        badBpms = self.dictionaryPlaneToBpms[plane].removeBadBpmAndGetListGoodBpmIndicesAndListBadBpm(badBpmIndices, reason)
        self.badBpmFile.addBadBpms(badBpms)
    
    def writeFile(self):
        self.badBpmFile.writeBadBpms()
        print "writing data to temp file: " + self.pathToSddsFile + ".tmp_svd_clean"
        timeStart = time.time()
        fileClean = open(self.pathToSddsFile + ".tmp_svd_clean", "w")
        fileClean.write(self.header)
        if "NTURNS calculated" not in self.header:
            fileClean.write("#NTURNS calculated: " + str(self.numberOfTurns) + "\n")
        #fileClean.write("#Modified: {0} By: {1}\n".format(time.strftime("%Y-%m-%d#%H:%M:%S"), __file__)) # to get only the filename: os.path.basename(__file__)
        # The previous line requires at least python 2.6 (tbach)
        fileClean.write("#Modified: %s By: %s\n" % (time.strftime("%Y-%m-%d#%H:%M:%S"), __file__)) # to get only the filename: os.path.basename(__file__)
        self.writeBpmDataToFile(planeX, fileClean)
        self.writeBpmDataToFile(planeY, fileClean)
        fileClean.close()
        if printTimes: print ">>Time for writeFile:", time.time() - timeStart, "s"
        
        print "writing done. rename to:  ", newFile
        os.rename(options.file + ".tmp_svd_clean", newFile)
        
    def writeBpmDataToFile(self, plane, fileToWrite):
        currentBpms = self.dictionaryPlaneToBpms[plane]
        for i, bpmsNameLocationPlaneItem in enumerate(currentBpms.bpmsNameLocationPlane):
            # fileToWrite.write("{0} {1[0]:<15} {1[1]:>12.5f} ".format(plane, bpmsNameLocationPlaneItem)) # This line requires at least python 2.6 (tbach)
            # {0} is the first argument, {1[0]} from the second argument the first entry, {1[1] from the second argument the second entry (tbach)
            # :>15 means right aligned, filled up to 15 characters. example: "   BPMYA.4L1.B1" (tbach)
            # :>12.5 means right aligned float, filled up to 12 characters and fixed precision of 5. Example: " 23347.14262" (tbach)
            fileToWrite.write("%s %15s %12.5f " % (plane, bpmsNameLocationPlaneItem[0], bpmsNameLocationPlaneItem[1]))
            currentBpms.bpmData[i].tofile(fileToWrite, sep=" ", format="% 8.5f")
            # % 8.5f is a float, filled up to 8 characters and fixed precision of 5. if negative, preceded by a sign, if positive by a space (tbach)
            # Example1: " -0.94424", example2: "  1.25630" (tbach)
            fileToWrite.write("\n")  


class BadBpmFile(object):
    def __init__(self, pathToSddsFile):
        self.pathToSddsFile = pathToSddsFile
        self.pathToBadBpmFile = pathToSddsFile.replace(".new", "") + ".bad"
        self.linesToWrite = []
        
    def __writeHeader(self, fileHandleBadBpm):
        fileHandleBadBpm.write("@  FILE %s " + self.pathToSddsFile + "\n")
        fileHandleBadBpm.write("*  Plane NAME   S     TurnData\n")
        fileHandleBadBpm.write("$  %s    %s     %le   %le\n")
        
    def addBadBpm(self, badBpm):
        self.linesToWrite.append("%s %15s %12.5f %s \n# %s\n" % \
                                 (badBpm.plane, badBpm.name, badBpm.location, " ".join(["%12.5f" % x for x in badBpm.data]), badBpm.reason))
        # for explanation of string formatting, see writeBpmDataToFile (tbach)
        # requires at least python 2.6:
        # "{0} {1:<15} {2:>12.5f} {3} \n# {4}\n".format(badBpm.plane, badBpm.name, badBpm.location, " ".join(["{0:>10.5f}".format(x) for x in badBpm.data])
        # --tbach
    
    def addBadBpms(self, listBadBpm):
        for badBpm in listBadBpm:
            self.addBadBpm(badBpm)
        
    def writeBadBpms(self):
        timeStart = time.time()
        print "Create file for bad BPMs: ", self.pathToBadBpmFile
        fileBadBpms = open(self.pathToBadBpmFile, "w")
        self.__writeHeader(fileBadBpms)
        for linesToWriteItem in self.linesToWrite:
            fileBadBpms.write(linesToWriteItem)
        fileBadBpms.close()
        if printTimes: print ">>Time for writeBadBpms:", time.time() - timeStart, "s"

class Bpms(object):
    """This represents some BPMs. It is used to distinguish between BPMs for X and Y"""
    def __init__(self):
        self.bpmData = []
        self.bpmsNameLocationPlane = []
    
    def getNumberOfBpms(self):
        return len(self.bpmsNameLocationPlane)
    
    def removeBadBpmAndGetListGoodBpmIndicesAndListBadBpm(self, badBpmIndices, reason):
        badBpms = []
        goodBpmIndices = range(self.getNumberOfBpms())
        for index in badBpmIndices:
            badBpms.append(BadBpm(self.bpmsNameLocationPlane[index], self.bpmData[index], reason))
            goodBpmIndices.remove(index)
            
        # we do not have to update bpmData, because we do not read them again (tbach)
        
        # And do not forget to update the bpmsNameLocationPlane all with good BPM indices (tbach)
        self.bpmsNameLocationPlane = [self.bpmsNameLocationPlane[x] for x in goodBpmIndices]
            
        return badBpms



class BadBpm(object):
    def __init__(self, bpmsNameLocationPlane, data, reason):
        self.bpmsNameLocationPlane = bpmsNameLocationPlane
        self.data = data
        self.reason = reason # The reason why this is a bad bpm (tbach)
        
    @property
    def name(self):
        return self.bpmsNameLocationPlane[0]
        
    @property
    def location(self):
        return self.bpmsNameLocationPlane[1]
        
    @property
    def plane(self):
        return self.bpmsNameLocationPlane[2]


        
class svdHandler(object):

    def __init__(self):
        self.sddsFile = SddsFile(options.file)
        self.sddsFile.init()
        print ""
        
        timeStart = time.time()
        self.doSvdClean(planeX)
        self.doSvdClean(planeY)
        if printTimes: print ">>Time for svdClean (all):", time.time() - timeStart, "s"
        print ""
        
        self.sddsFile.writeFile()
        print ""
    
    def doSvdClean(self, plane):
        timeStart = time.time()
        print "(plane", plane, ") removing noise floor with SVD"
        
        A = numpy.array(self.sddsFile.dictionaryPlaneToBpms[plane].bpmData)
        numberOfBpms = A.shape[0]
        
        if numberOfBpms <= 10:
            sys.exit("Number of bpms <= 10")
        
        # normalise the matrix
        sqrtNumberOfTurns = numpy.sqrt(A.shape[1])
        A_mean = numpy.mean(A)
        A = (A - A_mean) / sqrtNumberOfTurns
        
        if printTimes: print ">>Time for svdClean (before SVD call):", time.time() - timeStart, "s"
        USV = self.getSingularValueDecomposition(A)
        if printTimes: print ">>Time for svdClean (after SVD call): ", time.time() - timeStart, "s"
        
        # remove bad BPM by SVD (tbach)
        goodBpmIndices = self.getBadBpmIndices(USV, plane)
        USV = (USV[0][goodBpmIndices], USV[1], USV[2])
        numberOfBpms = len(goodBpmIndices)
        
        #----SVD cut for noise floor
        if sing_val < numberOfBpms:
            print "(plane", plane, ") svdcut:", sing_val
            USV[1][sing_val:] = 0
        else:
            print "requested more singular values than available"
        

        A = matrixmultiply(USV[0], matrixmultiply(numpy.diag(USV[1]), USV[2]))
        # A0 * (A1 * A2) should require less operations than (A0 * A1) * A2, because of the different sizes
        # A0 has (M, K), A1 has (K, K) and A2 has (K, N) with K=min(M,N)
        # Most of the time, the number of turns is greater then the number of BPM, so M > N
        # --tbach 
        A = (A * sqrtNumberOfTurns) + A_mean
        self.sddsFile.dictionaryPlaneToBpms[plane].bpmData = A
        if printTimes: print ">>Time for doSvdClean:", time.time() - timeStart, "s"
        
    
    def getSingularValueDecomposition(self, matrix):
        return numpy.linalg.svd(matrix, full_matrices=False) # full matrices do not have any interesting value for us (tbach)

        
    def getBadBpmIndices(self, USV, plane):
        timeStart = time.time()
        
        # A is [u,s,v] with u * np.diag(s) * v = original matrix (tbach)
        U_t_abs = numpy.transpose(abs(USV[0])) # This creates a view, which is nice and fast (tbach)
        # What happens here?
        # From the SDDS ASCII file, We have a BPM x Turns matrix. Let B = Number of BPM, T = Number of Turns, then matrix size is (B,T)
        # If we do the SVD, we have U,S,V. With U has the size (B,x) and V (x,T)
        # If we transpose, U^t has the size (x,B)
        # Now, we can look in every row, what is the maximum value for every BPM. If one BPM is dominating, we remove it as a bad BPM
        # --tbach
        
        badBpmIndices = set() # we do not want duplicates (tbach)
        
        for rowIndex in range(len(U_t_abs)):
            maxIndex = numpy.argmax(U_t_abs[rowIndex])
            maxValue = U_t_abs[rowIndex][maxIndex]
            if (maxValue > sumsquare):
                badBpmIndices.add(maxIndex)
            
        if printDebug: print "Bad BPM indices: ", badBpmIndices
        
        numberOfBadBpms = len(badBpmIndices)
        if numberOfBadBpms > 0:
            print "(plane", plane, ") Bad BPMs from SVD detected. Number of BPMs removed:", len(badBpmIndices) 

        # add the bad BPMs to the general list of bad BPMs (tbach)
        reasonForBadBpm = "Detected from SVD, single peak value is greater then " + str(sumsquare)
        self.sddsFile.removeBpmsWithIndicesForPlaneWithReason(badBpmIndices, plane, reasonForBadBpm)
        
        numberOfBpms = USV[0].shape[0]
        goodBpmIndices = range(numberOfBpms)
        for value in badBpmIndices:
            goodBpmIndices.remove(value)
        
        if printTimes: print ">>Time for removeBadBpms:", time.time() - timeStart, "s"
        return goodBpmIndices    


timeStartGlobal = time.time()
svdHandler()
print "Global Time:", time.time() - timeStartGlobal, "s"
