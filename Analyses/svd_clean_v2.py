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
parser = OptionParser()
parser.add_option("-t", "--turn",
                help="Turn number to start",
                metavar="startTurn", default="1",dest="startTurn")
parser.add_option("-p", "--pk-2-pk",
                help="Peak to peak amplitude cut",
                metavar="peak", default="0.00000001", dest="peak")
parser.add_option("-s", "--sumsquare",
                help="threshold for single BPM dominating a mode. Should be > 0.9",
                metavar="sun", default="0.925", dest="sum")
parser.add_option("-v", "--sing_val",
                help="keep svdcut number of singular values in decreasing order",
                metavar="sing", default="100000", dest="sing")
parser.add_option("-f", "--file",
                help="file to clean",
                metavar="file", default="./", dest="file")
parser.add_option("-m", "--maxturns",
                help="Maximum number of turns to be analysed",
                metavar="MAXTURNS", default="4500", dest="maxturns")

(options, args) = parser.parse_args()


#----INPUT VARIABLES
startTurnHuman = int(options.startTurn) # startTurn number to start
startTurn = startTurnHuman - 1          # startTurn number to start for internal usage. indices start with 0, not with 1
pk_pk_cut = float(options.peak)         # peak to peak amplitude 
sumsquare = float(options.sum)          # sqroot(sumsquare of the 4 bpm values) should be > 0.9
sing_val = int(float(options.sing))     # keep svdcut number of singular values in decreasing order
maxturns = int(options.maxturns)        # maximum number of turns which should be parsed
#----------------------------------------
#-----Input variables check
if len(options.file) <= 2:
    print "Missing file name"
    sys.exit(1)
if startTurn >= maxturns :
    print "startTurn >= maxturns", startTurn, ">=", maxturns
    sys.exit(1)
#----------------------------------------

# internal options
printTimes = True  # If True, the execution >>Time for each part of the script and overall time are printed
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
        self.bpmsX = Bpms(planeX)
        self.bpmsY = Bpms(planeY)
            
    def init(self):
        """Parse the current SDDS file and sets all member variables"""
        if (self.parsed):
            return
        timeStart = time.time()
        
        lastNumberOfTurns = 0
        
        fileSdds = open(self.pathToSddsFile, "r")
        print "Extracting data from file..."
        for fileSddsItem in fileSdds: # Iterator over all lines (tbach)
            if fileSddsItem.startswith("#"): # then it is a comment line (tbach)
                self.header += fileSddsItem
                continue
            
            # from here, we have a data line (tbach)
            listSplittedLineValues = fileSddsItem.split()
            if (len(listSplittedLineValues) < 3): # this is not a valid data line (tbach)
                continue
            plane = listSplittedLineValues.pop(0)
            bpmName = listSplittedLineValues.pop(0)
            location = float(listSplittedLineValues.pop(0))
            
            # pop removed the first 3 entries, rest should be turn data (tbach)
            self.numberOfTurns = min(maxturns, len(listSplittedLineValues))
            
            if lastNumberOfTurns > 0:
                if lastNumberOfTurns != self.numberOfTurns:
                    print "BPM has a different number of turns then previous. BPM:", bpmName, "previous turns:", lastNumberOfTurns, "current turns:", self.numberOfTurns
                    sys.exit(1)
            else: # this will happen only once, when lastNumberOfTurns <= 0 (tbach)
                if (startTurnHuman > self.numberOfTurns):
                    print "startTurn > numberOfTurns. startTurn:", startTurnHuman, "numberOfTurns:", self.numberOfTurns
                    sys.exit(1)
            lastNumberOfTurns = self.numberOfTurns
            
            currentBpmsForPlane = self.getBpmsForPlane(plane)
            currentBpmsForPlane.bpmData.append(numpy.array(listSplittedLineValues, dtype=numpy.float64))
            currentBpmsForPlane.bpmsNameLocationPlane.append((bpmName, location, plane))
            # be careful with the bpmsNameLocationPlane format, it is important for other parts of the program
            # sadly, it is a lot slower if we create an object for every BPM here instead of using the list
            # -- tbach
        
        fileSdds.close()
        self.parsed = True
        if printTimes: print ">>Time for init (read in file):", time.time() - timeStart, "s"
        print "Startturn:", startTurnHuman, "Maxturn:", maxturns
        print "Number of turns:", self.numberOfTurns
        print "Horizontal BPMs:", self.getBpmsForPlane(planeX).getNumberOfBpms(), 
        print "Vertical BPMs:", self.getBpmsForPlane(planeY).getNumberOfBpms(), "\n"
                
    
    def getBpmsForPlane(self, plane):
        if plane == planeX:
            return self.bpmsX
        elif plane == planeY:
            return self.bpmsY
        else:
            print "Unknown plane:", plane
            sys.exit(1)
    
    def removeBpmsWithIndicesForPlaneWithReason(self, badBpmIndices, plane, reason):
        badBpms = []
        if len(badBpmIndices) == 0:
            return badBpms
        
        currentBpms = self.getBpmsForPlane(plane)
        
        badBpms = currentBpms.getListBadBpm(badBpmIndices, reason)
        self.badBpmFile.writeBadBpms(badBpms)
        
        currentBpms.removeAllBpmsWithIndices(badBpmIndices)
        
    
    def writeFile(self):
        print "writing data to temporary file"
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
        
        print "writing done. replace original file"
        os.rename(options.file + ".tmp_svd_clean", options.file)
        
    def writeBpmDataToFile(self, plane, fileToWrite):
        currentBpms = self.getBpmsForPlane(plane)
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
        self.wroteHeader = False
        
    def writeHeader(self):
        if (self.wroteHeader): # Do it only once (tbach)
            return
        print "Create file for bad BPMs: ", self.pathToBadBpmFile
        fileBadBpms = open(self.pathToBadBpmFile, "w")
        fileBadBpms.write("@  FILE %s " + self.pathToSddsFile + "\n")
        fileBadBpms.write("*  Plane NAME   S     TurnData\n")
        fileBadBpms.write("$  %s    %s     %le   %le\n")
        fileBadBpms.close()
        self.wroteHeader = True
        
    def writeBadBpms(self, listBadBpm):
        if (len(listBadBpm) is 0):
            return
        if (not self.wroteHeader): self.writeHeader() # Only write the header and therefore the file, if needed (tbach)
        
        fileBadBpms = open(self.pathToBadBpmFile, "a")
        for badBpm in listBadBpm:
            # fileBadBpms.write("{0} {1:<15} {2:>12.5f} ".format(badBpm.plane, badBpm.name, badBpm.location)) # This line requires at least python 2.6
            fileBadBpms.write("%s %15s %12.5f " % (badBpm.plane, badBpm.name, badBpm.location))
#            for data in badBpm.data:
                #fileBadBpms.write("{0:>10.5f}".format(data)) # This line requires at least python 2.6
            fileBadBpms.write(" ".join(["%12.5f" % x for x in badBpm.data]))
            fileBadBpms.write("\n# " + badBpm.reason + "\n")
        fileBadBpms.close()
        
        

class Bpms(object):
    """This represents some BPMs. It is used to distinguish between BPMs for X and Y"""
    def __init__(self, plane):
        self.bpmData = []
        self.bpmsNameLocationPlane = []
        self.plane = plane
        self.numberOfBpms = 0
    
    def getNumberOfBpms(self):
        return len(self.bpmsNameLocationPlane)
    
    def removeAllBpmsWithIndices(self, badBpmIndices):
        if len(badBpmIndices) == 0:
            return
        
        # First, we need all good BPM indices (tbach)
        goodBpmIndices = range(self.getNumberOfBpms())
        for value in badBpmIndices:
            goodBpmIndices.remove(value)
        
        # Update data with all good BPM indices (tbach)
        self.bpmData = [self.bpmData[x] for x in goodBpmIndices]
        
        # And do not forget to update the bpmsNameLocationPlane all with good BPM indices (tbach)
        self.bpmsNameLocationPlane = [self.bpmsNameLocationPlane[x] for x in goodBpmIndices]
        
    def getListBadBpm(self, badBpmIndices, reason):
        badBpms = []
        for index in badBpmIndices:
            badBpms.append(BadBpm(self.bpmsNameLocationPlane[index], self.bpmData[index], reason))
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


        
class svdBadBpm(object):

    def __init__(self):
        self.sddsFile = SddsFile(options.file)
        self.sddsFile.init()
        
        self.removeBadBpms(planeX)
        print ""
        self.removeBadBpms(planeY)
        print ""
        self.svdClean(planeX)
        self.svdClean(planeY)
        
        self.sddsFile.writeFile()
    
    
    def removeBadBpms(self, plane):
        self.removeFlatBpms(plane)
        self.removeBadBpmsFromSVD(plane)
        
        
    def removeFlatBpms(self, plane):
        timeStart = time.time()
        bpms_x_TurnsMatrix = numpy.array(self.sddsFile.getBpmsForPlane(plane).bpmData)
        
        badBpmIndices = []
        for rowIndex in range(len(bpms_x_TurnsMatrix)):
            pk_pk = numpy.max(bpms_x_TurnsMatrix[rowIndex]) - numpy.min(bpms_x_TurnsMatrix[rowIndex])
            if pk_pk <= pk_pk_cut:
                badBpmIndices.append(rowIndex)
        
        if printDebug: print "plane, len(badBpmIndices), badBpmIndices: ", plane, len(badBpmIndices), badBpmIndices
        print "Number of bad BPMs from peak to peak cut:", len(badBpmIndices), "plane:", plane
        
        reasonForBadBpm = "Flat BPM, difference between values is smaller then " + str(pk_pk_cut)
        self.sddsFile.removeBpmsWithIndicesForPlaneWithReason(badBpmIndices, plane, reasonForBadBpm)
        
        if printTimes: print ">>Time for removeFlatBpms:", time.time() - timeStart, "s"
    
    def removeBadBpmsFromSVD(self, plane):
        timeStart = time.time()
        normalizedBpmDataMatrix = self.getNormalizedMatrix(plane)
        USV = self.getSingularValueDecomposition(normalizedBpmDataMatrix)
        
        # A is [u,s,v] with u * np.diag(s) * v = original matrix (tbach)
        U_t_abs = numpy.transpose(abs(USV[0])) # This creates a view, which is nice and fast (tbach)
        # What happens here?
        # From the SDDS ASCII file, We have a BPM x Turns matrix. Let B = Number of BPM, T = Number of Turns, then matrix size is (B,T)
        # If we do the SVD, we have U,S,V. With U has the size (B,x) and V (x,T)
        # If we transpose, U^t has the size (x,B)
        # Now, we can look in every row, what is the maximum value for every BPM. If one BPM is dominating, we remove it as a bad BPM
        # --tbach
        
        badBpmIndices = set()
        
        for rowIndex in range(len(U_t_abs)):
            maxIndex = numpy.argmax(U_t_abs[rowIndex])
            maxValue = U_t_abs[rowIndex][maxIndex]
            if (maxValue > sumsquare):
                badBpmIndices.add(maxIndex)
            
        if printDebug: print "Bad BPM indices: ", badBpmIndices
        
        print "Number of bad BPMs from SVD:", len(badBpmIndices), "plane:", plane 

        reasonForBadBpm = "Detected from SVD, single peak value is greater then " + str(sumsquare)
        self.sddsFile.removeBpmsWithIndicesForPlaneWithReason(badBpmIndices, plane, reasonForBadBpm)
        
        if printTimes: print ">>Time for removeBadBpms:", time.time() - timeStart, "s"


    def getNormalizedMatrix(self, plane):
        timeStart = time.time()
        B = numpy.array(self.sddsFile.getBpmsForPlane(plane).bpmData)
        numberOfBpms = B.shape[0]
        numberOfTurns = B.shape[1]
        
        if numberOfTurns <= 0 or numberOfBpms <= 10:
            print "No turns or BPMs left. Turns:", numberOfTurns, "Number of BPMs:", numberOfBpms
            sys.exit()
        B = (B - numpy.mean(B)) / numpy.sqrt(numberOfTurns)
        if printTimes: print ">>Time for getNormalizedMatrix:", time.time() - timeStart, "s"
        return B
    
       
    def getSingularValueDecomposition(self, matrix):
        return numpy.linalg.svd(matrix, full_matrices=False) #singular_value_decomposition
    
            
    def svdClean(self, plane):
        timeStart = time.time()
        print "removing noise floor for plane:", plane
        
        B = numpy.array(self.sddsFile.getBpmsForPlane(plane).bpmData)
        numberOfBpms = B.shape[0]
        
        if numberOfBpms <= 10:
            sys.exit("Number of bpms <= 10")
            
        sqrtNumberOfTurns = numpy.sqrt(B.shape[1])
        B_mean = numpy.mean(B)
        B = (B - B_mean) / sqrtNumberOfTurns
        
        if printTimes: print ">>Time for svdClean1:", time.time() - timeStart, "s"
        USV = self.getSingularValueDecomposition(B)
        if printTimes: print ">>Time for svdClean2:", time.time() - timeStart, "s"
        
        
        #----SVD cut for noise floor
        svdcut = sing_val
        if svdcut > numberOfBpms:
            svdcut = numberOfBpms
            print "requested more singular values than available"
            
        print "svdcut:", svdcut, "for plane:", plane
        USV[1][svdcut:] = 0

        B = matrixmultiply(USV[0], matrixmultiply(numpy.diag(USV[1]), USV[2]))
        # A0 * (A1 * A2) should require less operations than (A0 * A1) * A2, because of the different sizes
        # A0 has (M, K), A1 has (K, K) and A2 has (K, N) with K=min(M,N)
        # Most of the time, the number of turns is greater then the number of BPM, so M > N
        # --tbach 
        
        B = (B * sqrtNumberOfTurns) + B_mean
        self.sddsFile.getBpmsForPlane(plane).bpmData = B
        if printTimes: print ">>Time for svdClean:", time.time() - timeStart, "s"
        
timeStartGlobal = time.time()
svdBadBpm()
print "Global Time:", time.time() - timeStartGlobal, "s"
