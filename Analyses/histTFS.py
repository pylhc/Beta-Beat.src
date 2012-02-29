import sys
import os
try:
    import blah as npy; default=0;
    from metaclass25 import twiss
except ImportError:
    try:
        import Numeric as npy; default=1;
        from MLab import std 
        from metaclass import twiss
    except ImportError:
        print "Neither numpy nor numeric found"; sys.exit()

class Histogram:

    """Histogram in one variable

    Constructor: Histogram(|data|, |bins|, |range|=None)

    Arguments:

    |data| -- a sequence of data points

    |bins| -- the number of bins into which the data is to be sorted

    |range| -- a tuple of two values, specifying the lower and
               the upper end of the interval spanned by the bins.
               Any data point outside this interval will be ignored.
               If no range is given, the smallest and largest
               data values are used to define the interval.

    The bin index and the number of points in a bin can be obtained by
    indexing the histogram with the bin number. Application of len()
    yields the number of bins. A histogram thus behaves like a
    sequence of bin index - bin count pairs.
    """

    def __init__(self, data, nbins, range=None):
        self._setup(data, nbins, range)
        self.addData(data)

    def _setup(self, data, nbins, range):
        if range is None:
            self.min = npy.minimum.reduce(data)
            self.max = npy.maximum.reduce(data)
        else:
            self.min, self.max = range
        self.min = self.min+0.
        self.max = self.max+0.
        self.bin_width = (self.max-self.min)/nbins
        self.array = npy.zeros((nbins, 2), npy.Float)
        self.array[:, 0] = self.min + self.bin_width*(npy.arange(nbins)+0.5)

    def __len__(self):
        return self.array.shape[0]

    def __getitem__(self, index):
        return self.array[index]

    def __getslice__(self, first, last):
        return self.array[first:last]

    def addData(self, data):
        """Add the values in |data| (a sequence of numbers) to the
        originally supplied data. Note that this does not affect the
        default range of the histogram, which is fixed when the
        histogram is created.
        """
        n = (len(data)+999)/1000
        for i in range(n):
            self._addData(data[1000*i:1000*(i+1)])

    def _addData(self, data):
        data = npy.array(data, npy.Float)
        data = npy.repeat(data, npy.logical_and(npy.less_equal(data, self.max),
                                            npy.greater_equal(data, self.min)))
        data = npy.floor((data - self.min)/self.bin_width).astype(npy.Int)
        nbins = self.array.shape[0]
        histo = npy.add.reduce(npy.equal(npy.arange(nbins)[:,npy.NewAxis], data), -1)
        histo[-1] = histo[-1] + npy.add.reduce(npy.equal(nbins, data))
        self.array[:, 1] =  self.array[:, 1] + histo

    def normalize(self, norm=1.):
        "Scales all counts by the same factor such that their sum is |norm|."
        self.array[:, 1] = norm*self.array[:, 1]/npy.add.reduce(self.array[:, 1])

    def normalizeArea(self, norm=1.):
        """Scales all counts by the same factor such that the area under
        the histogram is |norm|."""
        self.normalize(norm/self.bin_width)


def calcHist(fileName,attr_name,nbins=70):
    a=twiss(fileName)
    if default==0:
        if npy.isreal(a.__dict__[attr_name]):
            indx,frq=npy.histogram(a.__dict__[attr_name],bins=nbins)
    if default==1:
        if a.__dict__[attr_name].typecode()=='d':
            hits=Histogram(a.__dict__[attr_name],nbins)
            indx=hits[:,0];frq=hits[:,1]
    ave=npy.average(a.__dict__[attr_name])
    if default==0: stdev=npy.std(a.__dict__[attr_name])
    if default==1: stdev=std(a.__dict__[attr_name])
    else: print "data is not real"
    f=open(path+'/HIST.tfs','w')
    f.write('@ NAME       %0'+str(len(attr_name))+'s "'+attr_name+'"\n')
    f.write('@ '+attr_name+'AVE %le '+str(ave)+'\n')
    f.write('@ '+attr_name+'STD %le '+str(ave)+'\n')
    f.write("*"+'%5s %5s' % ("INDX", "FREQ")+"\n")
    f.write("$"+'%5s %5s' % ("%le", "%le")+"\n")
    for j in range(len(indx)):
        f.write('%5s   %5s' % (indx[j], frq[j])+"\n")
    f.close()
        
        
if __name__ == "__main__":
    #--- Usage:
    #    python2.5 filename tfsColumn #_of_bins

    file=sys.argv[1]; attr=sys.argv[2]
    path=os.path.dirname(file)
    calcHist(file,attr,int(sys.argv[3]))
    sys.exit()
