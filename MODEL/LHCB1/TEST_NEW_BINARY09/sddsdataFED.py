#ATTEMPT EDIT, F.R. 30-06-09:
# -add header interpretation for getting parameter names
 
#########################################################
#TODO pages not supported
#TODO tables not supported
#TODO multi dimensional array in binary are not supported


# struct format
#     >: big
#     <: little
#     |: machine
#     x: pad byte (no data);
#     c:char;
#     b:signed byte;
#     B:unsigned byte;
#     h:short;
#     H:unsigned short;
#     i:int;
#     I:unsigned int;
#     l:long;
#     L:unsigned long;
#     f:float;
#     d:double.
#     s:string (array of char);
#     p:pascal string (with count byte).
#     P:an integer type that is wide enough to hold a pointer.
#     q:long long;
#     Q:unsigned long long

import StringIO
import gzip
import numpy as n
import struct


def myopen(fn):
  try:
    if fn.endswith('.gz'):
      return gzip.open(fn)
    else:
      return open(fn)
  except IOError:
    return StringIO.StringIO(fn)



def iterheader(o):
  c=1
  while c:
    c=o.read(1)
    if c=='!':
      while c not in '\n\r':
        c=o.read(1)
    elif c not in '\n\t\r':
       yield c

def iterbinarydata(o):
  c=1
  while c:
    c=o.read(1)
    if c=='!':
      while c not in '\n\r':
        c=o.read(1)
      if c in '\n\r': # possibly a bug
        c=o.read(1)
    else:
      yield c

def readtoken(o):
  buf=[]
  for i in o:
    buf.append(i)
    if ''.join(buf[-4:])=='&end':
      yield ''.join(buf)
      buf=[]

def myreadline(o):
  buf=o.readline()
  while buf[0]=='!':
    buf=o.readline()
  return buf

def parseheader(l):
  t,data=l.split(' ',1)
  data=data.replace(' ','')
  data=[d.split('=') for d in data.split(',')]
  data.pop()
  data=dict(data)
  data['header']=t
  return data



sddstypes={
    'short' : 'i2',
    'long'  : 'i4',
    'llong' : 'u8',
    'string': 'S',
    'float': 'f4',
    'double': 'f8',
}

def myarray(fh,typ,count,endian):
#  print typ
  typ=n.dtype(endian+sddstypes.get(typ,typ))
  size=typ.itemsize*count
  ss=fh.read(size)
  if len(ss)==size:
#  print typ,count,size,repr(ss[:16])
    return n.fromstring(ss,dtype=typ,count=count)
  else:
    return None
#  return s

def mystruct(fh,typ,count,endian):
  typ='%s%d%s' % (endian,count,typ)
  size=struct.calcsize(typ)
  ss=fh.read(size)
  if len(ss)==size:
    return struct.unpack(typ,ss)
  else:
    return None

def mysplit(fh,count):
  out=[]
  while len(out)<count:
   l=fh.readline()
   out.extend(l.split())
  return out


class sddsdata(object):
  def __init__(self,filename,endian='little',full=True):
    self.endian={'little':'<','big':'>'}[endian]
    self.filename=filename
    self.fh=myopen(filename)
    self.version=self.fh.readline()
    # read headear
#    print 'read header'
    it=readtoken(iterheader(self.fh))
    header=[]
    for i in it:
     header.append(parseheader(i))
     if header[-1]['header']=='&data':
       break
    header2=[]; istable=True
    for i in header:
      if i['header']=='&column':
        if istable==True:
          header2.append({'header':'&table'})
          header2[-1]['columns']=[i]
          istable=False
        else:
          header2[-1]['columns'].append(i)
      else:
        header2.append(i)
    self.header=header2
#    print self.header
    # read data
    if full:
      self.fh.read(1)
      self.data=[]
      if self.header[-1]['mode']=='ascii':
        data={}
        self.data.append(data)
        for i in self.header:
          if 'type' in i:
            typ=i['type']
            typ=n.dtype(sddstypes.get(typ,typ))
            if i['header']=='&parameter':
              ss=myreadline(self.fh)
              print ss, typ, i['type']
              #Before d=n.array(ss,typ)
              d=n.array(ss)
              print d, "again"
              i['value']=d
            elif i['header']=='&array':
              dims=map(int,myreadline(self.fh).split())
              i['shape']=dims
              cnt=reduce(lambda a,b:a*b,dims)
  #            ss=myreadline(self.fh)
  #            print dims, len(ss)
  #Before d=n.array(mysplit(self.fh,cnt),typ).reshape(dims)
              d=n.array(mysplit(self.fh,cnt)).reshape(dims)
            data[i['name']]=d
      elif self.header[-1]['mode']=='binary':
        while 1:
          row=myarray(self.fh,'long',1,self.endian)
          if row is None:
            break
          data={}
          self.data.append(data)
          for i in self.header:
            if 'type' in i:
              typ=i['type']
              if i['header']=='&parameter':
                if typ=='string':
                  smax=0
                  subcount=myarray(self.fh,'long',1,self.endian)[0]
                  smax=subcount<smax and smax or subcount
                  d=mystruct(self.fh,'s',subcount,self.endian)[0]
                else:
                  d=myarray(self.fh,typ,1,self.endian)
              elif i['header']=='&array':
                count=myarray(self.fh,'long',1,self.endian)[0]
                if typ=='string':
                  d=[];smax=0
                  for r in range(count):
                    subcount=myarray(self.fh,'long',1,self.endian)[0]
                    smax=subcount<smax and smax or subcount
    #                d.append(myarray(self.fh,'>S1',subcount,self.endian)[0])
                    d.append(mystruct(self.fh,'s',subcount,self.endian)[0])
                  d=n.array(d,n.dtype('S%s'%smax))
                else:
                  d=myarray(self.fh,typ,count,self.endian)
              data[i['name']]=d
    self.fh.close()
  def __str__(self):
    out=['%s: %s' % (self.filename,self.version)]
    for i in self.header:
      oo=[]
      for n,k in i.items():
        if n is not 'header':
          oo.append('%s=%s' % (n,k))
      out.append(i['header'][1:]+' '+', '.join(oo))
    return '\n'.join(out)
  __repr__=__str__


 #### F.R. 30-06-09
  def treat_header(self):
    self.par_list = []
    lines = self.header
    

    #        lines = open(filename).read().split('\n')
    #    except IOError:            
    #        print 'cannot open file ' + filename
    #        return

    self.parameters = []
    self.records = []
        

    idx = 0
    nparams = 0
    narrays = 0

    for line in lines:

      if line['header']=='&parameter':
        name  = line['name']
        type  = line['type']
        value = line['value']
        print name,type,value
        if (type == 'long') | (type == 'double'):
          type = '%le'
        elif type =='int':
          type = '%d'
        elif type =='string':
          if len(value)<10:
            type = '%0' + str(len(value)) + s
          else:
            type = '%' + str(len(value)) + s

        self.par_list.append([name,value,type])
        
#        QUAAAAAAAAA da qua: preparare header per TFS table
'''
    #        words = line.replace('=',' ').replace(',',' ').split()

            if len(words)<1:
                continue        

            if line[0] == '!':
                continue

            if words[0] == '&parameter':
                #print 'adding parameter ' + words[2] + ' type ' + words[4]
                self.parameters.append([words[2], nparams])

                if words[2] == 'date':
                    date_id = nparams
                
                if words[2] == 'time':
                    time_id = nparams            
            
                idx = idx + 1
                nparams = nparams + 1
                continue

            if words[0] == '&array':
                #print 'adding array ' + words[2] + ' type ' + words[4]
                self.parameters.append([words[2], nparams + 1 + narrays*2])

                 
                idx = idx + 1
                narrays = narrays + 1
                continue

            if line[0] == '&':
                continue

            # start data block
            self.records.append(line.split())




    def get_parameter_names(self):
        parnames = []
        for p in  self.parameters:
            parnames.append(p[0])
        return parnames
'''
