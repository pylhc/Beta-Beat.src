#########################
#
#
# IO functions
# online model
#
#########################

import gzip

#########################
# sdds 
#########################


# reader 

class SDDSReader:

    def __init__(self, filename):

        try:
            lines = open(filename).read().split('\n')
        except IOError:            
            print 'cannot open file ' + filename
            return

        self.parameters = []
        self.records = []
        

        idx = 0
        nparams = 0
        narrays = 0

        for line in lines:

            if line == 'SDDS1':
                continue


            words = line.replace('=',' ').replace(',',' ').split()

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
        


    def get_parameter(self,name):

        idx = -1

        for p in self.parameters:
            if p[0] == name:
                idx = p[1]

        if idx >= 0:
            return self.records[idx]




#########################
# tfs 
#########################


# reader

class TFSReader:

    def __init__(self, filename):

        try:
            if filename.endswith('.gz'):
                lines = gzip.open(filename).read().split('\n')
            else:
                lines = open(filename).read().split('\n')
        except IOError:            
            print 'cannot open file ' + filename
            return

        self.scalar_parameters = []
        self.vector_parameters = [] # name, idx
        self.vector_parameter_types = [] # name, idx
        self.records = []

        self.filters = {}

        self.row_count = 0


        for line in lines:
            words = line.replace('\"','').split()

            if len(words)<1:
                continue

            # TODO : not correct when attributes are strings with whitespaces
            if words[0] == '@':
                if len(words) >= 4:
                    self.scalar_parameters.append([words[1],words[3]])
                    #print 'adding scalar ' + words[1] 
                else:                
                    print 'WARNING: parse error : ' + line + ' ' + str(len(words))

                continue


            if words[0] == '*':
                vector_par_id = 0
                for p in words[1:]:
                    self.vector_parameters.append([p,vector_par_id])
                    vector_par_id = vector_par_id + 1
                continue
                    
                    
            if words[0] == '$':
                for p in words[1:]:
                    self.vector_parameter_types.append(p)
                continue


            self.records.append(words)


    def rewind(self):
        self.row_count = 0

    def pop_record(self):

        if self.row_count >= len(self.records):
            return []
        else:
            self.row_count = self.row_count + 1
            return self.records[self.row_count - 1]


    def get_parameter_names(self):
        names = []

        for p in self.scalar_parameters:
            names.append(p[0])

        for p in self.vector_parameters:
            names.append(p[0])


        return names

    
    def get_vector_parameter_names(self):
        names = []

        for p in self.vector_parameters:
            names.append(p[0])


        return names

    def get_vector_parameter_types(self):
        return self.vector_parameter_types

    
    def get_scalar_parameter_names(self):
        names = []

        for p in self.scalar_parameters:
            names.append(p[0])


        return names

    def get_parameter_id(self,name):
        idx = -1
        for i in xrange(0,len(self.vector_parameters)):
            if self.vector_parameters[i][0].lower().strip() == name.lower():
                idx = self.vector_parameters[i][1]


        return idx

    def get_parameter_ids(self,names):
        ids = []

        for n in names:
            id = self.get_parameter_id(n)
            if id >= 0:
                ids.append(id)

        return ids


    def get_scalar_parameter(self,name):

        for p in self.scalar_parameters:
            if p[0].lower() == name.lower():
                return p[1]
            

    def get_parameter(self,name):

        res = []
        
        idx = self.get_parameter_id(name)
        f_idx  = {}

        for f in self.filters.keys():
            i = self.get_parameter_id(f)
            f_idx[i] = self.filters[f]

        if idx >=0:
            for r in self.records:
                doAppend = True
                for i in f_idx.keys():
                    if r[i] not in f_idx[i]:
                        doAppend = False
                if doAppend:
                    res.append(r[idx])


        return res
    

        

# writer

class TFSWriter:

    def __init__(self):

        self.scalar_parameters = []
        self.vector_parameters = [] # name, idx, type

        self.records = []

        vector_par_id = 0


    def set_parameters(self, names, types):

        idx = 0

        for i in xrange(0,len(names)):
            self.vector_parameters.append([names[i], idx, types[i]])
            idx = idx + 1


    def add_scalar_parameter_value(self,par,value,type):
        self.scalar_parameters.append([par,value,type])


    def get_parameter_id(self,name):
        idx = -1
        for i in xrange(0,len(self.vector_parameters)):
            if self.vector_parameters[i][0].lower().strip() == name.lower():
                idx = self.vector_parameters[i][1]


        return idx


    def add_record(self,record):
        if len(record) != len(self.vector_parameters):
            print 'add row: length error'
            print record
            print self.vector_parameters
            return

        self.records.append(record)

    def write(self, filename):

        of = open(filename,'w')

        of.write('@ NAME             %05s \"user\"\n' + \
                     '@ TYPE             %05s \"user\"\n' + \
                     '@ TITLE            %09s \"gradients\"\n')

        for i in xrange(0,len(self.scalar_parameters)):
            p = self.scalar_parameters[i]
            of.write('@ ' + str(p[0]) + ' ' + str(p[2]) + ' ' + str(p[1]) + '\n')

        of.write('* ')

        for i in xrange(0,len(self.vector_parameters)):
            of.write(' ' + self.vector_parameters[i][0])
        of.write('\n')

        of.write('$ ')

        for i in xrange(0,len(self.vector_parameters)):
            of.write(' ' + self.vector_parameters[i][2])
        of.write('\n')

        for r in self.records:
            for i in xrange(0,len(r)):
                if self.vector_parameters[i][2].find('s')>=0:
                    of.write('\"' + str(r[i]) + '\" ' )
                else:
                    of.write(str(r[i]) + ' ' )
            of.write('\n')

            
#########################
# misc specific file format converters
#########################

# merge tfs files into one; 'keys' columns for both file should coincide
# so far the header from the first file is kept
def merge_tfs_files(f1, f2, f, keys):
    
    r1 = TFSReader(f1)
    r2 = TFSReader(f2)

    pars1 = r1.get_scalar_parameter_names()
    pars2 = r2.get_scalar_parameter_names()    

    w = TFSWriter() 

    for i in xrange(0,len(pars1)):
        name = pars1[i]
        val = r1.get_scalar_parameter(name)

        #print name + '  ----> ' + str(val)

        try:
            eval(val.replace('"',''))
            dval = eval(val)
            ptype = '%le'
            w.add_scalar_parameter_value(name,dval,ptype)
        except:
            #print 'cannot evaluate ' + str(val)
            ptype = '%s'
            w.add_scalar_parameter_value(name,val,ptype)



    # assuming that we have x xp in f1 and y yp in y yp

    #id_s_1 = r1.get_parameter_id('s')
    #id_name_1 = r1.get_parameter_id('name')

    #id_s_2 = r2.get_parameter_id('s')
    #id_name_2 = r2.get_parameter_id('name')

    #id_x = r1.get_parameter_id('x')
    #id_xp = r1.get_parameter_id('xp')

    #id_y = r2.get_parameter_id('y')
    #id_yp = r2.get_parameter_id('yp')

    
    pars1 = r1.get_vector_parameter_names()
    typs1 = r1.get_vector_parameter_types()

    pars2 = r2.get_vector_parameter_names()
    typs2 = r2.get_vector_parameter_types()


    if not is_subset(keys,pars1):
        print 'keys are not contained in table 1'
        print keys
        print pars1
        return 

    if not is_subset(keys,pars2):
        print 'keys are not contained in table 2'
        print keys
        print pars2
        return

    pars1 = exclude(pars1,keys)
    pars2 = exclude(pars2,keys)
    pars2 = exclude(pars2,pars1)
    pars = keys + pars1+pars2

    key_ids1 = r1.get_parameter_ids(keys)
    key_ids2 = r2.get_parameter_ids(keys)
    par_ids1 = r1.get_parameter_ids(pars1)
    par_ids2 = r2.get_parameter_ids(pars2)

    #TODO: add actual types
    typs = []
    for i in xrange(0,len(pars)):
        typs.append('%le')



    #if id_s_1 < 0 or id_s_2 < 0 or id_name_1 < 0 or id_name_2 < 0 or id_x < 0 or id_xp < 0 or id_y < 0 or id_yp < 0:
    #    print 'WARNING: not enough parameters for complete orbit table'
    #    return 


    w.set_parameters(pars, typs)
  
    while True:
        rc1 = r1.pop_record()
        rc2 = r2.pop_record()


        if len(rc1) < 1 or len(rc2)<1:
            break

        #s1 = eval(rc[id_s_1])
        #s2 = eval(rc[id_s_2])

        #name1 = eval(rc[id_name_1])
        #name2 = eval(rc[id_name_2])

        kvals1 = []

        for k1 in key_ids1:
            kvals1.append(rc1[k1].replace('\"',''))

        kvals2 = []

        for k2 in key_ids2:
            kvals2.append(rc2[k2].replace('\"',''))


        if not kvals1 == kvals2:
            print 'rows do not match, quitting'
            print kvals1
            print kvals2
            return
        

        ps = []

        for i in par_ids1:
            ps.append(rc1[i])

        for i in par_ids2:
            ps.append(rc2[i])
        

        rc = kvals1+ps

        #print ' adding ' + str(ps)
        
        w.add_record(rc)
            

    print 'saving file...'
    w.write(f)







def tfs2yasp(old_file, new_file):
    # yasp format
    #             * NAME  PLANE  BEAM  POS  RMS  SUM  HW-STATUS STATUS STATUS-TAG
    #             * NAME  PLANE  BEAM  STRENGTH-NAME  KICK
    # mad units - rad

    lines = open(old_file).read().split('\n')

    f = open(new_file,"w")

    f.write("@ YASP %s V2.0 \n")
    f.write("# CORRECTOR %s microrad \n")
    f.write("* NAME  PLANE  BEAM  STRENGTH-NAME  KICK \n")

    plane = 'H'
    beam = '1'

    for line in lines:
        words = line.split()

        if len(words) < 3:
            continue

        if words[0] in ['@','#','$','*']:
            continue

        px_cor = words[3]
        py_cor = words[4]

        str_name = mad2dev(words[0])
        name = words[0].replace('"','')

        
        if abs(eval(px_cor)) > abs(eval(py_cor)):
            plane = 'H'
            f.write(name + '  ' + plane + '  ' + beam + '  ' + str_name + '  '+   '%(pxc)22f \n' % {'pxc': eval(px_cor)*1.e+6 })
        if abs(eval(px_cor)) < abs(eval(py_cor)):
            plane = 'V'
            f.write(name + '  ' + plane + '  ' + beam + '  ' + str_name + '  '+  '%(pyc)22f \n' % {'pyc': eval(py_cor)*1.e+6 })


def yasp2tfs(file, outf_x, outf_y):

    try:
        scale = 1

        if file.endswith('.gz'):
            lines = gzip.open(file).read().split('\n')
        else:
            lines = open(file).read().split('\n')

        fx = open(outf_x,'w')
        fy = open(outf_y,'w')
        
        fx.write('@ NAME %s "ORBIT_X"\n')
        fx.write('@ TYPE %s "ORBIT"\n')
        fx.write('@ DELTAP %le 0\n')
        fx.write('* NAME X Y \n')
        fx.write('$ %s %le %le \n')
        
        fy.write('@ NAME %s "ORBIT_Y"\n')
        fy.write('@ TYPE %s "ORBIT"\n')
        fy.write('* NAME X Y \n')
        fy.write('$ %s %le %le  \n')
    

        reading = False

        for line in lines:
            words = line.split()

            if len(words) < 2:
                continue

            if words[0] == '#' and words[1] == 'MONITOR':
                reading = True
                continue

            if reading:
                if words[0] == '*':
                    continue

                if words[0] == '#' and words[1] != 'MONITOR':
                    reading = False
                    continue

                # otherwise monitor readings

                if words[1] == 'H' and words[8] == 'OK':
                    fx.write(words[0] + ' ' +  str(scale*eval(words[3])*1.e-6 ) + ' 0 '+ '\n')

                if words[1] == 'V' and words[8] == 'OK':
                    fy.write(words[0] + ' 0 ' + str(scale*eval(words[3])*1.e-6) + '\n')
            
        

        fx.close()
        fy.close()

    except Exception, e:
        print e

def yasp2tfsP(file, outf_x, outf_y, bpmModel):

    try:
        scale = 1
    
        lines = open(file).read().split('\n')

        fx = open(outf_x,'w')
        fy = open(outf_y,'w')
        
        fx.write('@ NAME %s "ORBIT_X"\n')
        fx.write('@ TYPE %s "ORBIT"\n')
        fx.write('@ DELTAP %le 0\n')
        fx.write('* NAME X Y \n')
        fx.write('$ %s %le %le \n')
        
        fy.write('@ NAME %s "ORBIT_Y"\n')
        fy.write('@ TYPE %s "ORBIT"\n')
        fy.write('* NAME X Y \n')
        fy.write('$ %s %le %le  \n')
    

        reading = False

        for line in lines:
            words = line.split()

            if len(words) < 2:
                continue

            if words[0] == '#' and words[1] == 'MONITOR':
                reading = True
                continue

            if reading:
                if words[0] == '*':
                    continue

                if words[0] == '#' and words[1] != 'MONITOR':
                    reading = False
                    continue

                # otherwise monitor readings

                if words[1] == 'H' and words[8] == 'OK':
                    name = words[0]
                    r = bpmModel.orbit(name.lower(), eval(words[3]),0.0)
                    if len(r) > 0:
                        fx.write(name + ' ' +  str(scale*r[0][0]*1.e-6 ) + ' 0 '+ '\n')

                if words[1] == 'V' and words[8] == 'OK':
                    name = words[0]
                    r = bpmModel.orbit(name.lower(), 0.0,eval(words[3]))
                    if len(r) > 0:
                        fy.write(name + ' 0 ' +  str(scale*r[1][0]*1.e-6 ) + '\n')
            
        

        fx.close()
        fy.close()

    except Exception, e:
        print 'exception'
        print e




def mad2dev(mad_name):

    #name_map = [['MDH.11605' ,'logical.MDHD11832'],['MDH.61804' ,'logical.MDHA6803']]

    mad_name=mad_name.replace('"','').replace('.','')

    #for name in name_map:
    #    if mad_name == name[0]:
    #        return name[1]

    name = 'logical.' + mad_name
    
    return name



# a subset of b
def is_subset(a,b):
    for e in a:
        if e not in b:
            return False
    return True

# res = a/b
def exclude(a,b):
    res = []
    for e in a:
        if e not in b:
            res.append(e)
    return res



#########################
# output drivers
#########################



        
# test
if __name__ == "__main__":
    import sys
    #r = SDDSReader(sys.argv[1])
    r = TFSReader(sys.argv[1])

    print r.get_parameter_names()
    print r.get_parameter('dx')

    w = TFSWriter()

    w.set_parameters(['x','y','z'],['%le','%le','%le'])
    w.add_record([1,2,3])
    w.add_record([4,5,6])

    w.write('test.tfs')
