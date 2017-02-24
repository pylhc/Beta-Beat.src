'''
Created April, 2014

:maintainer: Yngve Inntjore Levinsen

:author: Yngve Inntjore Levinsen (based on awk script by Rogelio Tomas)

'''

from __future__ import print_function
from numpy.random import normal, randint

def convert_files( nparticles=1, infile='trackone', outfile='ALLBPMs',
                   x_error=0.0, y_error=0.0, n_faulty=0., turnrange=None):
    '''
    Similar to the old awk script ALLBPMs.

    The main difference is that this can handle track files with more than one particle.
    In this case it will write each particle track to different files

    :param int nparticles: Number of particles in the track file
    :param string infile: Name of input file
    :param string outfile: Name of output file
    :param float x_error: Sigma of noise in horizontal plane
    :param float y_error: Sigma of noise in vertical plane
    :param int n_faulty: Number of faulty BPMs
    '''

    # Print information:
    print("Parsing Mad-X track file %s to %s"%(infile,outfile))
    if x_error or y_error:
        print("Adding random errors with sigma %f horizontal, and %f vertical"%(x_error,y_error))
    if n_faulty:
        print("Adding %i faulty BPM's"%n_faulty)

    bpm_name=''
    x=[{} for i in xrange(nparticles)]
    y=[{} for i in xrange(nparticles)]
    bpms=[]
    bpm_loc=[]


    for l in file(infile):
        if l.strip()[0] in ['@','*','$']:
            continue
        lsp=l.split()
        if lsp[0]=='#segment':
                bpm_name=lsp[-1].upper()
                if bpm_name not in x[0]:
                    for i in xrange(len(x)):
                        x[i][bpm_name]=[]
                        y[i][bpm_name]=[]
        elif 'BPM' in bpm_name or 'PICK' in bpm_name:  # PICK for ESRF
                pid = int(lsp[0]) - 1
                if bpm_name not in bpms:
                    bpms.append(bpm_name)
                    bpm_loc.append(lsp[-2])
                if x_error and y_error:
                    x[pid][bpm_name].append(str(float(lsp[2]) * 1000 + normal(0.0, x_error)))
                    y[pid][bpm_name].append(str(float(lsp[4]) * 1000 + normal(0.0, y_error)))
                else:
                    x[pid][bpm_name].append(str(float(lsp[2]) * 1000))
                    y[pid][bpm_name].append(str(float(lsp[4]) * 1000))


    if n_faulty:
        failing_bpms=randint(0,len(bpms)*2,n_faulty)
    else:
        failing_bpms=[]

    for pid in xrange(nparticles):
        print("Writing particle",pid+1)
        if nparticles>1:
            fout=file('%s_%i'%(outfile,pid+1),'w')
        else:
            fout=file(outfile,'w')
        fout.write('# title\n')
        for i in xrange(len(bpms)):
            bpm=bpms[i]

            s=bpm_loc[i]

            if i not in failing_bpms:
                if turnrange is None:
                    fout.write(' '.join(['0',bpm,s]+x[pid][bpm])+'\n')
                else:
                    fout.write('0 {0:s} {1:s} '.format(bpm, str(s)))
                    for j in range(max(0, turnrange[0]), min(turnrange[1], len(x[pid][bpm]))):
                        fout.write(str(x[pid][bpm][j]) + " ")
                    fout.write("\n")
            if i+len(bpms) not in failing_bpms:
                if turnrange is None:
                    fout.write(' '.join(['1',bpm,s]+x[pid][bpm])+'\n')
                else:
                    fout.write('1 {0:s} {1:s} '.format(bpm, str(s)))
                    for j in range(max(0, turnrange[0]), min(turnrange[1], len(y[pid][bpm]))):
                        fout.write(str(y[pid][bpm][j]) + " ")
                    fout.write("\n")


if __name__=="__main__":
    convert_files()
