import sys
import __init__  # @UnusedImport used for appending paths
from Python_Classes4MAD.metaclass import twiss


def main():

    IPN = 1
    if len(sys.argv) > 1:
        IPN = int(sys.argv[1])

    if IPN > 8 or IPN < 1:
        IPN = 1

    IPNO = str(IPN)
    print IPN

    fulldataxb1 = twiss('Beam1/getphasex.out')
    print 'B1: Q1=', fulldataxb1.Q1
    print 'B1: Q2=', fulldataxb1.Q2

    fulldataxb2 = twiss('Beam2/getphasex.out')
    print 'B2: Q1=', fulldataxb2.Q1
    print 'B2: Q2=', fulldataxb2.Q2

    #driven is closer to what it should be than the free
    B1QX = fulldataxb1.Q1
    B1QY = fulldataxb1.Q2
    B2QX = fulldataxb2.Q1
    B2QY = fulldataxb2.Q2

    dataxb1 = twiss('Beam1/sbs/sbsphasext_IP' + IPNO + '.out')
    datayb1 = twiss('Beam1/sbs/sbsphaseyt_IP' + IPNO + '.out')
    dataxb2 = twiss('Beam2/sbs/sbsphasext_IP' + IPNO + '.out')
    datayb2 = twiss('Beam2/sbs/sbsphaseyt_IP' + IPNO + '.out')

    print dataxb2.NAME[0], dataxb2.PHASEX[0], dataxb2.PHASEXT[0]

    foutb1 = open('constraintsb1.seqx', 'w')
    foutb2 = open('constraintsb2.seqx', 'w')

    fdump1 = open("dumpb1.seqx", 'w')
    fdump2 = open("dumpb2.seqx", 'w')

    foutb1.write('\n!!!! BEAM 1 H !!!!!\n\n')
    fdump1.write('delete,table=dmuxb1;\n')
    fdump1.write('create,table=dmuxb1,column=sss,cdmux,tdmux;\n')

    i = 0
    for name in dataxb1.NAME:
        w = 1

        ckstr = '+ 0.0'
        if ((dataxb1.S[i] > 214.0 and dataxb1.S[i] < 681.0) and IPN == 2.0):
            dataxb1.PHASEXT[i] = dataxb1.PHASEXT[i] + B1QX
            ckstr = '+ B1 Qx(' + str(B1QX) + ')'

        if (abs(dataxb1.PHASEXT[i]) > 0.25):
            w = 1e-6

        foutb1.write('   constraint, weight = ' + str(w) + ' , ')
        foutb1.write('expr =  dmux' + name + ' = ' + str(dataxb1.PHASEXT[i]) + '; ')
        foutb1.write('!   S = ' + str(dataxb1.S[i]) + ' ' + ckstr)
        foutb1.write(';\n')

        fdump1.write('   sss = table(twiss, ' + name + ', s); cdmux = dmux' + name + ';  tdmux =  ' + str(dataxb1.PHASEXT[i]) + ';\n')
        fdump1.write('   fill,table=dmuxb1;\n')

        i = i + 1

    fdump1.write('write,table=dmuxb1, file=dmuxb1.tfs;\n')

    foutb1.write('\n!!!! BEAM 1 V !!!!!\n\n')
    fdump1.write('delete,table=dmuyb1;\n')
    fdump1.write('create,table=dmuyb1,column=sss,cdmuy,tdmuy;\n')

    i = 0
    for name in datayb1.NAME:
        w = 1

        ckstr = '+ 0.0'
        if ((datayb1.S[i] > 214.0 and datayb1.S[i] < 681.0) and IPN == 2.0):
            datayb1.PHASEYT[i] = datayb1.PHASEYT[i] + B1QY
            ckstr = '+ B1 Qy(' + str(B1QY) + ')'

        if (abs(datayb1.PHASEYT[i]) > 0.25):
            w = 1e-6

        foutb1.write('   constraint, weight = ' + str(w) + ' , ')
        foutb1.write('expr =  dmuy' + name + ' = ' + str(datayb1.PHASEYT[i]) + '; ')
        foutb1.write('!   S = ' + str(datayb1.S[i]) + ' ' + ckstr)
        foutb1.write(';\n')

        fdump1.write('   sss = table(twiss, ' + name + ', s); cdmuy = dmuy' + name + ';  tdmuy =  ' + str(datayb1.PHASEYT[i]) + ';\n')
        fdump1.write('   fill,table=dmuyb1;\n')

        i = i + 1

    fdump1.write('write,table=dmuyb1, file=dmuyb1.tfs;\n')

    foutb2.write('\n!!!! BEAM 2 H !!!!!\n\n')
    fdump2.write('delete,table=dmuxb2;\n')
    fdump2.write('create,table=dmuxb2,column=sss,cdmux,tdmux;\n')

    i = 0
    for name in dataxb2.NAME:
        w = 1
        ckstr = '+ 0.0'
        if ((dataxb2.S[i] < 460.0 or dataxb2.S[i] > 26532.0) and IPN == 8.0):
            dataxb2.PHASEXT[i] = dataxb2.PHASEXT[i] + B2QX
            ckstr = '+ B2 Qx(' + str(B2QX) + ')'

        if (abs(dataxb2.PHASEXT[i]) > 0.25):
            w = 1e-6

        foutb2.write('   constraint, weight = ' + str(w) + ' , ')
        foutb2.write('expr =  dmux' + name + ' = ' + str(dataxb2.PHASEXT[i]) + ';')

        foutb2.write('!   S = ' + str(dataxb2.S[i]) + ' ' + ckstr)
        foutb2.write(';\n')

        fdump2.write('   sss = table(twiss, ' + name + ', s); cdmux = dmux' + name + ';  tdmux =  ' + str(dataxb2.PHASEXT[i]) + ';\n')
        fdump2.write('   fill,table=dmuxb2;\n')

        i = i + 1

    fdump2.write('write,table=dmuxb2, file=dmuxb2.tfs;\n')

    foutb2.write('\n!!!! BEAM 2 V !!!!!\n\n')
    fdump2.write('delete,table=dmuyb2;\n')
    fdump2.write('create,table=dmuyb2,column=sss,cdmuy,tdmuy;\n')

    i = 0
    for name in datayb2.NAME:
        w = 1
        ckstr = '+ 0.0'
        if ((datayb2.S[i] < 460.0 or datayb2.S[i] > 26532.0) and IPN == 8.0):
            datayb2.PHASEYT[i] = datayb2.PHASEYT[i] + B2QY
            ckstr = '+ B2 Qy(' + str(B2QY) + ')'

        if (abs(datayb2.PHASEYT[i]) > 0.25):
            w = 1e-6

        foutb2.write('   constraint, weight = ' + str(w) + ' , ')
        foutb2.write('expr =  dmuy' + name + ' = ' + str(datayb2.PHASEYT[i]) + ';')

        foutb2.write('!   S = ' + str(datayb2.S[i]) + ' ' + ckstr)
        foutb2.write(';\n')

        fdump2.write('   sss = table(twiss, ' + name + ', s); cdmuy = dmuy' + name + ';  tdmuy =  ' + str(datayb2.PHASEYT[i]) + ';\n')
        fdump2.write('   fill,table=dmuyb2;\n')

        i = i + 1

    foutb1.close()
    foutb2.close()

    fdump2.write('write,table=dmuyb2, file=dmuyb2.tfs;\n')
    fdump2.close()

    fdump1.close()

    print ' !!!! END GENCONSTRAINTS !!!!!\n'


if __name__ == "__main__":
    main()
