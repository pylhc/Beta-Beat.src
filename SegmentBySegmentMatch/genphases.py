import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Python_Classes4MAD")
from metaclass import *
import numpy


def main():

    IPN = 1
    if len(sys.argv) > 1:
        IPN = int(sys.argv[1])

    if IPN > 8 or IPN < 1:
        IPN = 1

    IPNO = str(IPN)
    print IPN

    data_beam1 = twiss('Beam1/sbs/sbsphasext_IP' + IPNO + '.out')

    data_beam2 = twiss('Beam2/sbs/sbsphasext_IP' + IPNO + '.out')

    beam1_s_list = data_beam1.MODEL_S
    sorted_index_beam1 = numpy.argsort(beam1_s_list)
    print sorted_index_beam1

    beam2_s_list = data_beam2.MODEL_S
    sorted_index_beam2 = numpy.argsort(beam2_s_list)
    print sorted_index_beam2

    print 'Start H B1: ', data_beam1.NAME[sorted_index_beam1[0]]
    print 'Start H B2: ', data_beam2.NAME[sorted_index_beam2[0]]

    phases_file = open("phases.seqx", 'w')
    phases0_beam1_file = open("phases0b1.seqx", 'w')
    phases0_beam2_file = open("phases0b2.seqx", 'w')

    print data_beam1.NAME[0], data_beam1.NAME[len(data_beam1.NAME) - 1]
    print data_beam2.NAME[0], data_beam2.NAME[len(data_beam2.NAME) - 1]

    phases_file.write('\n!!!! BEAM 1 H !!!!!\n\n')
    phases0_beam1_file.write('\n!!!! BEAM 1 H !!!!!\n\n')

    for name in data_beam1.NAME:
        phases0_beam1_file.write('mux0' + name + ' = ')
        phases0_beam1_file.write('table(twiss, ' + name + ', mux) - ')
        phases0_beam1_file.write('table(twiss, ' + data_beam1.NAME[sorted_index_beam1[0]] + ', mux);\n')

        phases_file.write('mux' + name + ' := ')
        phases_file.write('table(twiss, ' + name + ', mux) - ')
        phases_file.write('table(twiss, ' + data_beam1.NAME[sorted_index_beam1[0]] + ', mux);\n')

    phases_file.write('\n')

    for name in data_beam1.NAME:
        phases_file.write('dmux' + name + ' := ')
        phases_file.write('mux' + name + ' - ' 'mux0' + name + ';\n')

    phases_file.write('\n!!!! BEAM 2 H !!!!!\n\n')

    phases0_beam2_file.write('\n!!!! BEAM 2 H !!!!!\n\n')

    for idx in sorted_index_beam2:
        name = data_beam2.NAME[idx]

        phases0_beam2_file.write('mux0' + name + ' = ')
        phases0_beam2_file.write('table(twiss, ' + name + ', mux) - ')
        phases0_beam2_file.write('table(twiss, ' + data_beam2.NAME[sorted_index_beam2[0]] + ', mux);\n')

        phases_file.write('mux' + name + ' := ')
        phases_file.write('table(twiss, ' + name + ', mux) - ')
        phases_file.write('table(twiss, ' + data_beam2.NAME[sorted_index_beam2[0]] + ', mux);\n')

    phases_file.write('\n')

    for name in data_beam2.NAME:
        phases_file.write('dmux' + name + ' := ')
        phases_file.write('mux' + name + ' - ' 'mux0' + name + ';\n')

    data_beam1 = twiss('Beam1/sbs/sbsphaseyt_IP' + IPNO + '.out')
    data_beam2 = twiss('Beam2/sbs/sbsphaseyt_IP' + IPNO + '.out')

    beam1_s_list = data_beam1.MODEL_S
    sorted_index_beam1 = numpy.argsort(beam1_s_list)
    print sorted_index_beam1

    beam2_s_list = data_beam2.MODEL_S
    sorted_index_beam2 = numpy.argsort(beam2_s_list)
    print sorted_index_beam2

    print 'Start V B1: ', data_beam1.NAME[sorted_index_beam1[0]]
    print 'Start V B2: ', data_beam2.NAME[sorted_index_beam2[0]]

    phases_file.write('\n!!!! BEAM 1 V !!!!!\n\n')

    phases0_beam1_file.write('\n!!!! BEAM 1 V !!!!!\n\n')

    for name in data_beam1.NAME:

        phases0_beam1_file.write('muy0' + name + ' = ')
        phases0_beam1_file.write('table(twiss, ' + name + ', muy) - ')
        phases0_beam1_file.write('table(twiss, ' + data_beam1.NAME[sorted_index_beam1[0]] + ', muy);\n')

        phases_file.write('muy' + name + ' := ')
        phases_file.write('table(twiss, ' + name + ', muy) - ')
        phases_file.write('table(twiss, ' + data_beam1.NAME[sorted_index_beam1[0]] + ', muy);\n')

    phases_file.write('\n')

    for name in data_beam1.NAME:
        phases_file.write('dmuy' + name + ' := ')
        phases_file.write('muy' + name + ' - ' 'muy0' + name + ';\n')

    phases_file.write('\n!!!! BEAM 2 V !!!!!\n\n')

    phases0_beam2_file.write('\n!!!! BEAM 2 V !!!!!\n\n')

    for name in data_beam2.NAME:
        print name

        phases0_beam2_file.write('muy0' + name + ' = ')
        phases0_beam2_file.write('table(twiss, ' + name + ', muy) - ')
        phases0_beam2_file.write('table(twiss, ' + data_beam2.NAME[sorted_index_beam2[0]] + ', muy);\n')

        phases_file.write('muy' + name + ' := ')
        phases_file.write('table(twiss, ' + name + ', muy) - ')
        phases_file.write('table(twiss, ' + data_beam2.NAME[sorted_index_beam2[0]] + ', muy);\n')

    phases_file.write('\n')

    for name in data_beam2.NAME:
        phases_file.write('dmuy' + name + ' := ')
        phases_file.write('muy' + name + ' - ' 'muy0' + name + ';\n')

    phases0_beam1_file.close()
    phases0_beam2_file.close()
    phases_file.close()


if __name__ == "__main__":
    main()
