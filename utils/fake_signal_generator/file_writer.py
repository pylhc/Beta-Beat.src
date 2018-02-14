import os.path
import time
import datetime

HEADER = '''#SDDSASCIIFORMAT v1
#Beam: Test
#Created: %s By: Drive Fake Signal Generator
#bunchid :0
#number of turns :%d
#number of monitors :%d
#NTURNS calculated: %d
'''


def write_bpm_file(filePath, bpmData, bpms):
    header = HEADER % (
        datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d#%H-%M-%S'),
        len(bpmData[bpmData.keys()[0]]),
        len(bpms),
        len(bpmData[bpmData.keys()[0]])
    )
    with open(filePath, 'w') as out:
        out.write(header)
        for bpm in bpmData:
            if bpms[bpm]['plane'] == 'x':
                plane = 0
            elif bpms[bpm]['plane'] == 'y':
                plane = 1
            else:
                raise ValueError('Unknown plane: %s' % bpms[bpm]['plane'])
            out.write('%d %s %f %s\n' % (plane, bpm, bpms[bpm]['pos'], ' '.join(format(x.real, '.5f') for x in bpmData[bpm])))
    if os.path.isfile(filePath):
        return True
    return False
