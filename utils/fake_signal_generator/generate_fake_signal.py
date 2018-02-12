# for creating and writing the fake BPM signals
from file_writer import write_bpm_file
import numpy as np
import random


def main():
    outFile = './test.sdds'
    bpm_map_path = './bpmMap_fake.dat'

    qx, qy = 0.28, 0.31
    turns = 4096

    bpms = {}
    for i in range(26):  # 100 BPMs
        bpms['BPMX.%d.BTEST' % i] = {'pos': i, 'plane': 'x', 'phaseOffset': 0}
        bpms['BPMY.%d.BTEST' % i] = {'pos': i, 'plane': 'y', 'phaseOffset': 0}

    # Resonances (random examples, may be adjusted)
    tunes = {
        'x': [qx],  # , 2*qy, -3*qx, -qx-qy, 2*qx-2*qy, -2*qy, qx-2*qy, -qx+3*qy, qx+2*qy, -2*qx+qy, qx+qy, 2*qx, -qx-2*qy, 3*qx],
        'y': [qy],  # , -qx-qy, -2*qx, qx-qy, -2*qy, -3*qy, -qx+qy, 2*qx+qy, -qx+3*qy, qx+qy, -qx+2*qy]
    }
#     amps   = {
#         'x':np.concatenate([[1],np.linspace(0.10,0.05,len(tunes['x'])-1)]),
#         'y':np.concatenate([[1],np.linspace(0.10,0.05,len(tunes['y'])-1)])
#     }
#     phases = {
#         'x':np.random.random(len(tunes['x']))*2*np.pi,
#         'y':np.random.random(len(tunes['y']))*2*np.pi
#     }
    amps = {
        'x': np.concatenate([[1], 0.1 * np.ones(len(tunes['x']) - 1)]),
        'y': np.concatenate([[1], 0.1 * np.ones(len(tunes['y']) - 1)])
    }
    phases = {
        'x': np.zeros(len(tunes['x'])),
        'y': np.zeros(len(tunes['y']))
    }

    write_fake_bpm_map(bpms, bpm_map_path, tunes, amps, phases)
    bpm_signal_dict = signal_generator(bpms, amps, tunes, phases, turns=turns, noise_amp=0.0)

    plot = False
    if plot:
        make_plot(bpm_signal_dict, bpms, fft=False)

    # write data to file
    wroteToFile = write_bpm_file(outFile, bpm_signal_dict, bpms)
    if wroteToFile:
        print('Wrote data file to: %s' % outFile)
    else:
        raise IOError('Something went wrong while writing the data to file!')
    print('Done!')


def noise(amp):
    return (random.random() * 2. - 1.) * amp


def signal(turn, amps, tunes, phases, noise_amp):
    assert len(amps) == len(tunes) and len(amps) == len(phases)
    pi2 = np.pi * 2

    signal = 0
    for i in range(len(tunes)):
        signal += amps[i] * np.cos(pi2 * tunes[i] * turn + phases[i])
    signal += noise(noise_amp / 2.)
    return signal


def signal_generator(bpms, amps, tunes, phases, turns=2000, noise_amp=.03):
    bpmDict = {}
    # initialize dict
    for bpm in bpms:
        bpmDict[bpm] = []
    for bpm in bpms:
        plane = bpms[bpm]['plane']
        for turn in range(turns):
            ph = [phases[plane][i] + bpms[bpm]['phaseOffset'] for i in range(len(phases[plane]))]
            bpmDict[bpm].append(signal(turn, amps[plane], tunes[plane], ph, noise_amp))
    return bpmDict


def make_plot(bpmSignalDict, bpms, fft=True):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    import numpy as np

    ax = plt.subplots()[1]
    if fft:
        x = np.fft.fftfreq(2000, d=1.)
        y = []
        for bpm in bpms:
            y.append(abs(np.fft.fft(bpmSignalDict[bpm])))
            plt.plot(x, y[-1], alpha=.7, label=bpm)
        plt.xlim(-.5, .5)
        # plt.ylim(1.)
        plt.xticks()
        plt.semilogy()
        majorLocator = MultipleLocator(0.05)
        minorLocator = MultipleLocator(0.01)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_minor_locator(minorLocator)
        plt.xlabel('Tune [1/Turns]')
        plt.ylabel('Amplitude [a.u.]')
    else:
        for bpm in bpms:
            plt.plot(bpmSignalDict[bpm], alpha=.7, label=bpm)
        plt.xlabel('Turns [1]')
        plt.ylabel('Amplitude [a.u.]')

    plt.subplots_adjust(left=0.04, right=0.97, top=0.97, bottom=0.04)
    plt.grid(which='both')
    plt.legend()
    plt.show()


def print_columns(outFile):
    import numpy as np
    from outFile.read import read_out_file
    data = read_out_file(outFile)
    for column in data:
        if 'AMP' in column:
            if np.mean(data[column]) != 0. and np.mean(data[column.replace('PHASE', 'AMP')]) > 1e-8:
                print(' ' + column + '\t: '),
                print([d * 360. for d in data[column]]),
                print(np.mean(data[column.replace('PHASE', 'AMP')]))


def write_fake_bpm_map(bpms, bpm_map_path, tunes, amps, phases):

    with open(bpm_map_path + 'x', 'w') as bpm_mapx:
        bpm_mapx.write('* NAME S AMPX AMP01 AMP_20 AMP02 AMP_30 AMP_1_1 AMP2_2 AMP0_2 AMP1_2 AMP_13 AMP12 AMP_21 AMP11 AMP20 AMP_1_2 AMP30 MUX PHASE01 PHASE_20 PHASE02 PHASE_30 PHASE_1_1 PHASE2_2 PHASE0_2 PHASE1_2 PHASE_13 PHASE12 PHASE_21 PHASE11 PHASE20 PHASE_1_2 PHASE30 \n')
        bpm_mapx.write('$ %s %le  %le %le   %le    %le   %le    %le     %le    %le    %le    %le    %le   %le    %le   %le   %le     %le   %le %le     %le      %le     %le      %le       %le      %le      %le      %le      %le     %le      %le     %le     %le       %le     \n')
        for bpm in bpms:
            phx = [phases['x'][i] + bpms[bpm]['phaseOffset'] for i in range(len(phases['x']))]
            if bpms[bpm]['plane'] == 'x':
                bpm_mapx.write(str(bpm) + ' ' + str(bpms[bpm]['pos']) + ' ' + " ".join(map(str, amps['x'])) + ' ' + " ".join(map(str, phx)) + '\n')
    with open(bpm_map_path + 'y', 'w') as bpm_mapy:
        bpm_mapy.write('* NAME S AMPY AMP10  AMP_1_1 AMP_20 AMP1_1 AMP0_2 AMP0_3 AMP_11 AMP21 AMP_13 AMP11 AMP_12 MUY PHASE10 PHASE_1_1 PHASE_20 PHASE1_1 PHASE0_2 PHASE0_3 PHASE_11 PHASE21 PHASE_13 PHASE11 PHASE_12 \n')
        bpm_mapy.write('$ %s %le %le  %le    %le     %le    %le    %le    %le    %le    %le   %le    %le   %le    %le %le     %le       %le      %le      %le      %le      %le      %le     %le      %le     %le      \n')
        for bpm in bpms:
            phy = [phases['y'][i] + bpms[bpm]['phaseOffset'] for i in range(len(phases['y']))]
            if bpms[bpm]['plane'] == 'y':
                bpm_mapy.write(str(bpm) + ' ' + str(bpms[bpm]['pos']) + ' ' + " ".join(map(str, amps['y'])) + ' ' + " ".join(map(str, phy)) + '\n')


if __name__ == '__main__':
    main()
