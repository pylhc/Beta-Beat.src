"""This module provides all the helpers for correct_non_linearities.py.

All correction coefficients and polynomials are taken from the note:
"Geometrical non-linearity correction procedure of LHC beam position monitors", A. Nosych (BE/BI) (CERN EDMS Id:1342295)

General note:
As the correction polynomial for x plane data is _old_P5(x, coefficients) and for y plane data _old_P5(y, coefficients)
this implementation will talk of u and v as generalized for x and y. So u can be either x or y, same for v.
"""
import numpy
import scipy.optimize

# Structure: YEAR:{(BPM_TYPE1, BPM_TYPE2, ...):(kf, A, B, C), ...}
old_coefficients = {
    2011:{
        ('BPMSW', 'BPMS',  'BPMCW'         ):(16.95, 6e-6   , 3e-4  , 1.00183),
        ('BPMSX', 'BPMSA', 'BPMSY'         ):(21.9 , 2e-6   , 2e-4  , 1.037  ),
        ('BPMD',  'BPMSE', 'BPMSB', 'BPMSI'):(35.55, 2.7e-7 , 8.8e-5, 1.014  ),
        ('BPMC',  'BPMCA'                  ):(12.25, 2e-6   , 2e-4  , 1.037  ),
        ('BPMWI', 'BPMWT', 'BPMWJ'         ):(21.7 , 0      , 0     , 1      ), # note from note: not sure
        ('BPMW',  'BPMWA', 'BPMWB', 'BPMWC',
         'BPMWE', 'BPMWF', 'BPMYA', 'BPMYB',
         'BPMWG', 'BPMCW'                  ):(16   , 8e-6   , 2e-4  , 1.022  ),
        ('BPM',   'BPM_A', 'BPMR',  'BPMRA',
         'BPMC', 'BPMCA'                   ):(12.58, 2.3e-5 , 3.7e-5, 1.031  ),
        ('BPMBH'                           ):(34   , 1.08e-5, 8e-6  , 1.033  ), # note from note: overscaled
        ('BPMBV'                           ):(20.75, 1.08e-5, 8e-6  , 1.033  ), # note from note: good within +/- 5mm
        ('BPMIA'                           ):(20   , 1.08e-5, 8e-6  , 1.033  ), # note from note: good within +/- 5mm
        ('BPMI'                            ):(15.34, 1.08e-5, 8e-6  , 1.033  )  # note from note: good within +/- 5mm
    },
    2012:{
        ('BPMSW', 'BPMS',  'BPMCW'         ):(15.25 , 5.36e-6, 8.05e-4, 0.98939), # note from note: underscaled by 10%
        ('BPMSX', 'BPMSA', 'BPMSY'         ):(20.25 , 1.77e-6, 4.38e-4, 0.99751), # note from note: underscaled by 10%
        ('BPMD',  'BPMSE', 'BPMSB', 'BPMSI'):(32.75 , 2.66e-7, 1.61e-5, 1.0057 ), # note from note: underscaled by 10%
        ('BPMC',  'BPMCA'                  ):(12.25 , 2e-6   , 2e-4   , 1.037  ), # note from note: not sure
        ('BPMWI', 'BPMWT', 'BPMWJ'         ):(20.704, 6e-6   , 1e-3   , 1.0915 ),
        ('BPMW',  'BPMWA', 'BPMWB', 'BPMWC',
         'BPMWE', 'BPMWF', 'BPMYA', 'BPMYB',
         'BPMWG', 'BPMCW'                  ):(15.25 , 4.65e-6, 1.13e-3, 1.0516 ),
        ('BPM',   'BPM_A', 'BPMR',  'BPMRA',
         'BPMC', 'BPMCA'                   ):(12.68 , 2.32e-5, 1.4e-5 , 1.031  ),
        ('BPMBH'                           ):(34    , 1.08e-5, 8e-6   , 1.033  ),
        ('BPMBV'                           ):(20.75 , 1.08e-5, 8e-6   , 1.033  ),
        ('BPMIA'                           ):(20    , 1.08e-5, 8e-6   , 1.033  ),
        ('BPMI'                            ):(15.34 , 1.08e-5, 8e-6   , 1.033  )
    },
    2013:{
        ('BPMSW', 'BPMS',  'BPMCW'         ):(1     , 4.2723 , 3.7608 , 16.985),
        ('BPMSX', 'BPMSA', 'BPMSY'         ):(1     , 5.3349 , 4.9154 , 22.031),
        ('BPMD',  'BPMSE', 'BPMSB', 'BPMSI'):(1     , 8.1936 , 7.7024 , 34.657),
        ('BPMC',  'BPMCA'                  ):(1     , 3.5781 , 3.1057 , 13.934),
        ('BPMWI', 'BPMWT', 'BPMWJ'         ):(20.704, 6e-6   , 1e-3   , 1.0915),
        ('BPMW',  'BPMWA', 'BPMWB', 'BPMWC',
         'BPMWE', 'BPMWF', 'BPMYA', 'BPMYB',
         'BPMWG', 'BPMCW'                  ):(15.25 , 4.65e-6, 1.13e-3, 1.0516),
        ('BPM',   'BPM_A', 'BPMR',  'BPMRA',
         'BPMC', 'BPMCA'                   ):(12.68 , 2.32e-5, 1.4e-5 , 1.031 ),
        ('BPMBH'                           ):(34    , 1.08e-5, 8e-6   , 1.033 ),
        ('BPMBV'                           ):(20.75 , 1.08e-5, 8e-6   , 1.033 ),
        ('BPMIA'                           ):(20    , 1.08e-5, 8e-6   , 1.033 ),
        ('BPMI'                            ):(15.34 , 1.08e-5, 8e-6   , 1.033 )
    }
}

# Structure: {(BPM_TYPE1, BPM_TYPE2, ...):(A, B, C, D, E, F), ...}
new_coefficients = {
    ('BPMSW', 'BPMS',  'BPMCW'         ):( 5.1631, 3.0715, 17.1037,  8.8287,  5.3447, 1.7684),
    ('BPMSX', 'BPMSA', 'BPMSY'         ):( 6.7291, 3.8788, 22.2015, 13.2101,  7.5656, 2.4416),
    ('BPMD',  'BPMSE', 'BPMSB', 'BPMSI'):(10.8535, 5.8018, 34.9546, 24.1361, 12.9804, 4.1678),
    ('BPMC',  'BPMCA'                  ):( 4.0763, 2.6828, 14.0142,  5.5974,  3.9145, 1.2896),
    ('BPMWI', 'BPMWT', 'BPMWJ'         ):( 5.3368, 4.7708, 20.7526,  2.2937,  5.1702, 1.8677),
    ('BPMW',  'BPMWA', 'BPMWB', 'BPMWC',
     'BPMWE', 'BPMWF', 'BPMYA', 'BPMYB',
     'BPMWG', 'BPMCW'                  ):( 3.8567, 3.9523, 16.0525, -0.1955,  3.0572, 1.3772),
    ('BPM',   'BPM_A', 'BPMR',  'BPMRA',
     'BPMC', 'BPMCA'                   ):( 3.0268, 3.4651, 13.0207, -1.8265,  2.0741, 0.5032),
    ('BPMBH'                           ):(11.9862, 6.8456, 39.5071, 24.1082, 14.3634, 4.8470),
    ('BPMBV'                           ):( 5.9129, 4.2735, 21.2507,  7.5930,  6.4179, 2.6188),
    ('BPMIA'                           ):( 5.6163, 4.2104, 20.5073,  6.8736,  6.1019, 2.5247),
    ('BPMI'                            ):( 4.0045, 3.5540, 15.6567,  2.2757,  3.7624, 1.7832)
}

# Structure: {(BPM_TYPE1, BPM_TYPE2, ...):(A, B, C, D, E, F), ...}
# This 1D fallback is used according to the note when only data of one plane is available.
new_coefficients_1D_fallback = {
    ('BPMSW', 'BPMS',  'BPMCW'         ):(4.2723, 3.7608, 16.9850, 0, 0, 0),
    ('BPMSX', 'BPMSA', 'BPMSY'         ):(5.3349, 4.9154, 22.0311, 0, 0, 0),
    ('BPMD',  'BPMSE', 'BPMSB', 'BPMSI'):(8.1936, 7.7024, 34.6566, 0, 0, 0),
    ('BPMC',  'BPMCA'                  ):(3.5781, 3.1057, 13.9338, 0, 0, 0),
    ('BPMWI', 'BPMWT', 'BPMWJ'         ):(5.4824, 4.8054, 20.7159, 0, 0, 0),
    ('BPMW',  'BPMWA', 'BPMWB', 'BPMWC',
     'BPMWE', 'BPMWF', 'BPMYA', 'BPMYB',
     'BPMWG', 'BPMCW'                  ):(4.1355, 3.8386, 16.0507, 0, 0, 0),
    ('BPM',   'BPM_A', 'BPMR',  'BPMRA',
     'BPMC', 'BPMCA'                   ):(3.5215, 3.2087, 13.0402, 0, 0, 0),
    ('BPMBH'                           ):(9.4602, 8.6719, 39.2139, 0, 0, 0),
    ('BPMBV'                           ):(5.2915, 4.7688, 21.1602, 0, 0, 0),
    ('BPMIA'                           ):(5.0663, 4.6539, 20.4250, 0, 0, 0),
    ('BPMI'                            ):(3.9801, 3.6268, 15.6328, 0, 0, 0)
}

def get_coefficients_of_bpm_type(bpm_type, coefficient_table):
    """Returns the coefficient tuple for the given bpm_type from the given coefficient_table.
    If bpm_type is not in coefficient_table will return None.
    """
    for bpm_group in coefficient_table:
        for coefficient_bpm_type in bpm_group:
            if bpm_type == coefficient_bpm_type:
                return coefficient_table[bpm_group]
    return None

def _old_P5(u_raw, coefficients):
    """Old correction polynomial for 2011, 2012 and 2013 data."""
    A, B, C = coefficients[1:]
    return A * (u_raw**5) + B * (u_raw**3) + C * (u_raw)

def revert_correction_polynomial(u_bpm, coefficients, polynomial=_old_P5):
    """Reverts given polynomial (default _old_P5()) for the given u_bpm point with the given coefficients."""
    u0 = 1 # some initial value for the numeric calculation
    return scipy.optimize.newton(lambda u_raw: polynomial(u_raw, coefficients) - u_bpm, u0, tol=1e-6, maxiter=1000)

def new_P5(u_raw, v_raw, coefficients):
    """New correction polynomial for 2015 data."""
    A, B, C, D, E, F = coefficients
    return A * (u_raw**5) + B * (u_raw**3) + C * (u_raw) + D * (u_raw**3 * v_raw**2) + E * (u_raw * v_raw**4) + F * (u_raw * v_raw**2)

def new_P5_helper(uv_raw, coefficients):
    """Wrapper for new_P5() to be able to call fmin on new_P5()."""
    return new_P5(uv_raw[0], uv_raw[1], coefficients)

def revert_correction_polynomial_2D(u_bpm, v_bpm, coefficients, polynomial=new_P5_helper):
    """Reverts the new correction polynomial, probably not needed. Implemented for testing."""
    return scipy.optimize.fmin_slsqp(lambda uv_raw: abs(polynomial(uv_raw, coefficients) - u_bpm), numpy.array([0, 0]), bounds=[(-50, 50), (-50, 50)], iprint=False)[0] # works

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x_inp = numpy.linspace(-.79, .79, 21) # raw in [-1,1]
    y_inp = numpy.linspace(-.79, .79, 21)

    x_raw = [] # raw position in [-1,1]
    y_raw = []
    x_bpm = [] # measured beam position in mm
    y_bpm = []
    x_old = [] # corrected beam position in mm (2012 old correction)
    y_old = []
    x_new = [] # corrected beam position in mm (2015 new correction)
    y_new = []
    x_undo_new = [] # new correction reverted
    y_undo_new = []

    old_coefficients = get_coefficients_of_bpm_type('BPM', old_coefficients[2012])
    new_coefficients = get_coefficients_of_bpm_type('BPM', new_coefficients)
    for x in x_inp:
        for y in y_inp:
            x_raw.append(x)
            y_raw.append(y)
            x_bpm.append(x*old_coefficients[0])
            y_bpm.append(y*old_coefficients[0])
            x_old.append(_old_P5(x*old_coefficients[0], old_coefficients))
            y_old.append(_old_P5(y*old_coefficients[0], old_coefficients))
            x_new.append(new_P5(x, y, new_coefficients))
            y_new.append(new_P5(y, x, new_coefficients))
            #x_undo_new.append(revert_correction_polynomial_2D(x_new[-1], y_new[-1], new_coefficients))
            #y_undo_new.append(revert_correction_polynomial_2D(y_new[-1], x_new[-1], new_coefficients))

    #corr_dif = (numpy.sqrt(numpy.array(x_new)**2+numpy.array(y_new))-numpy.sqrt(numpy.array(x_bpm)**2+numpy.array(y_bpm)))/numpy.sqrt(numpy.array(x_bpm)**2+numpy.array(y_bpm)**2)*100

    fig = plt.figure()
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111, aspect='equal')
    #ax = fig.add_subplot(111)

    ax.scatter(numpy.array(x_old), numpy.array(y_old), c='b', alpha=.5, marker='o', s=40, linewidth=0, label='Old correction')
    ax.scatter(numpy.array(x_new), numpy.array(y_new), c='g', alpha=.5, marker='o', s=40, linewidth=0, label='New correction')
    ax.scatter(numpy.array(x_bpm), numpy.array(y_bpm), c='r', alpha=.5, marker='o', s=40, linewidth=0, label='Measured position')
    #ax.scatter(numpy.array(x_undo_new)*old_coefficients[0], numpy.array(y_undo_new)*old_coefficients[0], c='b', alpha=.5, marker='+', label='Undo new correction')
    #ax.scatter(numpy.sqrt(numpy.array(x_bpm)**2+numpy.array(y_bpm)**2), numpy.array(corr_dif), c='r', alpha=.5, marker='+', label='Deviation')

    ax.set_xlabel('$x$ [mm]')
    ax.set_ylabel('$y$ [mm]')

    ax.legend(shadow=True, fancybox=True)

    plt.subplots_adjust(left=0.09, right=0.97, top=1, bottom=0.02)
    ax.grid(which='both')

    plt.gcf().set_size_inches(7, 7)
    plt.savefig('/afs/cern.ch/user/r/rwestenb/Desktop/old_vs_new_correction.png', dpi=400)
