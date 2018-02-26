import numpy as np


def get_r_terms(name, meas):
    betx = meas.beta_x.BETX[name]
    bety = meas.beta_y.BETY[name]
    alfx = meas.beta_x.ALFX[name]
    alfy = meas.beta_y.ALFY[name]
    f1001r = meas.coupling.F1001R[name]
    f1001i = meas.coupling.F1001I[name]
    f1010r = meas.coupling.F1010R[name]
    f1010i = meas.coupling.F1010I[name]
    if (f1001r, f1001i, f1010r, f1010i) == (0., 0., 0., 0.):
        return 0., 0., 0., 0.,

    ga11 = 1 / np.sqrt(betx)
    ga12 = 0
    ga21 = alfx / np.sqrt(betx)
    ga22 = np.sqrt(betx)
    Ga = np.reshape(np.array([ga11, ga12, ga21, ga22]), (2, 2))

    gb11 = 1 / np.sqrt(bety)
    gb12 = 0
    gb21 = alfy / np.sqrt(bety)
    gb22 = np.sqrt(bety)
    Gb = np.reshape(np.array([gb11, gb12, gb21, gb22]), (2, 2))

    J = np.reshape(np.array([0, 1, -1, 0]), (2, 2))

    absf1001 = np.sqrt(f1001r ** 2 + f1001i ** 2)
    absf1010 = np.sqrt(f1010r ** 2 + f1010i ** 2)

    gamma2 = 1. / (1. + 4. * (absf1001 ** 2 - absf1010 ** 2))
    c11 = (f1001i + f1010i)
    c22 = (f1001i - f1010i)
    c12 = -(f1010r - f1001r)
    c21 = -(f1010r + f1001r)
    Cbar = np.reshape(2 * np.sqrt(gamma2) *
                         np.array([c11, c12, c21, c22]), (2, 2))

    C = np.dot(np.linalg.inv(Ga), np.dot(Cbar, Gb))
    jCj = np.dot(J, np.dot(C, -J))
    c = np.linalg.det(C)
    r = -c / (c - 1)
    R = np.transpose(np.sqrt(1 + r) * jCj)
    r11, r12, r21, r22 = np.ravel(R)
    return r11, r12, r21, r22


def propagate_error_beta(errb0, erra0, dphi, bets, bet0, alf0):
    return np.sqrt(
        (
            bets*np.sin(4*np.pi*dphi)*alf0/bet0 +
            bets*np.cos(4*np.pi*dphi)/bet0
        )**2*errb0**2 +
        (bets*np.sin(4*np.pi*dphi))**2*erra0**2
    )


def propagate_error_alfa(errb0, erra0, dphi, alfs, bet0, alf0):
    return np.sqrt(
        (
            (alfs*((np.sin(4*np.pi*dphi)*alf0/bet0) + (np.cos(4*np.pi*dphi)/bet0))) -
            (np.cos(4*np.pi*dphi)*alf0/bet0) + (np.sin(4*np.pi*dphi)/bet0)
        )**2*errb0**2 +
        ((np.cos(4*np.pi*dphi)) - (alfs*np.sin(4*np.pi*dphi)))**2*erra0**2
    )


def propagate_error_phase(errb0, erra0, dphi, bet0, alf0):
    return np.sqrt(
        (
            ((1/2.*np.cos(4*np.pi*dphi)*alf0/bet0)-
             (1/2.*np.sin(4*np.pi*dphi)/bet0)
             -(1/2.*alf0/bet0)
            )*errb0
        )**2+
        ((-(1/2.*np.cos(4*np.pi*dphi))+(1/2.))*erra0)**2
    )/(2*np.pi)


def propagate_error_coupling_1001_re(f1001ab_ini, p1001_ini, phasex, phasey, f1001_std_ini, p1001_std_ini):
    return np.sqrt(
        (f1001_std_ini * np.cos(2 * np.pi * (p1001_ini - phasex + phasey)))**2 +
        (2 * np.pi * p1001_std_ini * f1001ab_ini *
         np.sin(2 * np.pi * (p1001_ini - phasex + phasey)))**2
    )


def propagate_error_coupling_1001_im(f1001ab_ini, p1001_ini, phasex, phasey, f1001_std_ini, p1001_std_ini):
    return np.sqrt(
        (f1001_std_ini * np.sin(2 * np.pi * (p1001_ini - phasex + phasey)))**2 +
        (2 * np.pi * p1001_std_ini * f1001ab_ini *
         np.cos(2 * np.pi * (p1001_ini - phasex + phasey)))**2
    )



def propagate_error_coupling_1010_re(f1010ab_ini, p1010_ini, phasex, phasey, f1010_std_ini, p1010_std_ini):
    return np.sqrt(
        (f1010_std_ini * np.cos(2 * np.pi * (p1010_ini - phasex - phasey)))**2 +
        (2 * np.pi * p1010_std_ini * f1010ab_ini *
         np.sin(2 * np.pi * (p1010_ini - phasex - phasey)))**2
    )



def propagate_error_coupling_1010_im(f1010ab_ini, p1010_ini, phasex, phasey, f1010_std_ini, p1010_std_ini):
    return np.sqrt(
        (f1010_std_ini * np.sin(2 * np.pi * (p1010_ini - phasex - phasey)))**2 +
        (2 * np.pi * p1010_std_ini * f1010ab_ini *
         np.cos(2 * np.pi * (p1010_ini - phasex - phasey)))**2
    )


def propagate_error_dispersion(std_D0, bet0, bets, dphi, alf0):
    return np.abs(
        std_D0 * np.sqrt(bets/bet0) *
        (np.cos(2*np.pi*dphi)+alf0*np.sin(2*np.pi*dphi))
    )


def weighted_average_for_SbS_elements(value1, sigma1, value2, sigma2):
    weighted_average = ((1/sigma1**2 * value1 + 1/sigma2**2 * value2) /
                        (1/sigma1**2 + 1/sigma2**2))
    uncertainty_of_average = np.sqrt(1 / (1/sigma1**2 + 1/sigma2**2))
    weighted_rms = np.sqrt(2 * (1/sigma1**2 * (value1 - weighted_average)**2 +
                                1/sigma2**2 * (value2 - weighted_average)**2) /
                           (1/sigma1**2 + 1/sigma2**2))
    final_error = np.sqrt(uncertainty_of_average**2 + weighted_rms**2)
    return weighted_average, final_error
