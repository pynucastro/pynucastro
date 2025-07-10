import numpy as np


class NeutrinoComponents:
    """A simple container that holds the individual components to the
    neutrino cooling."""

    def __init__(self):
        self.splas = None
        self.spair = None
        self.sphot = None
        self.sbrem = None
        self.sreco = None

    def __str__(self):
        return f"splas: {self.splas}; spair: {self.spair}; sphot: {self.sphot}; sbrem: {self.sbrem}; sreco: {self.sreco}"


def ifermi12(f):
    """apply a rational function expansion to get the inverse
    fermi-dirac integral of order 1/2 when it is equal to f.
    maximum error is 4.19e-9_rt.
    reference: Antia ApJS 84,101 1993"""

    # coefficients of the expansion from Table 8 of Antia

    m1 = 3
    k1 = 2
    m2 = 5
    k2 = 4

    a1 = np.array([1.999266880833e4, 5.702479099336e3, 6.610132843877e2, 3.818838129486e1, 1.0e0])

    b1 = np.array([1.771804140488e4, -2.014785161019e3, 9.130355392717e1, -1.670718177489e0])

    a2 = np.array([-1.277060388085e-2,
                   7.187946804945e-2,
                   -4.262314235106e-1,
                   4.997559426872e-1,
                   -1.285579118012e0,
                   -3.930805454272e-1,
                   1.0e0])

    b2 = np.array([-9.745794806288e-3,
                   5.485432756838e-2,
                   -3.299466243260e-1,
                   4.077841975923e-1,
                   -1.145531476975e0,
                   -6.067091689181e-2])

    if f < 4.0:

        # build sum_{i=0, m1} a_i x**i.  This is the numerator
        # in Eq. 4 of Antia.
        #
        # with our 1-based indexing, this expression is
        # a[1] + a[2] * f + a[3] * f**2 + ... a[m1+1] * f**m1
        #
        # we do the sum starting with the largest term and working
        # on a single power of f each iteration.
        #
        # in the starting rn here, the leading f is actually
        # a1(m1+1) * f, but that element is 1

        rn = f + a1[m1]

        for i in reversed(range(m1)):
            rn = rn * f + a1[i]

        # now we do the denominator in Eq. 4.  None of the coefficients
        # are 1, so we loop over all

        den = b1[k1 + 1]

        for i in reversed(range(k1 + 1)):
            den = den * f + b1[i]

        # Eq. 6 of Antia

        ifermi12r = np.log(f * rn / den)

    else:

        # this construction is the same as above, but using the
        # second set of coefficients

        an = 0.5
        ff = 1.0 / f ** (1.0 / (1.0 + an))

        rn = ff + a2[m2]

        for i in reversed(range(m2)):
            rn = rn * ff + a2[i]

        den = b2[k2 + 1]

        for i in reversed(range(k2 + 1)):
            den = den * ff + b2[i]

        ifermi12r = rn / (den * ff)

    return ifermi12r


def sneut5(rho, T, comp=None, *, abar=None, zbar=None,
           full_output=False):
    """compute thermal neutrino losses from the analytic fits of Itoh
    et al. ApJS 102, 411, 1996.  Note that either a Composition object
    of abar/zbar need to be provided.

    """

    # input:
    # T = temperature
    # rho  = density
    # abar = mean atomic weight
    # zbar = mean charge

    # output:
    # snu    = total neutrino loss rate in erg/g/sec

    if abar is None or zbar is None:
        abar = comp.abar
        zbar = comp.zbar

    # numerical constants

    fac1 = 5.0 * np.pi / 3.0
    fac2 = 10.0 * np.pi
    fac3 = np.pi / 5.0
    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0
    con1 = 1.0 / 5.9302e0
    sixth = 1.0 / 6.0

    # theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
    # xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
    # change theta and xnufam if need be, and the changes will automatically
    # propagate through the routine. cv and ca are the vector and axial currents.

    theta = 0.2319
    xnufam = 3.0
    cv = 0.5 + 2.0 * theta
    cvp = 1.0 - cv
    ca = 0.5
    cap = 1.0 - ca
    tfac1 = cv * cv + ca * ca + (xnufam - 1.0) * (cvp * cvp + cap * cap)
    tfac2 = cv * cv - ca * ca + (xnufam - 1.0) * (cvp * cvp - cap * cap)
    tfac3 = tfac2 / tfac1
    tfac4 = 0.5 * tfac1
    tfac5 = 0.5e0 * tfac2
    tfac6 = cv * cv + 1.5e0 * ca * ca + (xnufam - 1.0) * (cvp * cvp + 1.5 * cap * cap)

    snu = 0.0

    if T < 1.0e7:
        return snu

    # to avoid lots of divisions

    deni = 1.0 / rho
    tempi = 1.0 / T
    abari = 1.0 / abar

    # some composition variables

    ye = zbar * abari

    # some frequent factors

    t9 = T * 1.0e-9
    xl = t9 * con1
    xlp5 = np.sqrt(xl)
    xl2 = xl * xl
    xl3 = xl2 * xl
    xl4 = xl3 * xl
    xl5 = xl4 * xl
    xl6 = xl5 * xl
    xl7 = xl6 * xl
    xl8 = xl7 * xl
    xl9 = xl8 * xl
    xlm1 = 1.0 / xl
    xlm2 = xlm1 * xlm1
    xlm3 = xlm1 * xlm2

    rm = rho * ye
    rmi = 1.0 / rm

    a0 = rm * 1.0e-9
    a1 = a0**oneth
    zeta = a1 * xlm1
    zeta2 = zeta * zeta
    zeta3 = zeta2 * zeta

    # pair neutrino section
    # for reactions like e+ + e- => nu_e + nubar_e

    # equation 2.8

    gl = 1.0 - 13.04 * xl2 + 133.5 * xl4 + 1534.0 * xl6 + 918.6 * xl8

    # equation 2.7

    a1 = 6.002e19 + 2.084e20 * zeta + 1.872e21 * zeta2

    if t9 < 10.0:
        b1 = np.exp(-5.5924 * zeta)
    else:
        b1 = np.exp(-4.9924 * zeta)

    xnum = a1 * b1

    if t9 < 10.0:
        a1 = 9.383e-1 * xlm1 - 4.141e-1 * xlm2 + 5.829e-2 * xlm3
    else:
        a1 = 1.2383 * xlm1 - 8.141e-1 * xlm2

    xden = zeta3 + a1

    a1 = 1.0 / xden
    fpair = xnum * a1

    # equation 2.6

    a1 = 10.7480 * xl2 + 0.3967 * xlp5 + 1.005
    xnum = 1.0 / a1

    a1 = 7.692e7 * xl3 + 9.715e6 * xlp5

    c = 1.0 / a1
    b1 = 1.0 + rm * c

    xden = b1**-0.3

    qpair = xnum * xden

    # equation 2.5

    a1 = np.exp(-2.0 * xlm1)

    spair = a1 * fpair

    a1 = spair
    spair = gl * a1

    a1 = tfac4 * (1.0 + tfac3 * qpair)

    a3 = spair
    spair = a1 * a3

    # plasma neutrino section
    # for collective reactions like gamma_plasmon => nu_e + nubar_e
    # equation 4.6

    a1 = 1.019e-6 * rm
    a2 = a1**twoth
    b1 = np.sqrt(1.0 + a2)

    c00 = 1.0 / (T * T * b1)

    gl2 = 1.1095e11 * rm * c00

    gl = np.sqrt(gl2)
    gl12 = np.sqrt(gl)
    gl32 = gl * gl12
    gl72 = gl2 * gl32
    gl6 = gl2 * gl2 * gl2

    # equation 4.7

    ft = 2.4 + 0.6 * gl12 + 0.51 * gl + 1.25 * gl32

    # equation 4.8

    a1 = 8.6 * gl2 + 1.35 * gl72

    b1 = 225.0 - 17.0 * gl + gl2

    c = 1.0 / b1
    fl = a1 * c

    # equation 4.9 and 4.10

    cc = np.log10(2.0 * rm)
    xlnt = np.log10(T)

    xnum = sixth * (17.5 + cc - 3.0 * xlnt)
    xden = sixth * (-24.5 + cc + 3.0 * xlnt)

    # equation 4.11

    if np.abs(xnum) > 0.7 or xden < 0.0:
        fxy = 1.0
    else:
        a1 = 0.39 - 1.25 * xnum - 0.35 * np.sin(4.5 * xnum)

        b1 = 0.3 * np.exp(-(4.5 * xnum + 0.9)**2)

        c = min(0.0, xden - 1.6 + 1.25 * xnum)

        d = 0.57 - 0.25 * xnum
        a3 = c / d
        c00 = np.exp(-a3 * a3)

        fxy = 1.05 + (a1 - b1) * c00

    # equation 4.1 and 4.5

    splas = (ft + fl) * fxy

    a2 = np.exp(-gl)

    a1 = splas
    splas = a2 * a1

    a2 = gl6

    a1 = splas
    splas = a2 * a1

    a2 = 0.93153 * 3.0e21 * xl9

    a1 = splas
    splas = a2 * a1

    # photoneutrino process section for reactions like e- + gamma =>
    # e- + nu_e + nubar_e e+ + gamma => e+ + nu_e + nubar_e equation
    # 3.8 for tau, equation 3.6 for cc, and table 2 written out for
    # speed

    if T < 1.0e8:

        # note: we already bailed above for T < 1.e7, so this is
        # really 1.e7 <= T < 1.e8

        tau = np.log10(T * 1.0e-7)
        cc = 0.5654 + tau
        c00 = 1.008e11
        c01 = 0.0
        c02 = 0.0
        c03 = 0.0
        c04 = 0.0
        c05 = 0.0
        c06 = 0.0
        c10 = 8.156e10
        c11 = 9.728e8
        c12 = -3.806e9
        c13 = -4.384e9
        c14 = -5.774e9
        c15 = -5.249e9
        c16 = -5.153e9
        c20 = 1.067e11
        c21 = -9.782e9
        c22 = -7.193e9
        c23 = -6.936e9
        c24 = -6.893e9
        c25 = -7.041e9
        c26 = -7.193e9
        dd01 = 0.0
        dd02 = 0.0
        dd03 = 0.0
        dd04 = 0.0
        dd05 = 0.0
        dd11 = -1.879e10
        dd12 = -9.667e9
        dd13 = -5.602e9
        dd14 = -3.370e9
        dd15 = -1.825e9
        dd21 = -2.919e10
        dd22 = -1.185e10
        dd23 = -7.270e9
        dd24 = -4.222e9
        dd25 = -1.560e9

    elif 1.0e8 <= T < 1.0e9:

        tau = np.log10(T * 1.0e-8)
        cc = 1.5654
        c00 = 9.889e10
        c01 = -4.524e8
        c02 = -6.088e6
        c03 = 4.269e7
        c04 = 5.172e7
        c05 = 4.910e7
        c06 = 4.388e7
        c10 = 1.813e11
        c11 = -7.556e9
        c12 = -3.304e9
        c13 = -1.031e9
        c14 = -1.764e9
        c15 = -1.851e9
        c16 = -1.928e9
        c20 = 9.750e10
        c21 = 3.484e10
        c22 = 5.199e9
        c23 = -1.695e9
        c24 = -2.865e9
        c25 = -3.395e9
        c26 = -3.418e9
        dd01 = -1.135e8
        dd02 = 1.256e8
        dd03 = 5.149e7
        dd04 = 3.436e7
        dd05 = 1.005e7
        dd11 = 1.652e9
        dd12 = -3.119e9
        dd13 = -1.839e9
        dd14 = -1.458e9
        dd15 = -8.956e8
        dd21 = -1.548e10
        dd22 = -9.338e9
        dd23 = -5.899e9
        dd24 = -3.035e9
        dd25 = -1.598e9

    else:

        # T > 1.e9

        tau = np.log10(t9)
        cc = 1.5654
        c00 = 9.581e10
        c01 = 4.107e8
        c02 = 2.305e8
        c03 = 2.236e8
        c04 = 1.580e8
        c05 = 2.165e8
        c06 = 1.721e8
        c10 = 1.459e12
        c11 = 1.314e11
        c12 = -1.169e11
        c13 = -1.765e11
        c14 = -1.867e11
        c15 = -1.983e11
        c16 = -1.896e11
        c20 = 2.424e11
        c21 = -3.669e9
        c22 = -8.691e9
        c23 = -7.967e9
        c24 = -7.932e9
        c25 = -7.987e9
        c26 = -8.333e9
        dd01 = 4.724e8
        dd02 = 2.976e8
        dd03 = 2.242e8
        dd04 = 7.937e7
        dd05 = 4.859e7
        dd11 = -7.094e11
        dd12 = -3.697e11
        dd13 = -2.189e11
        dd14 = -1.273e11
        dd15 = -5.705e10
        dd21 = -2.254e10
        dd22 = -1.551e10
        dd23 = -7.793e9
        dd24 = -4.489e9
        dd25 = -2.185e9

    # equation 3.7, compute the expensive trig functions only one time

    cos1 = np.cos(fac1 * tau)
    cos2 = np.cos(fac1 * 2.0 * tau)
    cos3 = np.cos(fac1 * 3.0 * tau)
    cos4 = np.cos(fac1 * 4.0 * tau)
    cos5 = np.cos(fac1 * 5.0 * tau)
    last = np.cos(fac2 * tau)

    sin1 = np.sin(fac1 * tau)
    sin2 = np.sin(fac1 * 2.0 * tau)
    sin3 = np.sin(fac1 * 3.0 * tau)
    sin4 = np.sin(fac1 * 4.0 * tau)
    sin5 = np.sin(fac1 * 5.0 * tau)

    a0 = (0.5 * c00 +
          c01 * cos1 + dd01 * sin1 +
          c02 * cos2 + dd02 * sin2 +
          c03 * cos3 + dd03 * sin3 +
          c04 * cos4 + dd04 * sin4 +
          c05 * cos5 + dd05 * sin5 +
          0.5 * c06 * last)

    a1 = (0.5 * c10 +
          c11 * cos1 + dd11 * sin1 +
          c12 * cos2 + dd12 * sin2 +
          c13 * cos3 + dd13 * sin3 +
          c14 * cos4 + dd14 * sin4 +
          c15 * cos5 + dd15 * sin5 +
          0.5 * c16 * last)

    a2 = (0.5 * c20 +
          c21 * cos1 + dd21 * sin1 +
          c22 * cos2 + dd22 * sin2 +
          c23 * cos3 + dd23 * sin3 +
          c24 * cos4 + dd24 * sin4 +
          c25 * cos5 + dd25 * sin5 +
          0.5 * c26 * last)

    # equation 3.4

    dum = a0 + a1 * zeta + a2 * zeta2

    z = np.exp(-cc * zeta)

    xnum = dum * z

    xden = zeta3 + 6.290e-3 * xlm1 + 7.483e-3 * xlm2 + 3.061e-4 * xlm3

    dum = 1.0 / xden
    fphot = xnum * dum

    # equation 3.3

    a0 = 1.0 + 2.045 * xl
    xnum = 0.666 * a0**-2.066

    dum = 1.875e8 * xl + 1.653e8 * xl2 + 8.499e8 * xl3 - 1.604e8 * xl4

    z = 1.0 / dum
    xden = 1.0 + rm * z

    z = 1.0 / xden
    qphot = xnum * z

    # equation 3.2

    sphot = xl5 * fphot

    a1 = sphot
    sphot = rm * a1

    a1 = tfac4 * (1.0 - tfac3 * qphot)
    a3 = sphot
    sphot = a1 * a3

    sphot = max(0.0, sphot)

    # bremsstrahlung neutrino section
    # for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
    #                    n  + n     => n + n + nu + nubar
    #                    n  + p     => n + p + nu + nubar
    # equation 4.3

    den6 = rho * 1.0e-6
    t8 = T * 1.0e-8
    t812 = np.sqrt(t8)
    t832 = t8 * t812
    t82 = t8 * t8
    t83 = t82 * t8
    t85 = t82 * t83
    t86 = t85 * t8
    t8m1 = 1.0 / t8
    t8m2 = t8m1 * t8m1
    t8m3 = t8m2 * t8m1
    t8m5 = t8m3 * t8m2

    tfermi = 5.9302e9 * (np.sqrt(1.0 + 1.018 * (den6 * ye)**twoth) - 1.0)

    # "weak" degenerate electrons only

    if T > 0.3 * tfermi:

        # equation 5.3

        dum = 7.05e6 * t832 + 5.12e4 * t83

        z = 1.0 / dum
        eta = rm * z

        etam1 = 1.0 / eta
        etam2 = etam1 * etam1

        # equation 5.2

        a0 = 23.5 + 6.83e4 * t8m2 + 7.81e8 * t8m5
        xnum = 1.0 / a0

        dum = 1.0 + 1.47 * etam1 + 3.29e-2 * etam2

        c00 = 1.26 * (1.0 + etam1)

        z = 1.0 / dum
        xden = c00 * z

        fbrem = xnum + xden

        # equation 5.9

        a0 = 230.0 + 6.7e5 * t8m2 + 7.66e9 * t8m5

        z = 1.0 + rm * 1.0e-9
        dum = a0 * z

        xnum = 1.0 / dum

        c00 = 7.75e5 * t832 + 247.0 * t8**3.85
        c01 = 4.07 + 0.0240 * t8**1.4
        c02 = 4.59e-5 * t8**-0.110

        z = rho**0.656
        dum = c00 * rmi + c01 + c02 * z

        xden = 1.0 / dum

        gbrem = xnum + xden

        # equation 5.1
        dum = 0.5738 * zbar * ye * t86 * rho

        z = tfac4 * fbrem - tfac5 * gbrem
        sbrem = dum * z

        # liquid metal with c12 parameters (not too different for other elements)
        # equation 5.18 and 5.16

    else:

        u = fac3 * (np.log10(rho) - 3.0)

        # compute the expensive trig functions of equation 5.21 only once
        cos1 = np.cos(u)
        cos2 = np.cos(2.0 * u)
        cos3 = np.cos(3.0 * u)
        cos4 = np.cos(4.0 * u)
        cos5 = np.cos(5.0 * u)

        sin1 = np.sin(u)
        sin2 = np.sin(2.0 * u)
        sin3 = np.sin(3.0 * u)
        sin4 = np.sin(4.0 * u)
        # sin5 = std::sin(5.0*u);

        # equation 5.21

        fb = (0.5 * 0.17946 + 0.00945 * u + 0.34529 -
              0.05821 * cos1 - 0.04969 * sin1 -
              0.01089 * cos2 - 0.01584 * sin2 -
              0.01147 * cos3 - 0.00504 * sin3 -
              0.00656 * cos4 - 0.00281 * sin4 -
              0.00519 * cos5)

        # equation 5.22

        ft = (0.5 * 0.06781 - 0.02342 * u + 0.24819 -
              0.00944 * cos1 - 0.02213 * sin1 -
              0.01289 * cos2 - 0.01136 * sin2 -
              0.00589 * cos3 - 0.00467 * sin3 -
              0.00404 * cos4 - 0.00131 * sin4 -
              0.00330 * cos5)

        # equation 5.23

        gb = (0.5 * 0.00766 - 0.01259 * u + 0.07917 -
              0.00710 * cos1 + 0.02300 * sin1 -
              0.00028 * cos2 - 0.01078 * sin2 +
              0.00232 * cos3 + 0.00118 * sin3 +
              0.00044 * cos4 - 0.00089 * sin4 +
              0.00158 * cos5)

        # equation 5.24

        gt = (-0.5 * 0.00769 - 0.00829 * u + 0.05211 +
              0.00356 * cos1 + 0.01052 * sin1 -
              0.00184 * cos2 - 0.00354 * sin2 +
              0.00146 * cos3 - 0.00014 * sin3 +
              0.00031 * cos4 - 0.00018 * sin4 +
              0.00069 * cos5)

        dum = 2.275e-1 * zbar * zbar * t8m1 * (den6 * abari)**oneth

        gm1 = 1.0 / dum
        gm13 = gm1**oneth
        gm23 = gm13 * gm13

        # equation 5.25 and 5.26

        v = -0.05483 - 0.01946 * gm13 + 1.86310 * gm23 - 0.78873 * gm1

        w = -0.06711 + 0.06859 * gm13 + 1.74360 * gm23 - 0.74498 * gm1

        # equation 5.19 and 5.20

        fliq = v * fb + (1.0 - v) * ft

        gliq = w * gb + (1.0 - w) * gt

        # equation 5.17

        dum = 0.5738 * zbar * ye * t86 * rho

        z = tfac4 * fliq - tfac5 * gliq
        sbrem = dum * z

    # recombination neutrino section
    # for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
    # equation 6.11 solved for nu

    xnum = 1.10520e8 * rho * ye / (T * np.sqrt(T))

    # the chemical potential

    nu = ifermi12(xnum)

    nu2 = nu * nu
    nu3 = nu2 * nu

    # table 12

    if 0.0 > nu >= -20.0:
        a1 = 1.51e-2
        a2 = 2.42e-1
        a3 = 1.21
        b = 3.71e-2
        c = 9.06e-1
        d = 9.28e-1
        f1 = 0.0
        f2 = 0.0
        f3 = 0.0

    elif 0.0 <= nu <= 10.0:
        a1 = 1.23e-2
        a2 = 2.66e-1
        a3 = 1.30
        b = 1.17e-1
        c = 8.97e-1
        d = 1.77e-1
        f1 = -1.20e-2
        f2 = 2.29e-2
        f3 = -1.04e-3

    # equation 6.7, 6.13 and 6.14

    sreco = 0.0

    if -20.0 <= nu <= 10.0:

        zeta = 1.579e5 * zbar * zbar * tempi

        # pylint: disable-next=possibly-used-before-assignment
        c00 = 1.0 / (1.0 + f1 * nu + f2 * nu2 + f3 * nu3)
        dum = zeta * c00

        z = 1.0 / dum
        dd00 = dum**-2.25
        dd01 = dum**-4.55
        c00 = a1 * z + a2 * dd00 + a3 * dd01

        z = np.exp(c * nu)
        # pylint: disable-next=possibly-used-before-assignment
        dd00 = b * z * (1.0 + d * dum)
        gum = 1.0 + dd00

        z = np.exp(nu)
        a1 = 1.0 / gum

        bigj = c00 * z * a1

        # equation 6.5
        z = np.exp(zeta + nu)
        dum = 1.0 + z
        a1 = 1.0 / dum

        sreco = tfac6 * 2.649e-18 * ye * zbar**13.0 * rho * bigj * a1

    # convert from erg/cm^3/s to erg/g/s
    # comment these out to duplicate the itoh et al plots

    spair = spair * deni
    splas = splas * deni
    sphot = sphot * deni
    sbrem = sbrem * deni
    sreco = sreco * deni

    # the total neutrino loss rate
    snu = splas + spair + sphot + sbrem + sreco

    if full_output:
        comps = NeutrinoComponents()
        comps.splas = splas
        comps.spair = spair
        comps.sphot = sphot
        comps.sbrem = sbrem
        comps.sreco = sreco

        return snu, comps

    return snu
