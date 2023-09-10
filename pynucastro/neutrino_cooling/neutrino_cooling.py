
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

    a1 = np.array([1.999266880833e4,
                   5.702479099336e3,
                   6.610132843877e2,
                   3.818838129486e1,
                   1.0e0])

    b1 = np.array([1.771804140488e4,
                   -2.014785161019e3,
                   9.130355392717e1,
                   -1.670718177489e0])

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

        rn  = f + a1[m1]

        for i in reversed(range(m1)):
            rn = rn * f + a1[i]

        # now we do the denominator in Eq. 4.  None of the coefficients
        # are 1, so we loop over all

        den = b1[k1+1]

        for i in reversed(range(k1+1)):
            den = den * f + b1[i]

        # Eq. 6 of Antia

        ifermi12r = np.log(f * rn / den)

    else:

        # this construction is the same as above, but using the
        # second set of coefficients

        an = 0.5
        ff = 1.0 / f**(1.0 / (1.0 + an))

        rn = ff + a2[m2];

        for i in reversed(range(m2)):
            rn = rn * ff + a2[i]

        den = b2[k2+1]

        for i in reversed(range(k2+1)):
            den = den * ff + b2[i]

        ifermi12r = rn / (den * ff)

    return ifermi12r


def zfermim12(x):
    """apply a rational function expansion to get the fermi-dirac
    integral of order -1/2 evaluated at x. maximum error is 1.23e-12_rt.
    reference: antia apjs 84,101 1993"""

    # coefficients of the expansion from Table 2 of Antia

    m1 = 6
    k1 = 6
    m2 = 10
    k2 = 10

    a1 = np.array([1.71446374704454e7,
                   3.88148302324068e7,
                   3.16743385304962e7,
                   1.14587609192151e7,
                   1.83696370756153e6,
                   1.14980998186874e5,
                   1.98276889924768e3,
                   1.0])

    b1 = np.array([9.67282587452899e6,
                   2.87386436731785e7,
                   3.26070130734158e7,
                   1.77657027846367e7,
                   4.81648022267831e6,
                   6.13709569333207e5,
                   3.13595854332114e4,
                   4.35061725080755e2])

    a2 = np.array([-4.46620341924942e-15,
                   -1.58654991146236e-12,
                   -4.44467627042232e-10,
                   -6.84738791621745e-8,
                   -6.64932238528105e-6,
                   -3.69976170193942e-4,
                   -1.12295393687006e-2,
                   -1.60926102124442e-1,
                   -8.52408612877447e-1,
                   -7.45519953763928e-1,
                   2.98435207466372e0,
                   1.0])

    b2 = np.array([-2.23310170962369e-15,
                   -7.94193282071464e-13,
                   -2.22564376956228e-10,
                   -3.43299431079845e-8,
                   -3.33919612678907e-6,
                   -1.86432212187088e-4,
                   -5.69764436880529e-3,
                   -8.34904593067194e-2,
                   -4.78770844009440e-1,
                   -4.99759250374148e-1,
                   1.86795964993052e0,
                   4.16485970495288e-1])

    if x < 2.0:

       xx = np.exp(x)

       rn = xx + a1[m1];
       for i in reversed(range(m1)):
           rn = rn * xx + a1[i];

       den = b1[k1+1]
       for i in reversed(range(k1+1)):
           den = den * xx + b1[i]

       zfermim12r = xx * rn/den;

    else:

       xx = 1.0 / (x*x)

       rn = xx + a2[m2]
       for i in reversed(range(m2)):
           rn = rn * xx + a2[i]

       den = b2[k2+1]
       for i in reversed(range(k2+1)):
           den = den * xx + b2[i]

       zfermim12r = np.sqrt(x) * rn / den;

    return zfermim12r


def sneut5(temp, den, abar, zbar):
    """compute thermal neutrino losses from the analytic fits of
    Itoh et al. ApJS 102, 411, 1996"""

    # input:
    # temp = temperature
    # den  = density
    # abar = mean atomic weight
    # zbar = mean charge

    # output:
    # snu    = total neutrino loss rate in erg/g/sec

    # numerical constants

    fac1 = 5.0 * np.pi / 3.0
    fac2 = 10.0 * np.pi;
    fac3 = np.pi / 5.0
    oneth = 1.0 / 3.0
    twoth = 2.0 /3.0
    con1 = 1.0 / 5.9302e0
    sixth = 1.0 / 6.0
    iln10 = 4.342944819032518e-1

    # theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
    # xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
    # change theta and xnufam if need be, and the changes will automatically
    # propagate through the routine. cv and ca are the vector and axial currents.

    constexpr Real theta  = 0.2319e0_rt;
    constexpr Real xnufam = 3.0e0_rt;
    constexpr Real cv     = 0.5e0_rt + 2.0e0_rt * theta;
    constexpr Real cvp    = 1.0e0_rt - cv;
    constexpr Real ca     = 0.5e0_rt;
    constexpr Real cap    = 1.0e0_rt - ca;
    constexpr Real tfac1  = cv*cv + ca*ca + (xnufam-1.0e0_rt) * (cvp*cvp+cap*cap);
    constexpr Real tfac2  = cv*cv - ca*ca + (xnufam-1.0e0_rt) * (cvp*cvp - cap*cap);
    constexpr Real tfac3  = tfac2/tfac1;
    constexpr Real tfac4  = 0.5e0_rt * tfac1;
    constexpr Real tfac5  = 0.5e0_rt * tfac2;
    constexpr Real tfac6  = cv*cv + 1.5e0_rt*ca*ca + (xnufam - 1.0e0_rt)*(cvp*cvp + 1.5e0_rt*cap*cap);

    # initialize
    Real spair{0.0e0_rt};
    Real spairdt{0.0e0_rt};
    Real spairda{0.0e0_rt};
    Real spairdz{0.0e0_rt};

    Real splas{0.0e0_rt};
    Real splasdt{0.0e0_rt};
    Real splasda{0.0e0_rt};
    Real splasdz{0.0e0_rt};

    Real sphot{0.0e0_rt};
    Real sphotdt{0.0e0_rt};
    Real sphotda{0.0e0_rt};
    Real sphotdz{0.0e0_rt};

    Real sbrem{0.0e0_rt};
    Real sbremdt{0.0e0_rt};
    Real sbremda{0.0e0_rt};
    Real sbremdz{0.0e0_rt};

    Real sreco{0.0e0_rt};
    Real srecodt{0.0e0_rt};
    Real srecoda{0.0e0_rt};
    Real srecodz{0.0e0_rt};

    snu     = 0.0e0_rt;
    dsnudt  = 0.0e0_rt;

    if (temp < 1.0e7_rt) return;

    # to avoid lots of divisions
    deni  = 1.0e0_rt/den;
    tempi = 1.0e0_rt/temp;
    abari = 1.0e0_rt/abar;
    zbari = 1.0e0_rt/zbar;

    # some composition variables
    ye    = zbar*abari;

    # some frequent factors
    t9     = temp * 1.0e-9_rt;
    xl     = t9 * con1;
    xldt   = 1.0e-9_rt * con1;
    xlp5   = sqrt(xl);
    xl2    = xl*xl;
    xl3    = xl2*xl;
    xl4    = xl3*xl;
    xl5    = xl4*xl;
    xl6    = xl5*xl;
    xl7    = xl6*xl;
    xl8    = xl7*xl;
    xl9    = xl8*xl;
    xlmp5  = 1.0e0_rt/xlp5;
    xlm1   = 1.0e0_rt/xl;
    xlm2   = xlm1*xlm1;
    xlm3   = xlm1*xlm2;
    xlm4   = xlm1*xlm3;

    rm     = den*ye;
    rmda   = -rm*abari;
    rmdz   = den*abari;
    rmi    = 1.0e0_rt/rm;
;
    a0     = rm * 1.0e-9_rt;
    a1     = std::pow(a0, oneth);
    zeta   = a1 * xlm1;
    zetadt = -a1 * xlm2 * xldt;
    a2     = oneth * a1*rmi * xlm1;
    zetada = a2 * rmda;
    zetadz = a2 * rmdz;
;
    zeta2 = zeta * zeta;
    zeta3 = zeta2 * zeta;

    # pair neutrino section
    # for reactions like e+ + e- => nu_e + nubar_e

    # equation 2.8
    gl   = 1.0e0_rt - 13.04e0_rt*xl2 +133.5e0_rt*xl4 +1534.0e0_rt*xl6 +918.6e0_rt*xl8;
    gldt = xldt*(-26.08e0_rt*xl +534.0e0_rt*xl3 +9204.0e0_rt*xl5 +7348.8e0_rt*xl7);

    # equation 2.7

    a1     = 6.002e19_rt + 2.084e20_rt*zeta + 1.872e21_rt*zeta2;
    a2     = 2.084e20_rt + 2.0e0_rt*1.872e21_rt*zeta;

    if (t9 < 10.0_rt) {
       b1     = exp(-5.5924e0_rt*zeta);
       b2     = -b1*5.5924e0_rt;
    } else {
       b1     = exp(-4.9924e0_rt*zeta);
       b2     = -b1*4.9924e0_rt;
    }

    xnum   = a1 * b1;
    c      = a2*b1 + a1*b2;
    xnumdt = c*zetadt;
    xnumda = c*zetada;
    xnumdz = c*zetadz;

    if (t9 < 10.0_rt) {
       a1   = 9.383e-1_rt*xlm1 - 4.141e-1_rt*xlm2 + 5.829e-2_rt*xlm3;
       a2   = -9.383e-1_rt*xlm2 + 2.0e0_rt*4.141e-1_rt*xlm3 - 3.0e0_rt*5.829e-2_rt*xlm4;
    } else {
       a1   = 1.2383e0_rt*xlm1 - 8.141e-1_rt*xlm2;
       a2   = -1.2383e0_rt*xlm2 + 2.0e0_rt*8.141e-1_rt*xlm3;
    }

    b1   = 3.0e0_rt*zeta2;

    xden   = zeta3 + a1;
    xdendt = b1*zetadt + a2*xldt;
    xdenda = b1*zetada;
    xdendz = b1*zetadz;

    a1      = 1.0e0_rt/xden;
    fpair   = xnum*a1;
    fpairdt = (xnumdt - fpair*xdendt)*a1;
    fpairda = (xnumda - fpair*xdenda)*a1;
    fpairdz = (xnumdz - fpair*xdendz)*a1;

    # equation 2.6
    a1     = 10.7480e0_rt*xl2 + 0.3967e0_rt*xlp5 + 1.005e0_rt;
    a2     = xldt*(2.0e0_rt*10.7480e0_rt*xl + 0.5e0_rt*0.3967e0_rt*xlmp5);
    xnum   = 1.0e0_rt/a1;
    xnumdt = -xnum*xnum*a2;

    a1     = 7.692e7_rt*xl3 + 9.715e6_rt*xlp5;
    a2     = xldt*(3.0e0_rt*7.692e7_rt*xl2 + 0.5e0_rt*9.715e6_rt*xlmp5);

    c      = 1.0e0_rt/a1;
    b1     = 1.0e0_rt + rm*c;

    xden   = std::pow(b1, -0.3e0_rt);

    d      = -0.3e0_rt*xden/b1;
    xdendt = -d*rm*c*c*a2;
    xdenda = d*rmda*c;
    xdendz = d*rmdz*c;

    qpair   = xnum*xden;
    qpairdt = xnumdt*xden + xnum*xdendt;
    qpairda = xnum*xdenda;
    qpairdz = xnum*xdendz;

    # equation 2.5
    a1    = std::exp(-2.0e0_rt*xlm1);
    a2    = a1*2.0e0_rt*xlm2*xldt;

    spair   = a1*fpair;
    spairdt = a2*fpair + a1*fpairdt;
    spairda = a1*fpairda;
    spairdz = a1*fpairdz;

    a1      = spair;
    spair   = gl*a1;
    spairdt = gl*spairdt + gldt*a1;
    spairda = gl*spairda;
    spairdz = gl*spairdz;

    a1      = tfac4*(1.0e0_rt + tfac3 * qpair);
    a2      = tfac4*tfac3;

    a3      = spair;
    spair   = a1*a3;
    spairdt = a1*spairdt + a2*qpairdt*a3;
    spairda = a1*spairda + a2*qpairda*a3;
    spairdz = a1*spairdz + a2*qpairdz*a3;

    # plasma neutrino section
    # for collective reactions like gamma_plasmon => nu_e + nubar_e
    # equation 4.6

    a1   = 1.019e-6_rt*rm;
    a2   = std::pow(a1, twoth);
    a3   = twoth*a2/a1;

    b1   = std::sqrt(1.0e0_rt + a2);
    b2   = 1.0e0_rt/b1;

    c00  = 1.0e0_rt/(temp*temp*b1);

    gl2   = 1.1095e11_rt * rm * c00;

    gl2dt = -2.0e0_rt*gl2*tempi;
    d     = rm*c00*b2*0.5e0_rt*b2*a3*1.019e-6_rt;
    gl2da = 1.1095e11_rt * (rmda*c00  - d*rmda);
    gl2dz = 1.1095e11_rt * (rmdz*c00  - d*rmdz);

    gl    = std::sqrt(gl2);
    gl12  = std::sqrt(gl);
    gl32  = gl * gl12;
    gl72  = gl2 * gl32;
    gl6   = gl2 * gl2 * gl2;

    # equation 4.7
    ft   = 2.4e0_rt + 0.6e0_rt*gl12 + 0.51e0_rt*gl + 1.25e0_rt*gl32;
    gum  = 1.0e0_rt/gl2;
    a1   =(0.25e0_rt*0.6e0_rt*gl12 +0.5e0_rt*0.51e0_rt*gl +0.75e0_rt*1.25e0_rt*gl32)*gum;
    ftdt = a1*gl2dt;
    ftda = a1*gl2da;
    ftdz = a1*gl2dz;

    # equation 4.8
    a1   = 8.6e0_rt*gl2 + 1.35e0_rt*gl72;
    a2   = 8.6e0_rt + 1.75e0_rt*1.35e0_rt*gl72*gum;

    b1   = 225.0e0_rt - 17.0e0_rt*gl + gl2;
    b2   = -0.5e0_rt*17.0e0_rt*gl*gum + 1.0e0_rt;

    c    = 1.0e0_rt/b1;
    fl   = a1*c;

    d    = (a2 - fl*b2)*c;
    fldt = d*gl2dt;
    flda = d*gl2da;
    fldz = d*gl2dz;

    # equation 4.9 and 4.10
    cc   = std::log10(2.0e0_rt*rm);
    xlnt = std::log10(temp);

    xnum   = sixth * (17.5e0_rt + cc - 3.0e0_rt*xlnt);
    xnumdt = -iln10*0.5e0_rt*tempi;
    a2     = iln10*sixth*rmi;
    xnumda = a2*rmda;
    xnumdz = a2*rmdz;

    xden   = sixth * (-24.5e0_rt + cc + 3.0e0_rt*xlnt);
    xdendt = iln10*0.5e0_rt*tempi;
    xdenda = a2*rmda;
    xdendz = a2*rmdz;

    # equation 4.11
    if (std::abs(xnum) > 0.7e0_rt || xden < 0.0e0_rt) {

       fxy   = 1.0e0_rt;
       fxydt = 0.0e0_rt;
       fxydz = 0.0e0_rt;
       fxyda = 0.0e0_rt;

    } else {

       a1  = 0.39e0_rt - 1.25e0_rt*xnum - 0.35e0_rt*std::sin(4.5e0_rt*xnum);
       a2  = -1.25e0_rt - 4.5e0_rt*0.35e0_rt*std::cos(4.5e0_rt*xnum);

       b1  = 0.3e0_rt * std::exp(-std::pow(4.5e0_rt*xnum + 0.9e0_rt, 2.0_rt));
       b2  = -b1*2.0e0_rt*(4.5e0_rt*xnum + 0.9e0_rt)*4.5e0_rt;

       c   = amrex::min(0.0e0_rt, xden - 1.6e0_rt + 1.25e0_rt*xnum);
       if (c == 0.0_rt) {
          dumdt = 0.0e0_rt;
          dumda = 0.0e0_rt;
          dumdz = 0.0e0_rt;
       } else {
          dumdt = xdendt + 1.25e0_rt*xnumdt;
          dumda = xdenda + 1.25e0_rt*xnumda;
          dumdz = xdendz + 1.25e0_rt*xnumdz;
       }

       d   = 0.57e0_rt - 0.25e0_rt*xnum;
       a3  = c/d;
       c00 = std::exp(-a3*a3);

       f1  = -c00*2.0e0_rt*a3/d;
       c01 = f1*(dumdt + a3*0.25e0_rt*xnumdt);
       c03 = f1*(dumda + a3*0.25e0_rt*xnumda);
       c04 = f1*(dumdz + a3*0.25e0_rt*xnumdz);

       fxy   = 1.05e0_rt + (a1 - b1)*c00;
       fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01;
       fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03;
       fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04;

    }

    # equation 4.1 and 4.5
    splas   = (ft + fl) * fxy;
    splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt;
    splasda = (ftda + flda)*fxy + (ft+fl)*fxyda;
    splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz;

    a2      = std::exp(-gl);
    a3      = -0.5e0_rt*a2*gl*gum;

    a1      = splas;
    splas   = a2*a1;
    splasdt = a2*splasdt + a3*gl2dt*a1;
    splasda = a2*splasda + a3*gl2da*a1;
    splasdz = a2*splasdz + a3*gl2dz*a1;

    a2      = gl6;
    a3      = 3.0e0_rt*gl6*gum;

    a1      = splas;
    splas   = a2*a1;
    splasdt = a2*splasdt + a3*gl2dt*a1;
    splasda = a2*splasda + a3*gl2da*a1;
    splasdz = a2*splasdz + a3*gl2dz*a1;

    a2      = 0.93153e0_rt * 3.0e21_rt * xl9;
    a3      = 0.93153e0_rt * 3.0e21_rt * 9.0e0_rt*xl8*xldt;

    a1      = splas;
    splas   = a2*a1;
    splasdt = a2*splasdt + a3*a1;
    splasda = a2*splasda;
    splasdz = a2*splasdz;

    # photoneutrino process section
    # for reactions like e- + gamma => e- + nu_e + nubar_e
    #                    e+ + gamma => e+ + nu_e + nubar_e
    # equation 3.8 for tau, equation 3.6 for cc,
    # and table 2 written out for speed
    if (temp < 1.0e8_rt) {

        # note: we already bailed above for T < 1.e7, so this is really 1.e7 <= T < 1.e8

       tau  =  std::log10(temp * 1.0e-7_rt);
       cc   =  0.5654e0_rt + tau;
       c00  =  1.008e11_rt;
       c01  =  0.0e0_rt;
       c02  =  0.0e0_rt;
       c03  =  0.0e0_rt;
       c04  =  0.0e0_rt;
       c05  =  0.0e0_rt;
       c06  =  0.0e0_rt;
       c10  =  8.156e10_rt;
       c11  =  9.728e8_rt;
       c12  = -3.806e9_rt;
       c13  = -4.384e9_rt;
       c14  = -5.774e9_rt;
       c15  = -5.249e9_rt;
       c16  = -5.153e9_rt;
       c20  =  1.067e11_rt;
       c21  = -9.782e9_rt;
       c22  = -7.193e9_rt;
       c23  = -6.936e9_rt;
       c24  = -6.893e9_rt;
       c25  = -7.041e9_rt;
       c26  = -7.193e9_rt;
       dd01 =  0.0e0_rt;
       dd02 =  0.0e0_rt;
       dd03 =  0.0e0_rt;
       dd04 =  0.0e0_rt;
       dd05 =  0.0e0_rt;
       dd11 = -1.879e10_rt;
       dd12 = -9.667e9_rt;
       dd13 = -5.602e9_rt;
       dd14 = -3.370e9_rt;
       dd15 = -1.825e9_rt;
       dd21 = -2.919e10_rt;
       dd22 = -1.185e10_rt;
       dd23 = -7.270e9_rt;
       dd24 = -4.222e9_rt;
       dd25 = -1.560e9_rt;

    } else if (temp >= 1.0e8_rt && temp < 1.0e9_rt) {

       tau  =  std::log10(temp * 1.0e-8_rt);
       cc   =  1.5654e0_rt;
       c00  =  9.889e10_rt;
       c01  = -4.524e8_rt;
       c02  = -6.088e6_rt;
       c03  =  4.269e7_rt;
       c04  =  5.172e7_rt;
       c05  =  4.910e7_rt;
       c06  =  4.388e7_rt;
       c10  =  1.813e11_rt;
       c11  = -7.556e9_rt;
       c12  = -3.304e9_rt;
       c13  = -1.031e9_rt;
       c14  = -1.764e9_rt;
       c15  = -1.851e9_rt;
       c16  = -1.928e9_rt;
       c20  =  9.750e10_rt;
       c21  =  3.484e10_rt;
       c22  =  5.199e9_rt;
       c23  = -1.695e9_rt;
       c24  = -2.865e9_rt;
       c25  = -3.395e9_rt;
       c26  = -3.418e9_rt;
       dd01 = -1.135e8_rt;
       dd02 =  1.256e8_rt;
       dd03 =  5.149e7_rt;
       dd04 =  3.436e7_rt;
       dd05 =  1.005e7_rt;
       dd11 =  1.652e9_rt;
       dd12 = -3.119e9_rt;
       dd13 = -1.839e9_rt;
       dd14 = -1.458e9_rt;
       dd15 = -8.956e8_rt;
       dd21 = -1.548e10_rt;
       dd22 = -9.338e9_rt;
       dd23 = -5.899e9_rt;
       dd24 = -3.035e9_rt;
       dd25 = -1.598e9_rt;

    } else {

        # T > 1.e9

       tau  =  std::log10(t9);
       cc   =  1.5654e0_rt;
       c00  =  9.581e10_rt;
       c01  =  4.107e8_rt;
       c02  =  2.305e8_rt;
       c03  =  2.236e8_rt;
       c04  =  1.580e8_rt;
       c05  =  2.165e8_rt;
       c06  =  1.721e8_rt;
       c10  =  1.459e12_rt;
       c11  =  1.314e11_rt;
       c12  = -1.169e11_rt;
       c13  = -1.765e11_rt;
       c14  = -1.867e11_rt;
       c15  = -1.983e11_rt;
       c16  = -1.896e11_rt;
       c20  =  2.424e11_rt;
       c21  = -3.669e9_rt;
       c22  = -8.691e9_rt;
       c23  = -7.967e9_rt;
       c24  = -7.932e9_rt;
       c25  = -7.987e9_rt;
       c26  = -8.333e9_rt;
       dd01 =  4.724e8_rt;
       dd02 =  2.976e8_rt;
       dd03 =  2.242e8_rt;
       dd04 =  7.937e7_rt;
       dd05 =  4.859e7_rt;
       dd11 = -7.094e11_rt;
       dd12 = -3.697e11_rt;
       dd13 = -2.189e11_rt;
       dd14 = -1.273e11_rt;
       dd15 = -5.705e10_rt;
       dd21 = -2.254e10_rt;
       dd22 = -1.551e10_rt;
       dd23 = -7.793e9_rt;
       dd24 = -4.489e9_rt;
       dd25 = -2.185e9_rt;

    }

    taudt = iln10*tempi;

    # equation 3.7, compute the expensive trig functions only one time
    cos1 = std::cos(fac1*tau);
    cos2 = std::cos(fac1*2.0e0_rt*tau);
    cos3 = std::cos(fac1*3.0e0_rt*tau);
    cos4 = std::cos(fac1*4.0e0_rt*tau);
    cos5 = std::cos(fac1*5.0e0_rt*tau);
    last = std::cos(fac2*tau);

    sin1 = std::sin(fac1*tau);
    sin2 = std::sin(fac1*2.0e0_rt*tau);
    sin3 = std::sin(fac1*3.0e0_rt*tau);
    sin4 = std::sin(fac1*4.0e0_rt*tau);
    sin5 = std::sin(fac1*5.0e0_rt*tau);
    xast = std::sin(fac2*tau);

    a0 = 0.5e0_rt*c00
         + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2
         + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4
         + c05*cos5 + dd05*sin5 + 0.5e0_rt*c06*last;

    f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0e0_rt
         + dd02*cos2*2.0e0_rt - c03*sin3*3.0e0_rt + dd03*cos3*3.0e0_rt
         - c04*sin4*4.0e0_rt + dd04*cos4*4.0e0_rt
         - c05*sin5*5.0e0_rt + dd05*cos5*5.0e0_rt)
         - 0.5e0_rt*c06*xast*fac2*taudt;

    a1 = 0.5e0_rt*c10
         + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2
         + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4
         + c15*cos5 + dd15*sin5 + 0.5e0_rt*c16*last;

    f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0e0_rt
         + dd12*cos2*2.0e0_rt - c13*sin3*3.0e0_rt + dd13*cos3*3.0e0_rt
         - c14*sin4*4.0e0_rt + dd14*cos4*4.0e0_rt - c15*sin5*5.0e0_rt
         + dd15*cos5*5.0e0_rt) - 0.5e0_rt*c16*xast*fac2*taudt;

    a2 = 0.5e0_rt*c20
         + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2
         + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4
         + c25*cos5 + dd25*sin5 + 0.5e0_rt*c26*last;

    f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0e0_rt
         + dd22*cos2*2.0e0_rt - c23*sin3*3.0e0_rt + dd23*cos3*3.0e0_rt
         - c24*sin4*4.0e0_rt + dd24*cos4*4.0e0_rt - c25*sin5*5.0e0_rt
         + dd25*cos5*5.0e0_rt) - 0.5e0_rt*c26*xast*fac2*taudt;

    # equation 3.4
    dum   = a0 + a1*zeta + a2*zeta2;
    dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0e0_rt*a2*zeta*zetadt;

    z      = std::exp(-cc*zeta);

    xnum   = dum*z;
    xnumdt = dumdt*z - dum*z*cc*zetadt;

    xden   = zeta3 + 6.290e-3_rt*xlm1 + 7.483e-3_rt*xlm2 + 3.061e-4_rt*xlm3;

    dum    = 3.0e0_rt*zeta2;
    xdendt = dum*zetadt - xldt*(6.290e-3_rt*xlm2
         + 2.0e0_rt*7.483e-3_rt*xlm3 + 3.0e0_rt*3.061e-4_rt*xlm4);

    dum      = 1.0e0_rt/xden;
    fphot   = xnum*dum;
    fphotdt = (xnumdt - fphot*xdendt)*dum;

    # equation 3.3
    a0     = 1.0e0_rt + 2.045e0_rt * xl;
    xnum   = 0.666e0_rt*std::pow(a0, -2.066e0_rt);
    xnumdt = -2.066e0_rt*xnum/a0 * 2.045e0_rt*xldt;

    dum    = 1.875e8_rt*xl + 1.653e8_rt*xl2 + 8.499e8_rt*xl3 - 1.604e8_rt*xl4;
    dumdt  = xldt*(1.875e8_rt + 2.0e0_rt*1.653e8_rt*xl + 3.0e0_rt*8.499e8_rt*xl2
             - 4.0e0_rt*1.604e8_rt*xl3);

    z      = 1.0e0_rt/dum;
    xden   = 1.0e0_rt + rm*z;
    xdendt =  -rm*z*z*dumdt;

    z      = 1.0e0_rt/xden;
    qphot = xnum*z;
    qphotdt = (xnumdt - qphot*xdendt)*z;
    dum      = -qphot*z;

    # equation 3.2
    sphot   = xl5 * fphot;
    sphotdt = 5.0e0_rt*xl4*xldt*fphot + xl5*fphotdt;

    a1      = sphot;
    sphot   = rm*a1;
    sphotdt = rm*sphotdt;

    a1      = tfac4*(1.0e0_rt - tfac3 * qphot);
    a2      = -tfac4*tfac3;

    a3      = sphot;
    sphot   = a1*a3;
    sphotdt = a1*sphotdt + a2*qphotdt*a3;

    if (sphot <= 0.0_rt) {
       sphot   = 0.0e0_rt;
       sphotdt = 0.0e0_rt;
    }

    # bremsstrahlung neutrino section
    # for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
    #                    n  + n     => n + n + nu + nubar
    #                    n  + p     => n + p + nu + nubar
    # equation 4.3

    den6   = den * 1.0e-6_rt;
    t8     = temp * 1.0e-8_rt;
    t812   = std::sqrt(t8);
    t832   = t8 * t812;
    t82    = t8*t8;
    t83    = t82*t8;
    t85    = t82*t83;
    t86    = t85*t8;
    t8m1   = 1.0e0_rt/t8;
    t8m2   = t8m1*t8m1;
    t8m3   = t8m2*t8m1;
    t8m5   = t8m3*t8m2;
    t8m6   = t8m5*t8m1;


    tfermi = 5.9302e9_rt*(std::sqrt(1.0e0_rt+1.018e0_rt*std::pow(den6*ye, twoth))-1.0e0_rt);

    # "weak" degenerate electrons only
    if (temp > 0.3e0_rt * tfermi) {

       # equation 5.3
       dum   = 7.05e6_rt * t832 + 5.12e4_rt * t83;
       dumdt = (1.5e0_rt*7.05e6_rt*t812 + 3.0e0_rt*5.12e4_rt*t82)*1.0e-8_rt;

       z     = 1.0e0_rt/dum;
       eta   = rm*z;
       etadt = -rm*z*z*dumdt;
       etada = rmda*z;
       etadz = rmdz*z;

       etam1 = 1.0e0_rt/eta;
       etam2 = etam1 * etam1;
       etam3 = etam2 * etam1;

       # equation 5.2
       a0    = 23.5e0_rt + 6.83e4_rt*t8m2 + 7.81e8_rt*t8m5;
       f0    = (-2.0e0_rt*6.83e4_rt*t8m3 - 5.0e0_rt*7.81e8_rt*t8m6)*1.0e-8_rt;
       xnum  = 1.0e0_rt/a0;

       dum   = 1.0e0_rt + 1.47e0_rt*etam1 + 3.29e-2_rt*etam2;
       z     = -1.47e0_rt*etam2 - 2.0e0_rt*3.29e-2_rt*etam3;
       dumdt = z*etadt;

       c00   = 1.26e0_rt * (1.0e0_rt+etam1);
       z     = -1.26e0_rt*etam2;
       c01   = z*etadt;
       c03   = z*etada;
       c04   = z*etadz;

       z      = 1.0e0_rt/dum;
       xden   = c00*z;
       xdendt = (c01 - xden*dumdt)*z;

       fbrem   = xnum + xden;
       fbremdt = -xnum*xnum*f0 + xdendt;

       # equation 5.9
       a0    = 230.0e0_rt + 6.7e5_rt*t8m2 + 7.66e9_rt*t8m5;
       f0    = (-2.0e0_rt*6.7e5_rt*t8m3 - 5.0e0_rt*7.66e9_rt*t8m6)*1.0e-8_rt;

       z     = 1.0e0_rt + rm*1.0e-9_rt;
       dum   = a0*z;
       dumdt = f0*z;
       z     = a0*1.0e-9_rt;

       xnum   = 1.0e0_rt/dum;
       z      = -xnum*xnum;
       xnumdt = z*dumdt;

       c00   = 7.75e5_rt*t832 + 247.0e0_rt * std::pow(t8, 3.85e0_rt);
       dd00  = (1.5e0_rt*7.75e5_rt*t812 + 3.85e0_rt*247.0e0_rt*std::pow(t8, 2.85e0_rt))*1.0e-8_rt;

       c01   = 4.07e0_rt + 0.0240e0_rt * std::pow(t8, 1.4e0_rt);
       dd01  = 1.4e0_rt*0.0240e0_rt * std::pow(t8, 0.4e0_rt)*1.0e-8_rt;

       c02   = 4.59e-5_rt * std::pow(t8, -0.110e0_rt);
       dd02  = -0.11e0_rt*4.59e-5_rt * std::pow(t8, -1.11e0_rt)*1.0e-8_rt;

       z     = std::pow(den, 0.656e0_rt);
       dum   = c00*rmi  + c01  + c02*z;
       dumdt = dd00*rmi + dd01 + dd02*z;
       z     = -c00*rmi*rmi;

       xden  = 1.0e0_rt/dum;
       z      = -xden*xden;
       xdendt = z*dumdt;

       gbrem   = xnum + xden;
       gbremdt = xnumdt + xdendt;

       # equation 5.1
       dum    = 0.5738e0_rt*zbar*ye*t86*den;
       dumdt  = 0.5738e0_rt*zbar*ye*6.0e0_rt*t85*den*1.0e-8_rt;

       z       = tfac4*fbrem - tfac5*gbrem;
       sbrem   = dum * z;
       sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt);

       # liquid metal with c12 parameters (not too different for other elements)
       # equation 5.18 and 5.16

    } else {

       u     = fac3 * (std::log10(den) - 3.0e0_rt);
       #a0    = iln10*fac3*deni;

       # compute the expensive trig functions of equation 5.21 only once
       cos1 = std::cos(u);
       cos2 = std::cos(2.0e0_rt*u);
       cos3 = std::cos(3.0e0_rt*u);
       cos4 = std::cos(4.0e0_rt*u);
       cos5 = std::cos(5.0e0_rt*u);

       sin1 = std::sin(u);
       sin2 = std::sin(2.0e0_rt*u);
       sin3 = std::sin(3.0e0_rt*u);
       sin4 = std::sin(4.0e0_rt*u);
       #sin5 = std::sin(5.0e0_rt*u);

       # equation 5.21
       fb =  0.5e0_rt * 0.17946e0_rt  + 0.00945e0_rt*u + 0.34529e0_rt
            - 0.05821e0_rt*cos1 - 0.04969e0_rt*sin1
            - 0.01089e0_rt*cos2 - 0.01584e0_rt*sin2
            - 0.01147e0_rt*cos3 - 0.00504e0_rt*sin3
            - 0.00656e0_rt*cos4 - 0.00281e0_rt*sin4
            - 0.00519e0_rt*cos5;

       # equation 5.22
       ft =  0.5e0_rt * 0.06781e0_rt - 0.02342e0_rt*u + 0.24819e0_rt
            - 0.00944e0_rt*cos1 - 0.02213e0_rt*sin1
            - 0.01289e0_rt*cos2 - 0.01136e0_rt*sin2
            - 0.00589e0_rt*cos3 - 0.00467e0_rt*sin3
            - 0.00404e0_rt*cos4 - 0.00131e0_rt*sin4
            - 0.00330e0_rt*cos5;

       # equation 5.23
       gb =  0.5e0_rt * 0.00766e0_rt - 0.01259e0_rt*u + 0.07917e0_rt
            - 0.00710e0_rt*cos1 + 0.02300e0_rt*sin1
            - 0.00028e0_rt*cos2 - 0.01078e0_rt*sin2
            + 0.00232e0_rt*cos3 + 0.00118e0_rt*sin3
            + 0.00044e0_rt*cos4 - 0.00089e0_rt*sin4
            + 0.00158e0_rt*cos5;

       # equation 5.24
       gt =  -0.5e0_rt * 0.00769e0_rt  - 0.00829e0_rt*u + 0.05211e0_rt
            + 0.00356e0_rt*cos1 + 0.01052e0_rt*sin1
            - 0.00184e0_rt*cos2 - 0.00354e0_rt*sin2
            + 0.00146e0_rt*cos3 - 0.00014e0_rt*sin3
            + 0.00031e0_rt*cos4 - 0.00018e0_rt*sin4
            + 0.00069e0_rt*cos5;
       dum   = 2.275e-1_rt * zbar * zbar*t8m1 * std::pow(den6*abari, oneth);
       dumdt = -dum*tempi;

       gm1   = 1.0e0_rt/dum;
       gm2   = gm1*gm1;
       gm13  = std::pow(gm1, oneth);
       gm23  = gm13 * gm13;
       gm43  = gm13*gm1;
       gm53  = gm23*gm1;

       # equation 5.25 and 5.26
       v  = -0.05483e0_rt - 0.01946e0_rt*gm13 + 1.86310e0_rt*gm23 - 0.78873e0_rt*gm1;
       a0 = oneth*0.01946e0_rt*gm43 - twoth*1.86310e0_rt*gm53 + 0.78873e0_rt*gm2;

       w  = -0.06711e0_rt + 0.06859e0_rt*gm13 + 1.74360e0_rt*gm23 - 0.74498e0_rt*gm1;
       a1 = -oneth*0.06859e0_rt*gm43 - twoth*1.74360e0_rt*gm53 + 0.74498e0_rt*gm2;

       # equation 5.19 and 5.20
       fliq   = v*fb + (1.0e0_rt - v)*ft;
       fliqdt = a0*dumdt*(fb - ft);

       gliq   = w*gb + (1.0e0_rt - w)*gt;
       gliqdt = a1*dumdt*(gb - gt);

       # equation 5.17
       dum    = 0.5738e0_rt*zbar*ye*t86*den;
       dumdt  = 0.5738e0_rt*zbar*ye*6.0e0_rt*t85*den*1.0e-8_rt;

       z       = tfac4*fliq - tfac5*gliq;
       sbrem   = dum * z;
       sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt);

    }

    # recombination neutrino section
    # for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
    # equation 6.11 solved for nu
    xnum   = 1.10520e8_rt * den * ye /(temp*std::sqrt(temp));
    xnumdt = -1.50e0_rt*xnum*tempi;

    # the chemical potential
    nu   = ifermi12(xnum);

    # a0 is d(nu)/d(xnum)
    a0 = 1.0e0_rt/(0.5e0_rt*zfermim12(nu));
    nudt = a0*xnumdt;

    nu2  = nu * nu;
    nu3  = nu2 * nu;

    # table 12
    if (nu >= -20.0_rt && nu < 0.0_rt) {
       a1 = 1.51e-2_rt;
       a2 = 2.42e-1_rt;
       a3 = 1.21e0_rt;
       b  = 3.71e-2_rt;
       c  = 9.06e-1_rt;
       d  = 9.28e-1_rt;
       f1 = 0.0e0_rt;
       f2 = 0.0e0_rt;
       f3 = 0.0e0_rt;
    } else if (nu >= 0.0_rt && nu <= 10.0_rt) {
       a1 = 1.23e-2_rt;
       a2 = 2.66e-1_rt;
       a3 = 1.30e0_rt;
       b  = 1.17e-1_rt;
       c  = 8.97e-1_rt;
       d  = 1.77e-1_rt;
       f1 = -1.20e-2_rt;
       f2 = 2.29e-2_rt;
       f3 = -1.04e-3_rt;
    }

    # equation 6.7, 6.13 and 6.14
    if (nu >= -20.0_rt &&  nu <= 10.0_rt) {

       zeta   = 1.579e5_rt*zbar*zbar*tempi;
       zetadt = -zeta*tempi;

       c00    = 1.0e0_rt/(1.0e0_rt + f1*nu + f2*nu2 + f3*nu3);
       c01    = f1 + f2*2.0e0_rt*nu + f3*3.0e0_rt*nu2;
       dum    = zeta*c00;
       dumdt  = zetadt*c00 + zeta*c01*nudt;

       z      = 1.0e0_rt/dum;
       dd00   = std::pow(dum, -2.25_rt);
       dd01   = std::pow(dum, -4.55_rt);
       c00    = a1*z + a2*dd00 + a3*dd01;
       c01    = -(a1*z + 2.25_rt*a2*dd00 + 4.55_rt*a3*dd01)*z;

       z      = std::exp(c*nu);
       dd00   = b*z*(1.0e0_rt + d*dum);
       gum    = 1.0e0_rt + dd00;
       gumdt  = dd00*c*nudt + b*z*d*dumdt;

       z   = std::exp(nu);
       a1  = 1.0e0_rt/gum;

       bigj   = c00 * z * a1;
       bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt;

       # equation 6.5
       z     = std::exp(zeta + nu);
       dum   = 1.0e0_rt + z;
       a1    = 1.0e0_rt/dum;
       a2    = 1.0e0_rt/bigj;

       sreco   = tfac6 * 2.649e-18_rt * ye * std::pow(zbar, 13.0_rt) * den * bigj*a1;
       srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1);

    }

    # convert from erg/cm^3/s to erg/g/s
    # comment these out to duplicate the itoh et al plots

    spair   = spair*deni;
    spairdt = spairdt*deni;

    splas   = splas*deni;
    splasdt = splasdt*deni;

    sphot   = sphot*deni;
    sphotdt = sphotdt*deni;

    sbrem   = sbrem*deni;
    sbremdt = sbremdt*deni;

    sreco   = sreco*deni;
    srecodt = srecodt*deni;

    # the total neutrino loss rate
    snu    =  splas + spair + sphot + sbrem + sreco;
    dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt;

}

