cimport numpy as np
import numpy as np

ip = 0
ihe4 = 1
ic12 = 2
ic13 = 3
in13 = 4
in14 = 5
in15 = 6
io14 = 7
io15 = 8
nnuc = 9

A = np.zeros((nnuc), dtype=np.int)

A[ip] = 1
A[ihe4] = 4
A[ic12] = 12
A[ic13] = 13
A[in13] = 13
A[in14] = 14
A[in15] = 15
A[io14] = 14
A[io15] = 15

cdef class Tfactors:
    """ precompute temperature factors for speed """

    cdef public double T9, T9i, T913i, T913, T953, lnT9

    def __cinit__(self, double T):
        """ return the Tfactors object.  Here, T is temperature in Kelvin """
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)


def c12_pg_n13(Tfactors tf):
    # p + c12 --> n13
    cdef double rate = 0.0
    
    # ls09n
    rate += np.exp(  17.1482 + -13.692*tf.T913i + -0.230881*tf.T913
                  + 4.44362*tf.T9 + -3.15898*tf.T953 + -0.666667*tf.lnT9)
    # ls09r
    rate += np.exp(  17.5428 + -3.77849*tf.T9i + -5.10735*tf.T913i + -2.24111*tf.T913
                  + 0.148883*tf.T9 + -1.5*tf.lnT9)
    
    return rate

def c13_pg_n14(Tfactors tf):
    # p + c13 --> n14
    rate = 0.0
    
    # nacrn
    rate += np.exp(  18.5155 + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9637 + -5.78147*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  15.1825 + -13.5543*tf.T9i
                  + -1.5*tf.lnT9)
    
    return rate

def n13_c13(Tfactors tf):
    # n13 --> c13
    cdef double rate = 0.0
    
    # wc12w
    rate += np.exp(  -6.7601)
    
    return rate

def n13_pg_o14(Tfactors tf):
    # p + n13 --> o14
    cdef double rate = 0.0
    
    # lg06n
    rate += np.exp(  18.1356 + -15.1676*tf.T913i + 0.0955166*tf.T913
                  + 3.0659*tf.T9 + -0.507339*tf.T953 + -0.666667*tf.lnT9)
    # lg06r
    rate += np.exp(  10.9971 + -6.12602*tf.T9i + 1.57122*tf.T913i
                  + -1.5*tf.lnT9)
    
    return rate

def n14_pg_o15(Tfactors tf):
    # p + n14 --> o15
    cdef double rate = 0.0
    
    # im05n
    rate += np.exp(  17.01 + -15.193*tf.T913i + -0.161954*tf.T913
                  + -7.52123*tf.T9 + -0.987565*tf.T953 + -0.666667*tf.lnT9)
    # im05r
    rate += np.exp(  6.73578 + -4.891*tf.T9i
                  + 0.0682*tf.lnT9)
    # im05r
    rate += np.exp(  7.65444 + -2.998*tf.T9i
                  + -1.5*tf.lnT9)
    # im05n
    rate += np.exp(  20.1169 + -15.193*tf.T913i + -4.63975*tf.T913
                  + 9.73458*tf.T9 + -9.55051*tf.T953 + 0.333333*tf.lnT9)
    
    return rate

def n15_pa_c12(Tfactors tf):
    # p + n15 --> he4 + c12
    cdef double rate = 0.0
    
    # nacrn
    rate += np.exp(  27.4764 + -15.253*tf.T913i + 1.59318*tf.T913
                  + 2.4479*tf.T9 + -2.19708*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  -6.57522 + -1.1638*tf.T9i + 22.7105*tf.T913
                  + -2.90707*tf.T9 + 0.205754*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  20.8972 + -7.406*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  -4.87347 + -2.02117*tf.T9i + 30.8497*tf.T913
                  + -8.50433*tf.T9 + -1.54426*tf.T953 + -1.5*tf.lnT9)
    
    return rate

def o14_n14(Tfactors tf):
    # o14 --> n14
    cdef double rate = 0.0
    
    # wc12w
    rate += np.exp(  -4.62354)
    
    return rate

def o15_n15(Tfactors tf):
    # o15 --> n15
    cdef double rate = 0.0
    
    # wc12w
    rate += np.exp(  -5.17053)
    
    return rate

def rhs(double t, np.ndarray[double, ndim=1] Y, double rho, double T):

    cdef Tfactors tf = Tfactors(T)

    cdef double lambda_c12_pg_n13 = c12_pg_n13(tf)
    cdef double lambda_c13_pg_n14 = c13_pg_n14(tf)
    cdef double lambda_n13_c13 = n13_c13(tf)
    cdef double lambda_n13_pg_o14 = n13_pg_o14(tf)
    cdef double lambda_n14_pg_o15 = n14_pg_o15(tf)
    cdef double lambda_n15_pa_c12 = n15_pa_c12(tf)
    cdef double lambda_o14_n14 = o14_n14(tf)
    cdef double lambda_o15_n15 = o15_n15(tf)

    cdef np.ndarray[double, ndim=1] dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[ip] = (
       -rho*Y[ic12]*Y[ip]*lambda_c12_pg_n13
       -rho*Y[ic13]*Y[ip]*lambda_c13_pg_n14
       -rho*Y[in13]*Y[ip]*lambda_n13_pg_o14
       -rho*Y[in14]*Y[ip]*lambda_n14_pg_o15
       -rho*Y[in15]*Y[ip]*lambda_n15_pa_c12
       )

    dYdt[ihe4] = (
       +rho*Y[in15]*Y[ip]*lambda_n15_pa_c12
       )

    dYdt[ic12] = (
       -rho*Y[ic12]*Y[ip]*lambda_c12_pg_n13
       +rho*Y[in15]*Y[ip]*lambda_n15_pa_c12
       )

    dYdt[ic13] = (
       -rho*Y[ic13]*Y[ip]*lambda_c13_pg_n14
       +Y[in13]*lambda_n13_c13
       )

    dYdt[in13] = (
       -Y[in13]*lambda_n13_c13
       -rho*Y[in13]*Y[ip]*lambda_n13_pg_o14
       +rho*Y[ic12]*Y[ip]*lambda_c12_pg_n13
       )

    dYdt[in14] = (
       -rho*Y[in14]*Y[ip]*lambda_n14_pg_o15
       +rho*Y[ic13]*Y[ip]*lambda_c13_pg_n14
       +Y[io14]*lambda_o14_n14
       )

    dYdt[in15] = (
       -rho*Y[in15]*Y[ip]*lambda_n15_pa_c12
       +Y[io15]*lambda_o15_n15
       )

    dYdt[io14] = (
       -Y[io14]*lambda_o14_n14
       +rho*Y[in13]*Y[ip]*lambda_n13_pg_o14
       )

    dYdt[io15] = (
       -Y[io15]*lambda_o15_n15
       +rho*Y[in14]*Y[ip]*lambda_n14_pg_o15
       )

    return dYdt
