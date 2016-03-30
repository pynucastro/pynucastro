import numpy as np
import reaclib

ihe4 = 0
ic12 = 1
nnuc = 2

A = np.zeros((nnuc), dtype=np.int32)

A[ihe4] = 4
A[ic12] = 12

def c12_gaa_he4(tf):
    # c12 --> he4 + he4 + he4
    rate = 0.0
    
    # fy05rv
    rate += np.exp(  34.9561 + -85.4472*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + 0.83333*tf.lnT9)
    # fy05nv
    rate += np.exp(  45.7734 + -84.4227*tf.T9i + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + 1.66667*tf.lnT9)
    # fy05rv
    rate += np.exp(  22.394 + -88.5493*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -10.1653*tf.lnT9)
    
    return rate

def he4_aag_c12(tf):
    # he4 + he4 + he4 --> c12
    rate = 0.0
    
    # fy05n
    rate += np.exp(  -0.971052 + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + -1.33333*tf.lnT9)
    # fy05r
    rate += np.exp(  -24.3505 + -4.12656*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -13.1653*tf.lnT9)
    # fy05r
    rate += np.exp(  -11.7884 + -1.02446*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + -2.16667*tf.lnT9)
    
    return rate

def rhs(t, Y, rho, T):

    tf = reaclib.Tfactors(T)

    lambda_c12_gaa_he4 = c12_gaa_he4(tf)
    lambda_he4_aag_c12 = he4_aag_c12(tf)

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[ihe4] = (
       -3*0.166666666666667*rho**2*Y[ihe4]**3*lambda_he4_aag_c12
       +3*Y[ic12]*lambda_c12_gaa_he4
       )

    dYdt[ic12] = (
       -Y[ic12]*lambda_c12_gaa_he4
       +0.166666666666667*rho**2*Y[ihe4]**3*lambda_he4_aag_c12
       )

    return dYdt
