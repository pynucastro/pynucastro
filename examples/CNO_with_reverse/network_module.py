import numpy as np
from pynucastro.rates import Tfactors
from pynucastro.nucdata import PartitionFunction, PartitionFunctionCollection

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

A = np.zeros((nnuc), dtype=np.int32)

A[ip] = 1
A[ihe4] = 4
A[ic12] = 12
A[ic13] = 13
A[in13] = 13
A[in14] = 14
A[in15] = 15
A[io14] = 14
A[io15] = 15

species_names = ["p", "he4", "c12", "c13", "n13", "n14", "n15", "o14", "o15"]

partfun_collection = PartitionFunctionCollection()

pf_p = partfun_collection.get_partition_function('p')
pf_he4 = partfun_collection.get_partition_function('he4')
pf_c12 = partfun_collection.get_partition_function('c12')
pf_c13 = partfun_collection.get_partition_function('c13')
pf_n13 = partfun_collection.get_partition_function('n13')
pf_n14 = partfun_collection.get_partition_function('n14')
pf_n15 = partfun_collection.get_partition_function('n15')
pf_o14 = partfun_collection.get_partition_function('o14')
pf_o15 = partfun_collection.get_partition_function('o15')

def n13__c13(tf):
    # n13 --> c13
    rate = 0.0
    
    # wc12w
    rate += np.exp(  -6.7601)
    
    return rate

def o14__n14(tf):
    # o14 --> n14
    rate = 0.0
    
    # wc12w
    rate += np.exp(  -4.62354)
    
    return rate

def o15__n15(tf):
    # o15 --> n15
    rate = 0.0
    
    # wc12w
    rate += np.exp(  -5.17053)
    
    return rate

def n13__p_c12(tf):
    # n13 --> p + c12
    rate = 0.0
    
    # ls09n
    rate += np.exp(  40.0408 + -22.5475*tf.T9i + -13.692*tf.T913i + -0.230881*tf.T913
                  + 4.44362*tf.T9 + -3.15898*tf.T953 + 0.833333*tf.lnT9)
    # ls09r
    rate += np.exp(  40.4354 + -26.326*tf.T9i + -5.10735*tf.T913i + -2.24111*tf.T913
                  + 0.148883*tf.T9)
    
    rate = rate * pf_p(tf.T9*1.0e9)
    rate = rate * pf_c12(tf.T9*1.0e9)
    rate = rate / pf_n13(tf.T9*1.0e9)
    
    return rate

def n14__p_c13(tf):
    # n14 --> p + c13
    rate = 0.0
    
    # nacrr
    rate += np.exp(  37.1528 + -93.4071*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953)
    # nacrr
    rate += np.exp(  38.3716 + -101.18*tf.T9i)
    # nacrn
    rate += np.exp(  41.7046 + -87.6256*tf.T9i + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + 0.833333*tf.lnT9)
    
    rate = rate * pf_p(tf.T9*1.0e9)
    rate = rate * pf_c13(tf.T9*1.0e9)
    rate = rate / pf_n14(tf.T9*1.0e9)
    
    return rate

def o14__p_n13(tf):
    # o14 --> p + n13
    rate = 0.0
    
    # lg06n
    rate += np.exp(  42.4234 + -53.7053*tf.T9i + -15.1676*tf.T913i + 0.0955166*tf.T913
                  + 3.0659*tf.T9 + -0.507339*tf.T953 + 0.833333*tf.lnT9)
    # lg06r
    rate += np.exp(  35.2849 + -59.8313*tf.T9i + 1.57122*tf.T913i)
    
    rate = rate * pf_p(tf.T9*1.0e9)
    rate = rate * pf_n13(tf.T9*1.0e9)
    rate = rate / pf_o14(tf.T9*1.0e9)
    
    return rate

def o15__p_n14(tf):
    # o15 --> p + n14
    rate = 0.0
    
    # im05n
    rate += np.exp(  44.1246 + -84.6757*tf.T9i + -15.193*tf.T913i + -4.63975*tf.T913
                  + 9.73458*tf.T9 + -9.55051*tf.T953 + 1.83333*tf.lnT9)
    # im05n
    rate += np.exp(  41.0177 + -84.6757*tf.T9i + -15.193*tf.T913i + -0.161954*tf.T913
                  + -7.52123*tf.T9 + -0.987565*tf.T953 + 0.833333*tf.lnT9)
    # im05r
    rate += np.exp(  30.7435 + -89.5667*tf.T9i
                  + 1.5682*tf.lnT9)
    # im05r
    rate += np.exp(  31.6622 + -87.6737*tf.T9i)
    
    rate = rate * pf_p(tf.T9*1.0e9)
    rate = rate * pf_n14(tf.T9*1.0e9)
    rate = rate / pf_o15(tf.T9*1.0e9)
    
    return rate

def c12__he4_he4_he4(tf):
    # c12 --> he4 + he4 + he4
    rate = 0.0
    
    # fy05n
    rate += np.exp(  45.7734 + -84.4227*tf.T9i + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + 1.66667*tf.lnT9)
    # fy05r
    rate += np.exp(  22.394 + -88.5493*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -10.1653*tf.lnT9)
    # fy05r
    rate += np.exp(  34.9561 + -85.4472*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + 0.83333*tf.lnT9)
    
    rate = rate * pf_he4(tf.T9*1.0e9)
    rate = rate * pf_he4(tf.T9*1.0e9)
    rate = rate * pf_he4(tf.T9*1.0e9)
    rate = rate / pf_c12(tf.T9*1.0e9)
    
    return rate

def p_c12__n13(tf):
    # c12 + p --> n13
    rate = 0.0
    
    # ls09r
    rate += np.exp(  17.5428 + -3.77849*tf.T9i + -5.10735*tf.T913i + -2.24111*tf.T913
                  + 0.148883*tf.T9 + -1.5*tf.lnT9)
    # ls09n
    rate += np.exp(  17.1482 + -13.692*tf.T913i + -0.230881*tf.T913
                  + 4.44362*tf.T9 + -3.15898*tf.T953 + -0.666667*tf.lnT9)
    
    return rate

def p_c13__n14(tf):
    # c13 + p --> n14
    rate = 0.0
    
    # nacrr
    rate += np.exp(  13.9637 + -5.78147*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  15.1825 + -13.5543*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  18.5155 + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + -0.666667*tf.lnT9)
    
    return rate

def p_n13__o14(tf):
    # n13 + p --> o14
    rate = 0.0
    
    # lg06r
    rate += np.exp(  10.9971 + -6.12602*tf.T9i + 1.57122*tf.T913i
                  + -1.5*tf.lnT9)
    # lg06n
    rate += np.exp(  18.1356 + -15.1676*tf.T913i + 0.0955166*tf.T913
                  + 3.0659*tf.T9 + -0.507339*tf.T953 + -0.666667*tf.lnT9)
    
    return rate

def p_n14__o15(tf):
    # n14 + p --> o15
    rate = 0.0
    
    # im05n
    rate += np.exp(  20.1169 + -15.193*tf.T913i + -4.63975*tf.T913
                  + 9.73458*tf.T9 + -9.55051*tf.T953 + 0.333333*tf.lnT9)
    # im05n
    rate += np.exp(  17.01 + -15.193*tf.T913i + -0.161954*tf.T913
                  + -7.52123*tf.T9 + -0.987565*tf.T953 + -0.666667*tf.lnT9)
    # im05r
    rate += np.exp(  6.73578 + -4.891*tf.T9i
                  + 0.0682*tf.lnT9)
    # im05r
    rate += np.exp(  7.65444 + -2.998*tf.T9i
                  + -1.5*tf.lnT9)
    
    return rate

def he4_c12__p_n15(tf):
    # c12 + he4 --> p + n15
    rate = 0.0
    
    # nacrr
    rate += np.exp(  -5.2319 + -59.6491*tf.T9i + 30.8497*tf.T913
                  + -8.50433*tf.T9 + -1.54426*tf.T953 + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  27.118 + -57.6279*tf.T9i + -15.253*tf.T913i + 1.59318*tf.T913
                  + 2.4479*tf.T9 + -2.19708*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  -6.93365 + -58.7917*tf.T9i + 22.7105*tf.T913
                  + -2.90707*tf.T9 + 0.205754*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  20.5388 + -65.034*tf.T9i
                  + -1.5*tf.lnT9)
    
    rate = rate * pf_p(tf.T9*1.0e9)
    rate = rate * pf_n15(tf.T9*1.0e9)
    rate = rate / pf_he4(tf.T9*1.0e9)
    rate = rate / pf_c12(tf.T9*1.0e9)
    
    return rate

def p_n15__he4_c12(tf):
    # n15 + p --> he4 + c12
    rate = 0.0
    
    # nacrr
    rate += np.exp(  -4.87347 + -2.02117*tf.T9i + 30.8497*tf.T913
                  + -8.50433*tf.T9 + -1.54426*tf.T953 + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  27.4764 + -15.253*tf.T913i + 1.59318*tf.T913
                  + 2.4479*tf.T9 + -2.19708*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  -6.57522 + -1.1638*tf.T9i + 22.7105*tf.T913
                  + -2.90707*tf.T9 + 0.205754*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  20.8972 + -7.406*tf.T9i
                  + -1.5*tf.lnT9)
    
    return rate

def he4_he4_he4__c12(tf):
    # he4 + he4 + he4 --> c12
    rate = 0.0
    
    # fy05r
    rate += np.exp(  -24.3505 + -4.12656*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -13.1653*tf.lnT9)
    # fy05r
    rate += np.exp(  -11.7884 + -1.02446*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + -2.16667*tf.lnT9)
    # fy05n
    rate += np.exp(  -0.971052 + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + -1.33333*tf.lnT9)
    
    return rate

def rhs(t, Y, rho, T):

    tf = Tfactors(T)

    lambda_n13__c13 = n13__c13(tf)
    lambda_o14__n14 = o14__n14(tf)
    lambda_o15__n15 = o15__n15(tf)
    lambda_n13__p_c12 = n13__p_c12(tf)
    lambda_n14__p_c13 = n14__p_c13(tf)
    lambda_o14__p_n13 = o14__p_n13(tf)
    lambda_o15__p_n14 = o15__p_n14(tf)
    lambda_c12__he4_he4_he4 = c12__he4_he4_he4(tf)
    lambda_p_c12__n13 = p_c12__n13(tf)
    lambda_p_c13__n14 = p_c13__n14(tf)
    lambda_p_n13__o14 = p_n13__o14(tf)
    lambda_p_n14__o15 = p_n14__o15(tf)
    lambda_he4_c12__p_n15 = he4_c12__p_n15(tf)
    lambda_p_n15__he4_c12 = p_n15__he4_c12(tf)
    lambda_he4_he4_he4__c12 = he4_he4_he4__c12(tf)

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[ip] = (
       -rho*Y[ip]*Y[ic12]*lambda_p_c12__n13
       -rho*Y[ip]*Y[ic13]*lambda_p_c13__n14
       -rho*Y[ip]*Y[in13]*lambda_p_n13__o14
       -rho*Y[ip]*Y[in14]*lambda_p_n14__o15
       -rho*Y[ip]*Y[in15]*lambda_p_n15__he4_c12
       +Y[in13]*lambda_n13__p_c12
       +Y[in14]*lambda_n14__p_c13
       +Y[io14]*lambda_o14__p_n13
       +Y[io15]*lambda_o15__p_n14
       +rho*Y[ihe4]*Y[ic12]*lambda_he4_c12__p_n15
       )

    dYdt[ihe4] = (
       -rho*Y[ihe4]*Y[ic12]*lambda_he4_c12__p_n15
       -3*1.66666666666667e-01*rho**2*Y[ihe4]**3*lambda_he4_he4_he4__c12
       +3*Y[ic12]*lambda_c12__he4_he4_he4
       +rho*Y[ip]*Y[in15]*lambda_p_n15__he4_c12
       )

    dYdt[ic12] = (
       -Y[ic12]*lambda_c12__he4_he4_he4
       -rho*Y[ip]*Y[ic12]*lambda_p_c12__n13
       -rho*Y[ihe4]*Y[ic12]*lambda_he4_c12__p_n15
       +Y[in13]*lambda_n13__p_c12
       +rho*Y[ip]*Y[in15]*lambda_p_n15__he4_c12
       +1.66666666666667e-01*rho**2*Y[ihe4]**3*lambda_he4_he4_he4__c12
       )

    dYdt[ic13] = (
       -rho*Y[ip]*Y[ic13]*lambda_p_c13__n14
       +Y[in13]*lambda_n13__c13
       +Y[in14]*lambda_n14__p_c13
       )

    dYdt[in13] = (
       -Y[in13]*lambda_n13__c13
       -Y[in13]*lambda_n13__p_c12
       -rho*Y[ip]*Y[in13]*lambda_p_n13__o14
       +Y[io14]*lambda_o14__p_n13
       +rho*Y[ip]*Y[ic12]*lambda_p_c12__n13
       )

    dYdt[in14] = (
       -Y[in14]*lambda_n14__p_c13
       -rho*Y[ip]*Y[in14]*lambda_p_n14__o15
       +Y[io14]*lambda_o14__n14
       +Y[io15]*lambda_o15__p_n14
       +rho*Y[ip]*Y[ic13]*lambda_p_c13__n14
       )

    dYdt[in15] = (
       -rho*Y[ip]*Y[in15]*lambda_p_n15__he4_c12
       +Y[io15]*lambda_o15__n15
       +rho*Y[ihe4]*Y[ic12]*lambda_he4_c12__p_n15
       )

    dYdt[io14] = (
       -Y[io14]*lambda_o14__n14
       -Y[io14]*lambda_o14__p_n13
       +rho*Y[ip]*Y[in13]*lambda_p_n13__o14
       )

    dYdt[io15] = (
       -Y[io15]*lambda_o15__n15
       -Y[io15]*lambda_o15__p_n14
       +rho*Y[ip]*Y[in14]*lambda_p_n14__o15
       )

    return dYdt
