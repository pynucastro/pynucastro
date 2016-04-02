import numpy as np
import reaclib

in = 0
ip = 1
id = 2
it = 3
ihe3 = 4
ihe4 = 5
ihe6 = 6
ili7 = 7
ib17 = 8
ic12 = 9
ic14 = 10
nnuc = 11

A = np.zeros((nnuc), dtype=np.int32)

A[in] = 1
A[ip] = 1
A[id] = 2
A[it] = 3
A[ihe3] = 3
A[ihe4] = 4
A[ihe6] = 6
A[ili7] = 7
A[ib17] = 17
A[ic12] = 12
A[ic14] = 14

def b17_nnn_c14(tf):
    # b17 --> n + n + n + c14
    rate = 0.0
    
    # wc12w
    rate += np.exp(  1.56352)
    
    return rate

def he3_he3pp_he4(tf):
    # he3 + he3 --> p + p + he4
    rate = 0.0
    
    # nacrn
    rate += np.exp(  24.7788 + -12.277*tf.T913i + -0.103699*tf.T913
                  + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -0.666667*tf.lnT9)
    
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

def he4_npahe3_li7(tf):
    # n + p + he4 + he4 --> he3 + li7
    rate = 0.0
    
    # mafonv
    rate += np.exp(  -14.8862 + -111.725*tf.T9i + -17.989*tf.T913i + -1.57523e-09*tf.T913
                  + 1.45934e-10*tf.T9 + -1.15341e-11*tf.T953 + -3.66667*tf.lnT9)
    
    return rate

def he4_pphe3_he3(tf):
    # p + p + he4 --> he3 + he3
    rate = 0.0
    
    # nacrnv
    rate += np.exp(  2.98257 + -149.222*tf.T9i + -12.277*tf.T913i + -0.103699*tf.T913
                  + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -2.16667*tf.lnT9)
    
    return rate

def he6_gnn_he4(tf):
    # he6 --> n + n + he4
    rate = 0.0
    
    # cf88rv
    rate += np.exp(  22.178 + -20.8994*tf.T9i + 0.694279*tf.T913i + -3.33326*tf.T913
                  + 0.507932*tf.T9 + -0.0427342*tf.T953 + 2.0*tf.lnT9)
    
    return rate

def li7_tnna_he4(tf):
    # t + li7 --> n + n + he4 + he4
    rate = 0.0
    
    # mafon
    rate += np.exp(  27.5043 + -5.31692e-12*tf.T9i + -11.333*tf.T913i + -2.24192e-09*tf.T913
                  + 2.21773e-10*tf.T9 + -1.83941e-11*tf.T953 + -0.666667*tf.lnT9)
    
    return rate

def n_p(tf):
    # n --> p
    rate = 0.0
    
    # wc12w
    rate += np.exp(  -6.78161)
    
    return rate

def p_ng_d(tf):
    # n + p --> d
    rate = 0.0
    
    # an06n
    rate += np.exp(  10.7548 + -2.30472*tf.T913
                  + -0.887862*tf.T9 + 0.137663*tf.T953)
    # an06n
    rate += np.exp(  8.84688 + -0.0102082*tf.T913
                  + -0.0893959*tf.T9 + 0.00696704*tf.T953 + 1.0*tf.lnT9)
    # an06n
    rate += np.exp(  12.3687 + -2.70618*tf.T913
                  + 0.11718*tf.T9 + -0.00312788*tf.T953 + 0.469127*tf.lnT9)
    
    return rate

def t_gn_d(tf):
    # t --> n + d
    rate = 0.0
    
    # nk06nv
    rate += np.exp(  28.869 + -72.6136*tf.T9i
                  + 1.575*tf.lnT9)
    # nk06nv
    rate += np.exp(  30.1124 + -72.6136*tf.T9i
                  + 2.5*tf.lnT9)
    
    return rate

def t_pn_he3(tf):
    # p + t --> n + he3
    rate = 0.0
    
    # v
    rate += np.exp(  20.3787 + -8.86352*tf.T9i + -0.332788*tf.T913
                  + -0.700485*tf.T9 + 0.0976521*tf.T953)
    # v
    rate += np.exp(  19.2762 + -8.86352*tf.T9i + 0.0438557*tf.T913
                  + -0.201527*tf.T9 + 0.0153433*tf.T953 + 1.0*tf.lnT9)
    
    return rate

def rhs(t, Y, rho, T):

    tf = reaclib.Tfactors(T)

    lambda_b17_nnn_c14 = b17_nnn_c14(tf)
    lambda_he3_he3pp_he4 = he3_he3pp_he4(tf)
    lambda_he4_aag_c12 = he4_aag_c12(tf)
    lambda_he4_npahe3_li7 = he4_npahe3_li7(tf)
    lambda_he4_pphe3_he3 = he4_pphe3_he3(tf)
    lambda_he6_gnn_he4 = he6_gnn_he4(tf)
    lambda_li7_tnna_he4 = li7_tnna_he4(tf)
    lambda_n_p = n_p(tf)
    lambda_p_ng_d = p_ng_d(tf)
    lambda_t_gn_d = t_gn_d(tf)
    lambda_t_pn_he3 = t_pn_he3(tf)

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[in] = (
       -0.5*rho**3*Y[ihe4]**2*Y[ip]*Y[in]*lambda_he4_npahe3_li7
       -Y[in]*lambda_n_p
       -rho*Y[ip]*Y[in]*lambda_p_ng_d
       +3*Y[ib17]*lambda_b17_nnn_c14
       +2*Y[ihe6]*lambda_he6_gnn_he4
       +2*rho*Y[ili7]*Y[it]*lambda_li7_tnna_he4
       +Y[it]*lambda_t_gn_d
       +rho*Y[ip]*Y[it]*lambda_t_pn_he3
       )

    dYdt[ip] = (
       -0.5*rho**3*Y[ihe4]**2*Y[ip]*Y[in]*lambda_he4_npahe3_li7
       -2*0.5*rho**2*Y[ihe4]*Y[ip]**2*lambda_he4_pphe3_he3
       -rho*Y[ip]*Y[in]*lambda_p_ng_d
       -rho*Y[ip]*Y[it]*lambda_t_pn_he3
       +2*0.5*rho*Y[ihe3]**2*lambda_he3_he3pp_he4
       +Y[in]*lambda_n_p
       )

    dYdt[id] = (
       +rho*Y[ip]*Y[in]*lambda_p_ng_d
       +Y[it]*lambda_t_gn_d
       )

    dYdt[it] = (
       -rho*Y[ili7]*Y[it]*lambda_li7_tnna_he4
       -Y[it]*lambda_t_gn_d
       -rho*Y[ip]*Y[it]*lambda_t_pn_he3
       )

    dYdt[ihe3] = (
       -2*0.5*rho*Y[ihe3]**2*lambda_he3_he3pp_he4
       +0.5*rho**3*Y[ihe4]**2*Y[ip]*Y[in]*lambda_he4_npahe3_li7
       +2*0.5*rho**2*Y[ihe4]*Y[ip]**2*lambda_he4_pphe3_he3
       +rho*Y[ip]*Y[it]*lambda_t_pn_he3
       )

    dYdt[ihe4] = (
       -3*0.166666666667*rho**2*Y[ihe4]**3*lambda_he4_aag_c12
       -2*0.5*rho**3*Y[ihe4]**2*Y[ip]*Y[in]*lambda_he4_npahe3_li7
       -0.5*rho**2*Y[ihe4]*Y[ip]**2*lambda_he4_pphe3_he3
       +0.5*rho*Y[ihe3]**2*lambda_he3_he3pp_he4
       +Y[ihe6]*lambda_he6_gnn_he4
       +2*rho*Y[ili7]*Y[it]*lambda_li7_tnna_he4
       )

    dYdt[ihe6] = (
       -Y[ihe6]*lambda_he6_gnn_he4
       )

    dYdt[ili7] = (
       -rho*Y[ili7]*Y[it]*lambda_li7_tnna_he4
       +0.5*rho**3*Y[ihe4]**2*Y[ip]*Y[in]*lambda_he4_npahe3_li7
       )

    dYdt[ib17] = (
       -Y[ib17]*lambda_b17_nnn_c14
       )

    dYdt[ic12] = (
       +0.166666666667*rho**2*Y[ihe4]**3*lambda_he4_aag_c12
       )

    dYdt[ic14] = (
       +Y[ib17]*lambda_b17_nnn_c14
       )

    return dYdt
