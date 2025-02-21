import numpy
from pynucastro.eos.fermi_integrals import fermi, dfermi_deta, dfermi_dbeta
from pynucastro.eos.eta_evaluation import (N_ele,
                                           N_pos,
                                           brentq_eta)
from pynucastro.constants import constants

m_e = constants.m_e
hbar = constants.hbar
c = constants.c_light

class ElectronPositron:

    def __init__(self, rho, temp, ye=1.0):

        self.rho = rho
        self.temp = temp
        self.beta = k*temp  (m_e*c**2)
        self.eta = self._compute_eta(rho, beta)
        self.eta_eff = -eta - 2/beta

        self.f = None
        self.df_dd = None
        self.df_dt = None
        self.d2f_dd2 = None
        self.d2f_dt2 = None
        self.d2f_ddt = None
        self.d3_dd2t = None
        self.d3_ddt2 = None
        self.d4_dd2t2 = None

        self.f = _compute_f(rho, eta, beta)
        self.df_dd = _compute_df_dd(rho, eta, beta)
        self.df_dt = _compute_df_dt(rho, eta, beta)
        self.d2f_dd2 = _compute_df_dd2(rho, eta, beta)
        self.d2f_dt2 = _compute_d2f_dt2(rho, eta, beta)
        self.d2f_ddt = _compute_d2f_ddt(rho, eta, beta)
        self.d3f_dd2t = _compute_d3f_dd2t(rho, eta, beta)
        self.d3_ddt2 = _compute_d3f_ddt2(rho, eta, beta)
        self.d4_dd2t2 = _compute_d3f_dd2t2(rho, eta, beta)

    def n_ele(self):
        eta = self.eta
        beta = self.beta
        quad = nconst * np.sqrt(beta)**3*(fermi(0.5, eta, beta) + beta*fermi(1.5, eta, beta))
        return quad

    def n_pos(self):
        eta = self.eta
        beta=self.beta
        eta_eff = -eta-2/beta
        integrand = n_ele(eta_eff, beta)
        return quad

    def _compute_n_der(self):
        eta = self.eta
        beta = self.beta
        eta_eff = -eta - 2/beta

        nconst = 8*np.pi*np.sqrt(2)*m_e**3*c**3 / h**3

        dn_ele_deta = nconst * np.sqrt(beta)**3 * (dfermi_deta(0.5, eta, beta) + beta*dfermi_deta(1.5, eta, beta))
        dn_pos_deta = -nconst * np.sqrt(beta)**3 * (dfermi_deta(0.5, eta_eff, beta) + beta*dfermi_deta(1.5, eta_eff, beta))
        dn_deta = dn_ele_deta - dn_pos_deta

        dn_ele_dbeta = nconst*(3/2)*np.sqrt(beta)*(fermi(0.5,eta,beta) + beta*fermi(1.5,eta,beta)) + \
                       nconst + np.sqrt(beta)**3 * (dfermi_dbeta(0.5, eta, beta) + fermi(1.5, eta, beta) + \
                       beta*dfermi(1.5, eta, beta))

        dn_pos_dbeta = nconst*(3/2)*np.sqrt(beta)*(fermi(0.5,eta_eff,beta) + beta*fermi(1.5,eta_eff,beta)) + \
                       nconst*np.sqrt(beta)**3*( (dfermi_deta(0.5, eta_eff, beta) + beta*dfermi_deta(1.5, eta_eff, beta))*(2/(beta**2)) + \
                       dfermi_dbeta(0.5, eta_eff, beta) + fermi(1.5, eta_eff, beta) + beta*dfermi_dbeta(1.5, eta_eff, beta))

        return dn_ele_deta, dn_pos_deta, dn_ele_dbeta, dn_pos_dbeta

    def _compute_eta_beta_der(self):
        eta = self.eta
        beta = self.beta

        dn_ele_deta, dn_pos_deta, dn_ele_dbeta, dn_pos_dbeta = self._compute_n_der(eta, beta)
        deta_dd = 1/(m_u) * 1/(dn_ele_deta - dn_pos_deta)
        deta_dt = -dbeta_dt * (dn_pos_dbeta - dn_ele_dbeta) / (dn_pos_deta - dn_ele_deta)
        dbeta_dd = 0
        dbeta_dt = k/(m_e*c**2)

        return deta_dd, deta_dt, dbeta_dd, dbeta_dt

    def target_function(self,eta):
        rho = self.rho
        beta = self.beta
        n_ele_matter = rho/m_u
        target = self.n_ele(eta, beta) - self.n_pos(eta,beta) - n_ele_matter
        return target

    def target_fprime(self, eta):
        beta = self.beta
        dn_ele_deta, dn_pos_deta, _, _ = self._compute_n_der(eta, beta)
        dn_deta = dn_ele_deta - dn_pos_deta
        return dn_deta

    def initial_eta_function(self):
        #should be improved
        return -10

    def newton_eta(self, initial_eta, max_iter):
        initial_eta = self.initial_eta_function()
        root = optimize.newton(func=lambda eta: self.target_function(eta), x0=initial_eta,
                               fprime=lambda eta: self.target_fprime(eta), tol=1.0e-15, maxiter=max_iter)
        return root

    def brentq_eta(self, low_eta, max_eta, max_iter):
        root = optimize.brentq(f=lambda eta: self.target_function(eta), a=low_eta, b=max_eta, xtol=1.0e-15, maxiter=max_iter)
        return root

    def _compute_eta(self):
        low_eta = -33
        max_eta = 25
        eta = self.brentq_eta(low_eta, max_eta, max_iter=150)
        return eta

    def _compute_pressure(self, details=False):
        eta = self.eta
        beta = self.beta

        pconst = 16*np.pi*np.sqrt(2)*m_e**4*c**5/(h**3)
        p_ele = pconst * np.sqrt(beta)**5*(fermi(1.5, eta, beta) + \
                       0.5*beta*fermi(2.5, eta, beta))/(3*hbar**3)
        p_pos = pconst * np.sqrt(beta)**5*(fermi(1.5,-eta-2/beta, beta) + \
                       0.5*beta*fermi(2.5,-eta-2/beta, beta))

        p = p_ele + p_pos

        if details:
            return p_ele, p_pos, p
        else:
            return p

    def _compute_energy(self):
        eta = self.eta
        beta = self.beta

        eta_eff = -eta - 2/beta
        econst = 8*np.pi*np.sqrt(2)*m_e**4*c**5/h**3
        e_ele = econst * np.sqrt(beta)**5*(fermi(1.5, eta, beta) + \
                beta*fermi(2.5, eta, beta))
        e_pos = econst * np.sqrt(beta)**5*(fermi(1.5, eta_eff, beta) + \
                beta*fermi(2.5, eta_eff, beta)) + 2*m_e*c**2*n_pos(eta, beta)
        return (e_ele + e_pos)/rho

    def _compute_entropy(self):
        rho = self.rho
        eta = self.eta
        beta = self.beta

        V= 1/rho
        temp = m_e*c**2*beta / k
        eta_eff = -eta - 2/beta
        p_ele, p_pos, _ = self._compute_pressure(details=True)
        s_ele = (rho*e_ele + p_ele*V - V*n_ele(eta, beta)*eta) / temp
        s_pos = (rho*e_pos + p_pos*V - V*n_pos(eta_eff, beta)*eta_eff) / temp
        return s_ele + s_pos

    def _compute_f(self):
        temp = self.temp
        f = self._compute_energy() - temp*self._compute_entropy()
        return f

    def _compute_df_dd(self):
        rho = self.rho
        pres = self._compute_pressure()
        return press/rho**2

    def _compute_df_dt(self, rho, eta, beta):
        entropy = self._compute_entropy()
        return -entropy

    def _compute_d2f_ddd(self):
        rho = self.rho
        eta = self.eta
        beta = self.beta
        eta_eff = self.eta_eff
        pconst = 16*np.pi*np.sqrt(2)*m_e**4*c**5/(h**3)
        deta_dd, _, _, _ = self._compute_eta_beta_der()

        dp_ele_deta = pconst*np.sqrt(beta)**5 + (dfermi_deta(1.5, eta, beta) + dfermi_deta(2.5, eta, beta))
        dp_pos_deta = - pconst*np.sqrt(beta)**5 + (dfermi_deta(1.5, eta_eff, beta) + dfermi_deta(2.5, eta_eff, beta))
        dp_deta = dp_ele_deta + dp_pos_deta
        dp_dd = dp_deta * deta_drho

        return dp_dd/(rho**2) - 2*p/(rho**3)

    def _compute_d2f_dt2(self, rho, eta, beta):
        temp = self.temp
        eta = self.eta
        beta = self.beta
        eta_eff = self.eta_eff
        _, deta_dt, _, dbeta_dt = self._compute_eta_beta_der()
        econst = 8*np.pi*np.sqrt(2)*m_e**4*c**5/(h**3)

        de_ele_deta = econst*np.sqrt(beta)**5*(dfermi_deta(1.5, eta, beta) + \
                                beta*dfermi_deta(2.5, eta, beta))
        de_ele_dbeta = econst*(5/2)*np.sqrt(beta)**3*(fermi(1.5, eta, beta) + \
                                beta*dfermi(2.5, eta, beta)) + \
                     econst*np.sqrt(beta)**5


        de_pos_deta = -econst*np.sqrt(beta)**5*(dfermi_deta(1.5, eta_eff, beta) + \
                                beta*dfermi_deta(2.5, eta_eff, beta))
        de_pos_dbeta = econst*(5/2)*np.sqrt(beta)**3*(fermi(1.5, eta_eff, beta) + \
                                beta*fermi(2.5, eta_eff, beta)) + \
                      econst*np.sqrt(beta)**3*((dfermi_deta(1.5, eta_eff, beta) + \
                                beta*dfermi_deta(2.5, eta_eff, beta))*(2/beta**2) + \
                                    (dfermi_dbeta(1.5, eta_eff, beta) + 0.5*fermi(2.5, eta_eff, beta)+
                                    beta*fermi(2.5, eta_eff, beta))) + \
                                    2*mc**2*dn_pos_deta(eta,beta)*(2/beta**2) + 2*mc**2*dn_pos_dbeta(eta,beta)

        de_dbeta = de_ele_dbeta + de_pos_dbeta
        de_deta = de_ele_deta + de_pos_deta
        de_dt = dE_deta * deta_dt + dE_dbeta * dbeta_dt

        de_dt = (1/temp) * de_dt

        return -ds_dt

    def _compute_d2f_dtdd(self):
        rho = self.rho
        eta = self.eta
        beta = self.beta
        eta_eff = self.eta_eff
        _, deta_dt, _, dbeta_dt = self._compute_eta_beta_der()
        pconst = 16*np.pi*np.sqrt(2)*m_e**4*c**5/(3*h**3)

        dp_ele_deta = pconst*np.sqrt(beta)**5*(dfermi_deta(1.5, eta, beta) + 0.5*beta*dfermi_deta(2.5, eta, beta))
        dp_ele_dbeta = pconst*(5/2)*np.sqrt(beta)**3*(fermi(1.5, eta, beta) + 0.5*beta*fermi(2.5, eta, beta)) + \
                       pconst*np.sqrt(beta)**5*(dfermi_dbeta(1.5, beta, eta) + 0.5*fermi(2.5,eta,beta) + \
                                                0.5*beta*dfermi_dbeta(2.5, eta, beta))

        dp_pos_deta = -pconst*np.sqrt(beta)**5*(dfermi_deta(1.5, eta, beta) + 0.5*beta*dfermi_deta(2.5, eta, beta))
        dp_pos_dbeta = pconst*(5/2)*np.sqrt(beta)**3*(fermi(1.5, eta_eff, beta) + 0.5*beta*fermi(2.5, eta_eff, beta)) + \
                       pconst*np.sqrt(beta)**5*((dfermi_deta(1.5, eta_eff, beta) + 0.5*beta*dfermi_deta(2.5, eta_eff, beta))* \
                                                (2/beta) + (dfermi_dbeta(2.5,eta_eff,beta) + 0.5*fermi(2.5, eta_eff, beta) + \
                                                 0.5*beta*dfermi_dbeta(2.5, eta_eff, beta)) )

        dp_deta = dp_ele_deta + dp_pos_dbeta
        dp_dt = dp_deta*deta_dt + dp_dbeta*dbeta_dt

        return (1/rho) * dp_dt

    def _compute_d3f_dd2_dt(self, rho, eta, beta):
        pass

    def _compute_d3f_dd_dt2(self, rho, eta, beta):
        pass

    def _compute_d4f_dd4(self, rho, eta, beta):
        pass

class ElectronPositronTable:

    def _log_per_decade(x_init, x_end, samples_per_decade):
        ndecades = np.log10(x_end) - np.log10(x_init)
        npoints = int(ndecades) * samples_per_decade
        x_points = np.logspace(np.log10(x_init), np.log10(x_end), num=npoints, endpoint=True, base=10)
        return x_points
