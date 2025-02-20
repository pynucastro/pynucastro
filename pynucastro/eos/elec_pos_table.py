import numpy
from pynucastro.eos.fermi_integrals import fermi, dfermi_deta, dfermi_dbeta
from pynucastro.eos.eta_evaluation import (N_ele,
                                           N_pos,
                                           brentq_eta)
from pynucastro.constants import constants

m_e = constants.m_e
hbar = constants.hbar
c = constants.c_light

class ElectronPositronTable:

    def __init__(self, rho_bounds, temp_bounds, ye=1.0, rho_samples_per_decade=20, temp_samples_per_decade=20):
        rho_l, rho_h = rho_bounds
        temp_l, temp_h = temp_bounds

        self.rho_points = self._log_per_decade(rho_l, rho_h, rho_samples_per_decade)
        self.temp_points = self._log_per_decade(temp_l, temp_h, temp_samples_per_decade)

        for rho in self.rho_points:
            for temp in self.temp_points:
                beta = k*temp/(m_e*c**2)
                eta = self._compute_eta(rho, beta)
                _compute_F(rho, eta, beta)
                _compute_dF_dd(rho, eta, beta)
                _compute_dF_dt(rho, eta, beta)
                _compute_d2F_dt2(rho, eta, beta)
                _compute_d2F_dr_dT(rho, eta, beta)
                _compute_d2F_dT2(rho, eta, beta)
                _compute_d3F_dr3(rho, eta, beta)
                _compute_d3F_dr2_dT(rho, eta, beta)
                _compute_d3F_dr_dT2(rho, eta, beta)
                _compute_d3F_dT3(rho, eta, beta)
                _compute_d4F_dr4(rho, eta, beta)
                _compute_d4F_dr3_dT(rho, eta, beta)
                _compute_d4F_dr2_dT2(rho, eta, beta)
                _compute_d4F_dr_dT3(rho, eta, beta)
                _compute_d4F_dT4(rho, eta, beta)

    def _log_per_decade(x_init, x_end, samples_per_decade):
        ndecades = np.log10(x_end) - np.log10(x_init)
        npoints = int(ndecades) * samples_per_decade
        x_points = np.logspace(np.log10(x_init), np.log10(x_end), num=npoints, endpoint=True, base=10)
        return x_points

    def n_ele(eta, beta):
        integrand = nconst * np.sqrt(beta)**3*(fermi(0.5, eta, beta) + beta*fermi(1.5, eta, beta))
        return integrand

    def n_pos(eta, beta):
        eta_eff = -eta-2/beta
        integrand = n_ele(eta_eff, beta)
        return integrand

    def _compute_n_der(eta, beta):
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

    def _compute_eta_beta_der(self, eta, beta):

        dn_ele_deta, dn_pos_deta, dn_ele_dbeta, dn_pos_dbeta = self._compute_n_der(eta, beta)

        deta_dd = 1/(m_u) * 1/(dn_ele_deta - dn_pos_deta)
        deta_dt = - dbeta_dt * (dn_pos_dbeta - dn_ele_dbeta) / (dn_pos_deta - dn_ele_deta)
        dbeta_dd = 0
        dbeta_dt = k/(m_e*c**2)

        return deta_dd, deta_dt, dbeta_dd, dbeta_dt

    def target_function(eta, beta):
        n_ele_matter = rho/m_u
        target = self.n_ele(eta, beta) - self.n_pos(eta,beta) - N_ele_matter
        return target

    def target_fprime(eta, beta):
        dn_ele_deta, dn_pos_deta, _, _ = self._compute_n_der(eta, beta)
        dn_deta = dn_ele_deta - dn_pos_deta
        return dn_deta

    def newton_eta(rho, beta, initial_eta, max_iter):
        root = optimize.newton(func=lambda eta: target_function(rho, eta, beta), x0=initial_eta, fprime=lambda eta: target_fprime(eta, beta), tol=1.0e-15, maxiter=max_iter)
        return root

    def brentq_eta(rho, beta, low_eta, max_eta, max_iter):
        root = optimize.brentq(f=lambda eta: target_function(rho, eta, beta), a=low_eta, b=max_eta, xtol=1.0e-15, maxiter=max_iter)
        return root

    def _compute_eta(rho, beta):
        low_eta = -33
        max_eta = 25
        eta = brentq_eta(rho, beta, low_eta, max_eta, max_iter=150)
        return eta

    def _compute_pressure(self, eta, beta):
        p_ele = 16*np.pi*np.sqrt(2)*m_e**4*c**5*np.sqrt(beta)**5*(fermi(1.5, eta, beta) + \
                       0.5*beta*fermi(2.5, eta, beta))/(3*hbar**3)
        p_pos = 16*np.pi*np.sqrt(2)*m_e**4*c**5*np.sqrt(beta)**5*(fermi(1.5,-eta-2/beta, beta) + \
                       0.5*beta*fermi(2.5,-eta-2/beta, beta))
        p = p_ele + p_pos

        return p_ele, p_pos, p

    def _compute_energy(self, eta, beta):
        eta_eff = -eta - 2/beta
        econst = 8*np.pi*np.sqrt(2)*m_e**4*c**5/h**3
        e_ele = econst * np.sqrt(beta)**5*(fermi(1.5, eta, beta) + \
                beta*fermi(2.5, eta, beta))
        e_pos = econst * np.sqrt(beta)**5*(fermi(1.5, eta_eff, beta) + \
                beta*fermi(2.5, eta_eff, beta)) + 2*m_e*c**2*n_pos(eta, beta)
        return (e_ele + e_pos)/rho

    def _compute_entropy(self, rho, eta, beta):
        V= 1/rho
        temp = m_e*c**2*beta / k
        eta_eff = -eta - 2/beta
        p_ele, p_pos,_ = self._compute_pressure(rho, eta, beta)
        s_ele = (rho*e_ele + p_ele*V - V*n_ele(eta, beta)*eta) / temp
        s_pos = (rho*e_pos + p_pos*V - V*n_pos(eta_eff, beta)*eta_eff) / temp
        return s_ele + s_pos

    def _compute_F(self, rho, eta, beta):
        temp = m_e*c**2*beta / k
        return self._compute_energy(eta,beta) - temp*self._compute_entropy(rho,eta,beta)

    def _compute_dF_dd(self, rho, eta, beta):
        pres = self._compute_pressure(rho,eta, beta)
        return press/rho**2

    def _compute_dF_dt(self, rho, eta, beta):
        entropy = self._compute_entropy(rho, eta, beta)
        return -entropy

    def _compute_d2F_ddd(self, rho, eta, beta):


        deta_dd, _, _, _ = self._compute(eta, beta)
        dp_dd = dp_deta * deta_drho


    def _compute_d2F_dt2(self, rho, eta, beta):
        temp = m_e*c**2*beta / k

        deta_dt =
        dbeta_dt = beta/temp

        econst = 8*np.pi*np.sqrt(2)*m_e**4*c**5/(hbar**3)

        dEle_deta = econst*np.sqrt(beta)**5*(dfermi_deta(1.5, eta, beta) + \
                                beta*dfermi_deta(2.5, eta, beta))
        dEle_dbeta = econst*(5/2)*np.sqrt(beta)**3*(fermi(1.5, eta, beta) + \
                                beta*dfermi(2.5, eta, beta)) + \
                     econst*np.sqrt(beta)**5

        eta_eff = -eta - 2/beta
        dEpos_deta = -econst*np.sqrt(beta)**5*(dfermi_deta(1.5, eta_eff, beta) + \
                                beta*dfermi_deta(2.5, eta_eff, beta))
        dEpos_dbeta = econst*(5/2)*np.sqrt(beta)**3*(fermi(1.5, eta_eff, beta) + \
                                beta*fermi(2.5, eta_eff, beta)) + \
                      econst*np.sqrt(beta)**3*((dfermi_deta(1.5, eta_eff, beta) + \
                                beta*dfermi_deta(2.5, eta_eff, beta))*(2/beta**2) + \
                                    (dfermi_dbeta(1.5, eta_eff, beta) + 0.5*fermi(2.5, eta_eff, beta)+
                                    beta*fermi(2.5, eta_eff, beta))) + \
                                2*mc**2*dn_pos_dbeta(eta,beta)

        dE_dbeta = dEele_dbeta + dEpos_dbeta
        dE_deta = dEele_deta + dEpos_deta
        dE_dt = dE_deta * deta_dt + dE_dbeta * dbeta_dt

        ds_dt = (1/temp) * dE_dt

        return -ds_dt

    def _compute_d2F_dtdd(self, rho, eta, beta):
        temp = m_e*c**2*beta / k
        eta_eff = -eta - 2/beta

        deta_dt = -eta/temp
        dbeta_dt = beta/temp

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


    def _compute_d3F_dr2_dT(self, rho, eta, beta):
        pass

    def _compute_d3F_dT3(self, rho, eta, beta):
        pass

    def _compute_d4F_dr4(self, rho, eta, beta):
        pass

    def _compute_d4F_dr3_dT(self, rho, eta, beta):
        pass

    def _compute_d4F_dr2_dT2(self, rho, eta, beta):
        pass

    def _compute_d4F_dr_dT3(self, rho, eta, beta):
        pass

    def _compute_d4F_dT4(self, rho, eta, beta):
        pass