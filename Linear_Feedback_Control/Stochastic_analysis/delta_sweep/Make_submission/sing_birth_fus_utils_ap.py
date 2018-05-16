
import numpy as np
import pdb
import re
from scipy.integrate import ode 

# Uses alternative parametrization of the linear feedback control



#########################################################
# Functions to characterise deterministic system
#########################################################


def get_ss_main(ms, xi, gamma, beta, kappa, b, mu, delta):
	ws =  (b*beta**2 + b**2*beta*kappa + b*beta*gamma*kappa + b**2*gamma*kappa**2 + 5*b*beta*mu - beta*gamma*mu + b**2*kappa*mu + \
         3*b*gamma*kappa*mu + 2*b**2*delta*ms*mu - 2*b*gamma*ms*mu - 2*b*delta*gamma*ms*mu + 2*gamma**2*ms*mu + 2*b*mu**2 + \
         2*gamma*mu**2 - b*beta*mu*xi + beta*gamma*mu*xi + b**2*kappa*mu*xi - b*gamma*kappa*mu*xi - 2*b**2*delta*ms*mu*xi + \
         2*b*gamma*ms*mu*xi + 2*b*delta*gamma*ms*mu*xi - 2*gamma**2*ms*mu*xi + 2*b*mu**2*xi - 2*gamma*mu**2*xi - \
         (beta + b*kappa + 2*mu)*np.sqrt(gamma**2*mu**2*(-1 + xi)**2 - \
            2*b*gamma*mu*(-1 + xi)*(-beta + gamma*kappa - 2*(-1 + delta)*gamma*ms + mu + mu*xi) + \
            b**2*(beta**2 - 4*(-1 + delta)*gamma*ms*mu*(-1 + xi) + 2*beta*(gamma*kappa + 3*mu - mu*xi) + (gamma*kappa + mu + mu*xi)**2)))/\
       (2.*(b - gamma)**2*mu*(-1 + xi))

	x1 = np.sqrt(gamma**2*mu**2*(-1 + xi)**2 - 2*b*gamma*mu*(-1 + xi)*(-beta + gamma*kappa - 2*(-1 + delta)*gamma*ms + mu + mu*xi) + \
        b**2*(beta**2 - 4*(-1 + delta)*gamma*ms*mu*(-1 + xi) + 2*beta*(gamma*kappa + 3*mu - mu*xi) + (gamma*kappa + mu + mu*xi)**2))
	
	wf = ((beta*x1)/(b*gamma*(beta + (-1 + delta)*gamma*ms)) - (ms*x1)/(b*(beta + (-1 + delta)*gamma*ms)) - \
         (gamma*(beta + gamma*kappa + 2*mu)**2)/((b - gamma)**2*mu*(-1 + xi)) + (2*x1)/((b - gamma)**2*(-1 + xi)) + \
         x1/((b - gamma)*gamma*(-1 + xi)) + (beta*x1)/((b - gamma)**2*mu*(-1 + xi)) + (kappa*x1)/((b - gamma)*mu*(-1 + xi)) + \
         (gamma*kappa*x1)/((b - gamma)**2*mu*(-1 + xi)) + \
         (-beta**2 - 2*(gamma*kappa + mu)*(gamma*kappa + 2*mu) + beta*(-3*gamma*kappa + 2*mu*(-3 + xi)) + \
            2*(-1 + delta)*gamma*ms*mu*(-1 + xi))/((b - gamma)*mu*(-1 + xi)) - \
         ((beta - gamma*ms)*mu*(-1 + xi))/(b*(beta + (-1 + delta)*gamma*ms)) - (x1*xi)/((b - gamma)*gamma*(-1 + xi)) - \
         (beta*kappa*(beta + gamma*kappa - mu*(-3 + xi)) - 2*(-1 + delta)*delta*gamma*ms**2*mu*(-1 + xi) + \
            ms*((-1 + delta)*gamma**2*kappa**2 + beta*((-1 + delta)*gamma*kappa - delta*mu*(-1 + xi)) + \
               gamma*kappa*mu*(-3 + 2*delta + xi) + delta*mu**2*(-1 + xi**2)))/((beta + (-1 + delta)*gamma*ms)*mu*(-1 + xi)))/2.
	 
	mf =  -(ms*(gamma*mu*(-1 + xi) + b*(beta - gamma*(kappa + 2*ms - 2*delta*ms) - mu*(1 + xi)) - \
            np.sqrt(gamma**2*mu**2*(-1 + xi)**2 - 2*b*gamma*mu*(-1 + xi)*(-beta + gamma*kappa - 2*(-1 + delta)*gamma*ms + mu + mu*xi) + \
              b**2*(beta**2 - 4*(-1 + delta)*gamma*ms*mu*(-1 + xi) + 2*beta*(gamma*kappa + 3*mu - mu*xi) + (gamma*kappa + mu + mu*xi)**2))\
            ))/(2.*b*(beta + (-1 + delta)*gamma*ms))

	return ws, wf, mf

def get_n(ms, xi, gamma, beta, kappa, b, mu, delta):
	ws, wf, mf = get_ss_main(ms, xi, gamma, beta, kappa, b, mu, delta)	
	n = ms + mf + ws + wf
	return n

def residual_n_manual(n_target, ms, xi, gamma, beta, kappa, b, mu, delta, kappa_guess):
	# expect kappa as a vector
	n_val = get_n(ms, xi, gamma, beta, kappa, b, mu, delta)
	n_val[np.isnan(n_val)] = 1e50 + 1e3*(kappa[np.isnan(n_val)]-kappa_guess)*(kappa[np.isnan(n_val)]-kappa_guess)
  	n_val[~np.isnan(n_val)] = (n_val[~np.isnan(n_val)] - n_target)*(n_val[~np.isnan(n_val)] - n_target)
  	return np.log10(n_val)

def get_ic(h_target, xi, gamma, beta, kappa, b, mu, delta):
	ms_space = np.linspace(0,2000,2001)
	ws_space, wf_space, mf_space = get_ss_main(ms_space, xi, gamma, beta, kappa, b, mu, delta)

	n = ms_space+mf_space+ws_space+wf_space
	n = n.astype(float)
	h_space = (ms_space+mf_space)/n

	ind = np.argmin(np.abs(h_space-h_target))

	return ws_space[ind], ms_space[ind], wf_space[ind], mf_space[ind]

