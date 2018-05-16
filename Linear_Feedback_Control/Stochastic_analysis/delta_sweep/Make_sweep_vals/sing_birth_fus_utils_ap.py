
import numpy as np
import pdb
import re
from scipy.integrate import ode 

# Uses alternative parametrization of the linear feedback control

###########################################
# Parameters
###########################################

# def get_initial_params():
# 	top_dir = '/home/juvid/Dropbox/Work/Mit_and_Metabolism/Networks_mtDNA_dynamics/'
# 	het = 30
# 	h0 = het/100.0
# 	datafilename = top_dir+'Stochastic_sim/Closed_loop_control/Singleton_birth_sing/1e4_iter/h_{}/output.txt'.format(het)
# 	datafile = open(datafilename, 'r')

# 	ics = []
# 	for i in range(15):
# 		line = datafile.readline()
# 		if i == 1:
# 			match = re.search('^NREPEATS = (.*\w)$', line)
# 			nrepeats = int(match.groups()[0])
# 		elif i == 4:
# 			vals = line.split(',')
# 			for item in vals:
# 				match = re.search('^.*\w = ([0-9]+)\n?$', item)
# 				ics.append(int(match.groups()[0]))
# 		elif i == 5:
# 			match = re.search('^xi = ([0-9]+.[0-9.]+)\n$', line)
# 			xi = float(match.groups()[0])
# 		elif i == 6:
# 			match = re.search('^beta = ([0-9]+.[0-9.]+)\n$', line)
# 			beta = float(match.groups()[0])
# 		elif i == 7:
# 			match = re.search('^gamma = ([0-9]+.[0-9.]+)\n$', line)
# 			gamma = float(match.groups()[0])
# 		elif i == 8:
# 			match = re.search('^a = ([0-9]+.[0-9.]+)\n$', line)
# 			a = float(match.groups()[0])
# 		elif i == 9:
# 			match = re.search('^b = ([0-9]+.[0-9.]+)\n$', line)
# 			b = float(match.groups()[0])
# 		elif i == 10:
# 			match = re.search('^mu = ([0-9]+.[0-9.]+)\n$', line)
# 			mu = float(match.groups()[0])
# 		elif i == 11:
# 			match = re.search('^TFinal = ([0-9]+.[0-9.]+)\n$', line)
# 			TFinal = float(match.groups()[0])

# 	ws_init, ms_init, wf_init, mf_init = ics # unpack initial conditions

# 	return xi, beta, gamma, a, b, mu, TFinal, ws_init, ms_init, wf_init, mf_init

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



def two_dim_network(t, y, beta, gamma, kappa, b, xi, mu, delta):
	ws, wf, ms, mf = y
	wsdot = -1.0 * gamma * ws * ws - gamma * ws * wf + beta * wf - (mu + b * (kappa - (ws + wf + delta * ms +  delta * mf))) * ws - mu * ws - gamma * mf * ws - gamma * ws * ms
	msdot = -1.0 * gamma * ms * ms - gamma * ms * mf + beta * mf - (mu + b * (kappa - (ws + wf + delta * ms +  delta * mf))) * ms - mu * ms - gamma * wf * ms - gamma * ws * ms
	
	wfdot = 1.0 * gamma * ws * ws + gamma * ws * wf - beta * wf + (mu + b * (kappa - (ws + wf + delta * ms +  delta * mf))) * (2.0*ws + wf) - xi * mu * wf + gamma * mf * ws + gamma * ws * ms
	mfdot = 1.0 * gamma * ms * ms + gamma * ms * mf - beta * mf + (mu + b * (kappa - (ws + wf + delta * ms +  delta * mf))) * (2.0*ms + mf) - xi * mu * mf + gamma * wf * ms + gamma * ws * ms

	return [wsdot, wfdot, msdot, mfdot]

def make_trajectory(axs, ws_init, wf_init, ms_init, mf_init, params, isfirst, TMAX):
	beta, gamma, kappa, b, xi, mu, delta = params


	y0, t0 = [ws_init, wf_init, ms_init, mf_init], 0.0
	r = ode(two_dim_network)
	r.set_initial_value(y0, t0).set_f_params(beta, gamma, kappa, b, xi, mu, delta)
	t1 = TMAX
	dt = 1

	ws_sol = [] 
	wf_sol = []
	ms_sol = []
	mf_sol = []


	t_sol = []
	t_sol.append(t0)
	ws_sol.append(y0[0])
	wf_sol.append(y0[1])
	ms_sol.append(y0[2])
	mf_sol.append(y0[3])

	
	while r.successful() and r.t < t1:
		r.integrate(r.t+dt)
		t_sol.append(r.t); ws_sol.append(r.y[0]); wf_sol.append(r.y[1]); ms_sol.append(r.y[2]); mf_sol.append(r.y[3])

	if r.successful() == False:
		print('integration error')
		return

	ws_final = ws_sol[-1] 
	wf_final = wf_sol[-1]
	ms_final = ms_sol[-1]
	mf_final = mf_sol[-1]

	h_final = (mf_final + ms_final) / float(ws_final + wf_final + ms_final + mf_final)

	ws_final_theory, wf_final_theory, mf_final_theory = get_ss_main(ms_final, xi, gamma, beta, kappa, b, mu, delta)

	n_final_theory = ms_final + ws_final_theory +  wf_final_theory + mf_final_theory	

	t_sol = np.array(t_sol); ws_sol = np.array(ws_sol); wf_sol = np.array(wf_sol); ms_sol = np.array(ms_sol); mf_sol = np.array(mf_sol)
	
	if isfirst:
		ax = axs[0]		
		ax.plot(ws_sol + wf_sol, ms_sol + mf_sol, '-b', alpha = 0.2, label = 'Trajectory')
		ax.plot(ws_init + wf_init, ms_init + mf_init, '.b', alpha = 0.2, label = 'Initial condition')
		ax.plot(ws_sol[-1] + wf_sol[-1], ms_sol[-1] + mf_sol[-1], 'sk', alpha = 0.2, label = 'Final condition')

		ax = axs[1]
		ax.plot(ws_sol + ms_sol, wf_sol + mf_sol, '-b', alpha = 0.2, label = 'Trajectory')
		ax.plot(ws_init + ms_init, wf_init + mf_init, '.b', alpha = 0.2, label = 'Initial condition')
		ax.plot(ws_sol[-1] + ms_sol[-1], wf_sol[-1] + mf_sol[-1], 'sk', alpha = 0.2, label = 'Final condition')

		ax = axs[2]
		ax.plot(t_sol, ws_sol + ms_sol + wf_sol + mf_sol, '-b', alpha = 0.2, label = 'Trajectory')
		ax.plot(t_sol[0], ws_init + ms_init + wf_init + mf_init, '.b', alpha = 0.2, label = 'Initial condition')
		ax.plot(t_sol[-1], ws_sol[-1] + ms_sol[-1] + wf_sol[-1] + mf_sol[-1], 'sk', alpha = 0.2, label = 'Final condition')

		ax.plot(t_sol[-1], n_final_theory, 'xr', alpha = 0.5, label = 'Final condition theory')
	else:
		ax = axs[0]		
		ax.plot(ws_sol + wf_sol, ms_sol + mf_sol, '-b', alpha = 0.2)
		ax.plot(ws_sol[-1] + wf_sol[-1], ms_sol[-1] + mf_sol[-1], 'sk', alpha = 0.2)
		ax.plot(ws_init + wf_init, ms_init + mf_init, '.b', alpha = 0.2)

		ax = axs[1]
		ax.plot(ws_sol + ms_sol, wf_sol + mf_sol, '-b', alpha = 0.2)
		ax.plot(ws_sol[-1] + ms_sol[-1], wf_sol[-1] + mf_sol[-1], 'sk', alpha = 0.2)
		ax.plot(ws_init + ms_init, wf_init + mf_init, '.b', alpha = 0.2)

		ax = axs[2]
		ax.plot(t_sol, ws_sol + ms_sol + wf_sol + mf_sol, '-b', alpha = 0.2)
		ax.plot(t_sol[-1], ws_sol[-1] + ms_sol[-1] + wf_sol[-1] + mf_sol[-1], 'sk', alpha = 0.2)
		ax.plot(t_sol[0], ws_init + ms_init + wf_init + mf_init, '.b', alpha = 0.2)

		ax.plot(t_sol[-1], n_final_theory, 'xr', alpha = 0.5)
	return h_final

def jacobian_matrix(xi, gamma, beta, kappa, b, mu, delta, ms, mf, ws, wf):
	J=np.array([
	[-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws,
	b*delta*ws - gamma*ws,
	beta + b*ws - gamma*ws,
	b*delta*ws - gamma*ws],
	[b*ms - gamma*ms,
	-(gamma*mf) + b*delta*ms - 2*gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws,
	b*ms - gamma*ms,
	beta + b*delta*ms - gamma*ms],
	[gamma*mf + gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + 2*gamma*ws - b*(wf + 2*ws),
	gamma*ws - b*delta*(wf + 2*ws),
	-beta + mu + b*(kappa - delta*mf - delta*ms - wf - ws) + gamma*ws - b*(wf + 2*ws) - mu*xi,
	gamma*ws - b*delta*(wf + 2*ws)],
	[gamma*ms - b*(mf + 2*ms),
	gamma*mf + 2*gamma*ms - b*delta*(mf + 2*ms) + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws,
	gamma*ms - b*(mf + 2*ms),
	-beta + gamma*ms - b*delta*(mf + 2*ms) + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi]])

	return J

##########################################
# System size expansion
###########################################
# See /home/juvid/Dropbox/Work/Mit_and_Metabolism/Networks_mtDNA_dynamics/SSE/Mathematica_nbs
def integrate_SSE_network_sing_birth_fus(t, y, xi, gamma, beta, kappa, b, mu, delta):

	ws,wf,ms,mf,cxi1xi1,cxi1xi2,cxi1xi3,cxi1xi4,cxi2xi1,cxi2xi2,cxi2xi3,cxi2xi4,cxi3xi1,cxi3xi2,cxi3xi3,cxi3xi4,cxi4xi1,cxi4xi2,cxi4xi3,cxi4xi4 = y




	dwsdt = beta*wf - gamma*mf*ws - gamma*ms*ws - mu*ws - gamma*wf*ws - (mu + b*(kappa - delta*mf - delta*ms - wf - ws))*ws - gamma*ws**2
	dwfdt = -(beta*wf) + wf*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*mf*ws + gamma*ms*ws + gamma*wf*ws + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws))*ws + gamma*ws**2 - mu*wf*xi
	dmsdt = beta*mf - gamma*mf*ms - gamma*ms**2 - ms*mu - gamma*ms*wf - ms*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - gamma*ms*ws
	dmfdt = -(beta*mf) + gamma*mf*ms + gamma*ms**2 + gamma*ms*wf + mf*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + 2*ms*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ms*ws - mf*mu*xi


	dcxi1xi1dt = beta*wf + 2*cxi1xi1*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + 2*cxi2xi1*(beta + b*ws - gamma*ws) + 2*cxi3xi1*(b*delta*ws - gamma*ws) + 2*cxi4xi1*(b*delta*ws - gamma*ws) + ws*(2*mu + b*(kappa - delta*(mf + ms) - wf - ws) + gamma*(mf + ms + wf + 2*ws))
	dcxi1xi2dt = -(beta*wf) + cxi1xi2*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi2xi2*(beta + b*ws - gamma*ws) + cxi3xi2*(b*delta*ws - gamma*ws) + cxi4xi2*(b*delta*ws - gamma*ws) + cxi3xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi1*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + ws*(-2*mu + 2*b*(-kappa + delta*(mf + ms) + wf + ws) - gamma*(mf + ms + wf + 2*ws)) + cxi2xi1*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi1xi3dt = cxi1xi1*(b*ms - gamma*ms) + cxi2xi1*(b*ms - gamma*ms) + cxi4xi1*(beta + b*delta*ms - gamma*ms) + gamma*ms*ws + cxi1xi3*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi3xi1*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi2xi3*(beta + b*ws - gamma*ws) + cxi3xi3*(b*delta*ws - gamma*ws) + cxi4xi3*(b*delta*ws - gamma*ws)
	dcxi1xi4dt = cxi1xi1*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi1*(-(b*mf) - 2*b*ms + gamma*ms) - gamma*ms*ws + cxi1xi4*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi2xi4*(beta + b*ws - gamma*ws) + cxi3xi4*(b*delta*ws - gamma*ws) + cxi4xi4*(b*delta*ws - gamma*ws) + cxi3xi1*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + cxi4xi1*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)
	dcxi2xi1dt = -(beta*wf) + cxi1xi2*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi2xi2*(beta + b*ws - gamma*ws) + cxi3xi2*(b*delta*ws - gamma*ws) + cxi4xi2*(b*delta*ws - gamma*ws) + cxi3xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi1*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + ws*(-2*mu + 2*b*(-kappa + delta*(mf + ms) + wf + ws) - gamma*(mf + ms + wf + 2*ws)) + cxi2xi1*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi2xi2dt = beta*wf + mu*wf + gamma*mf*ws + gamma*ms*ws + 4*mu*ws + gamma*wf*ws + 2*gamma*ws**2 + b*(kappa - delta*(mf + ms) - wf - ws)*(wf + 4*ws) + 2*cxi3xi2*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + 2*cxi4xi2*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + 2*cxi1xi2*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + mu*wf*xi + 2*cxi2xi2*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi2xi3dt = cxi1xi2*(b*ms - gamma*ms) + cxi2xi2*(b*ms - gamma*ms) + cxi4xi2*(beta + b*delta*ms - gamma*ms) - gamma*ms*ws + cxi3xi2*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi3xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi3*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + cxi2xi3*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi2xi4dt = cxi1xi2*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi2*(-(b*mf) - 2*b*ms + gamma*ms) + gamma*ms*ws + cxi3xi2*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + cxi3xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi4*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + cxi4xi2*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi) + cxi2xi4*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi3xi1dt = cxi1xi1*(b*ms - gamma*ms) + cxi2xi1*(b*ms - gamma*ms) + cxi4xi1*(beta + b*delta*ms - gamma*ms) + gamma*ms*ws + cxi1xi3*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi3xi1*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi2xi3*(beta + b*ws - gamma*ws) + cxi3xi3*(b*delta*ws - gamma*ws) + cxi4xi3*(b*delta*ws - gamma*ws)
	dcxi3xi2dt = cxi1xi2*(b*ms - gamma*ms) + cxi2xi2*(b*ms - gamma*ms) + cxi4xi2*(beta + b*delta*ms - gamma*ms) - gamma*ms*ws + cxi3xi2*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi3xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi3*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + cxi2xi3*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi3xi3dt = beta*mf + 2*cxi1xi3*(b*ms - gamma*ms) + 2*cxi2xi3*(b*ms - gamma*ms) + 2*cxi4xi3*(beta + b*delta*ms - gamma*ms) + 2*cxi3xi3*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + ms*(2*mu + b*(kappa - delta*(mf + ms) - wf - ws) + gamma*(mf + 2*ms + wf + ws))
	dcxi3xi4dt = -(beta*mf) + cxi1xi4*(b*ms - gamma*ms) + cxi2xi4*(b*ms - gamma*ms) + cxi4xi4*(beta + b*delta*ms - gamma*ms) + cxi1xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi3xi4*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi3xi3*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + ms*(-2*mu - gamma*(mf + 2*ms + wf + ws) + 2*b*(-kappa + delta*(mf + ms) + wf + ws)) + cxi4xi3*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)
	dcxi4xi1dt = cxi1xi1*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi1*(-(b*mf) - 2*b*ms + gamma*ms) - gamma*ms*ws + cxi1xi4*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi2xi4*(beta + b*ws - gamma*ws) + cxi3xi4*(b*delta*ws - gamma*ws) + cxi4xi4*(b*delta*ws - gamma*ws) + cxi3xi1*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + cxi4xi1*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)
	dcxi4xi2dt = cxi1xi2*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi2*(-(b*mf) - 2*b*ms + gamma*ms) + gamma*ms*ws + cxi3xi2*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + cxi3xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi4*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + cxi4xi2*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi) + cxi2xi4*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi4xi3dt = -(beta*mf) + cxi1xi4*(b*ms - gamma*ms) + cxi2xi4*(b*ms - gamma*ms) + cxi4xi4*(beta + b*delta*ms - gamma*ms) + cxi1xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi3xi4*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi3xi3*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + ms*(-2*mu - gamma*(mf + 2*ms + wf + ws) + 2*b*(-kappa + delta*(mf + ms) + wf + ws)) + cxi4xi3*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)
	dcxi4xi4dt = beta*mf + gamma*mf*ms + 2*gamma*ms**2 + 2*cxi1xi4*(-(b*mf) - 2*b*ms + gamma*ms) + 2*cxi2xi4*(-(b*mf) - 2*b*ms + gamma*ms) + mf*mu + 4*ms*mu + gamma*ms*wf + gamma*ms*ws - b*(mf + 4*ms)*(-kappa + delta*(mf + ms) + wf + ws) + 2*cxi3xi4*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + mf*mu*xi + 2*cxi4xi4*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)

	deriv = [dwsdt,dwfdt,dmsdt,dmfdt,
			 dcxi1xi1dt,dcxi1xi2dt,dcxi1xi3dt,dcxi1xi4dt,
			 dcxi2xi1dt,dcxi2xi2dt,dcxi2xi3dt,dcxi2xi4dt,
			 dcxi3xi1dt,dcxi3xi2dt,dcxi3xi3dt,dcxi3xi4dt,
			 dcxi4xi1dt,dcxi4xi2dt,dcxi4xi3dt,dcxi4xi4dt]
	return deriv

def partial_deriv_h(ws,wf,ms,mf):
     dhdws = -((mf + ms)/(mf + ms + wf + ws)**2)
     dhdwf = -((mf + ms)/(mf + ms + wf + ws)**2)
     dhdms = -((mf + ms)/(mf + ms + wf + ws)**2) + 1/(mf + ms + wf + ws)
     dhdmf = -((mf + ms)/(mf + ms + wf + ws)**2) + 1/(mf + ms + wf + ws)
     return dhdws,dhdwf,dhdms,dhdmf

def partial_deriv_fs(ws,wf,ms,mf):
     dfsdws = -((ms + ws)/(mf + ms + wf + ws)**2) + 1/(mf + ms + wf + ws)
     dfsdwf = -((ms + ws)/(mf + ms + wf + ws)**2)
     dfsdms = -((ms + ws)/(mf + ms + wf + ws)**2) + 1/(mf + ms + wf + ws)
     dfsdmf = -((ms + ws)/(mf + ms + wf + ws)**2)
     return dfsdws,dfsdwf,dfsdms,dfsdmf
	

#a = mu + b*kappa

# dwsdt = beta*wf - gamma*mf*ws - gamma*ms*ws - mu*ws - gamma*wf*ws - gamma*ws**2 - ws*(a - b*(delta*mf + delta*ms + wf + ws))
# dwfdt = -(beta*wf) + gamma*mf*ws + gamma*ms*ws + gamma*wf*ws + gamma*ws**2 + (2*ws + wf)*(a - b*(delta*mf + delta*ms + wf + ws)) - mu*wf*xi
# dmsdt =  beta*mf - gamma*mf*ms - gamma*ms**2 - ms*mu - gamma*ms*wf - gamma*ms*ws - ms*(a - b*(delta*mf + delta*ms + wf + ws))
# dmfdt = -(beta*mf) + gamma*mf*ms + gamma*ms**2 + gamma*ms*wf + gamma*ms*ws + (mf+ 2*ms)*(a - b*(delta*mf + delta*ms + wf + ws)) - mf*mu*xi

# dcxi1xi1dt =  beta*wf + gamma*mf*ws + gamma*ms*ws + mu*ws + gamma*wf*ws + 2*gamma*ws**2 + 2*cxi2xi1*(beta + b*ws - gamma*ws) + \
#        2*cxi3xi1*(b*delta*ws - gamma*ws) + 2*cxi4xi1*(b*delta*ws - gamma*ws) + \
#        2*cxi1xi1*(-a - gamma*mf - gamma*ms - mu - gamma*wf + b*ws - 2*gamma*ws + b*(delta*mf + delta*ms + wf + ws)) + \
#        ws*(a - b*(delta*(mf + ms) + wf + ws))
# dcxi1xi2dt = -(beta*wf) + cxi2xi2*(beta + b*ws - gamma*ws) + cxi3xi2*(b*delta*ws - gamma*ws) + cxi4xi2*(b*delta*ws - gamma*ws) + \
#        cxi3xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + \
#        cxi1xi2*(-a - gamma*mf - gamma*ms - mu - gamma*wf + b*ws - 2*gamma*ws + b*(delta*mf + delta*ms + wf + ws)) - \
#        ws*(2*a - 2*b*(delta*(mf + ms) + wf + ws) + gamma*(mf + ms + wf + 2*ws)) + \
#        cxi1xi1*(gamma*mf + gamma*ms - b*wf + gamma*wf - 2*b*ws + 2*gamma*ws + 2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi2xi1*(a - beta - b*wf - 2*b*ws + gamma*ws - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi1xi3dt =  cxi1xi1*(b*ms - gamma*ms) + cxi2xi1*(b*ms - gamma*ms) + cxi4xi1*(beta + b*delta*ms - gamma*ms) + gamma*ms*ws + \
#        cxi2xi3*(beta + b*ws - gamma*ws) + cxi3xi3*(b*delta*ws - gamma*ws) + cxi4xi3*(b*delta*ws - gamma*ws) + \
#        cxi1xi3*(-a - gamma*mf - gamma*ms - mu - gamma*wf + b*ws - 2*gamma*ws + b*(delta*mf + delta*ms + wf + ws)) + \
#        cxi3xi1*(-a - gamma*mf + b*delta*ms - 2*gamma*ms - mu - gamma*wf - gamma*ws + b*(delta*mf + delta*ms + wf + ws))
# dcxi1xi4dt = cxi1xi1*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi1*(-(b*mf) - 2*b*ms + gamma*ms) - gamma*ms*ws + cxi2xi4*(beta + b*ws - gamma*ws) + \
#        cxi3xi4*(b*delta*ws - gamma*ws) + cxi4xi4*(b*delta*ws - gamma*ws) + \
#        cxi1xi4*(-a - gamma*mf - gamma*ms - mu - gamma*wf + b*ws - 2*gamma*ws + b*(delta*mf + delta*ms + wf + ws)) + \
#        cxi3xi1*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + gamma*ws + \
#           2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi4xi1*(a - beta - b*delta*mf - 2*b*delta*ms + gamma*ms - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi2xi1dt = -(beta*wf) + cxi2xi2*(beta + b*ws - gamma*ws) + cxi3xi2*(b*delta*ws - gamma*ws) + cxi4xi2*(b*delta*ws - gamma*ws) + \
#        cxi3xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + \
#        cxi1xi2*(-a - gamma*mf - gamma*ms - mu - gamma*wf + b*ws - 2*gamma*ws + b*(delta*mf + delta*ms + wf + ws)) - \
#        ws*(2*a - 2*b*(delta*(mf + ms) + wf + ws) + gamma*(mf + ms + wf + 2*ws)) + \
#        cxi1xi1*(gamma*mf + gamma*ms - b*wf + gamma*wf - 2*b*ws + 2*gamma*ws + 2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi2xi1*(a - beta - b*wf - 2*b*ws + gamma*ws - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi2xi2dt = beta*wf + gamma*mf*ws + gamma*ms*ws + gamma*wf*ws + 2*gamma*ws**2 + 2*cxi3xi2*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + \
#        2*cxi4xi2*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + wf*(a - b*(delta*(mf + ms) + wf + ws)) + \
#        4*ws*(a - b*(delta*(mf + ms) + wf + ws)) + 2*cxi1xi2*\
#         (gamma*mf + gamma*ms - b*wf + gamma*wf - 2*b*ws + 2*gamma*ws + 2*(a - b*(delta*mf + delta*ms + wf + ws))) + mu*wf*xi + \
#        2*cxi2xi2*(a - beta - b*wf - 2*b*ws + gamma*ws - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi2xi3dt =  cxi1xi2*(b*ms - gamma*ms) + cxi2xi2*(b*ms - gamma*ms) + cxi4xi2*(beta + b*delta*ms - gamma*ms) - gamma*ms*ws + \
#        cxi3xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + \
#        cxi3xi2*(-a - gamma*mf + b*delta*ms - 2*gamma*ms - mu - gamma*wf - gamma*ws + b*(delta*mf + delta*ms + wf + ws)) + \
#        cxi1xi3*(gamma*mf + gamma*ms - b*wf + gamma*wf - 2*b*ws + 2*gamma*ws + 2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi2xi3*(a - beta - b*wf - 2*b*ws + gamma*ws - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi2xi4dt = cxi1xi2*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi2*(-(b*mf) - 2*b*ms + gamma*ms) + gamma*ms*ws + \
#        cxi3xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + \
#        cxi3xi2*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + gamma*ws + \
#           2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi1xi4*(gamma*mf + gamma*ms - b*wf + gamma*wf - 2*b*ws + 2*gamma*ws + 2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi4xi2*(a - beta - b*delta*mf - 2*b*delta*ms + gamma*ms - b*(delta*mf + delta*ms + wf + ws) - mu*xi) + \
#        cxi2xi4*(a - beta - b*wf - 2*b*ws + gamma*ws - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi3xi1dt = cxi1xi1*(b*ms - gamma*ms) + cxi2xi1*(b*ms - gamma*ms) + cxi4xi1*(beta + b*delta*ms - gamma*ms) + gamma*ms*ws + \
#        cxi2xi3*(beta + b*ws - gamma*ws) + cxi3xi3*(b*delta*ws - gamma*ws) + cxi4xi3*(b*delta*ws - gamma*ws) + \
#        cxi1xi3*(-a - gamma*mf - gamma*ms - mu - gamma*wf + b*ws - 2*gamma*ws + b*(delta*mf + delta*ms + wf + ws)) + \
#        cxi3xi1*(-a - gamma*mf + b*delta*ms - 2*gamma*ms - mu - gamma*wf - gamma*ws + b*(delta*mf + delta*ms + wf + ws))
# dcxi3xi2dt =  cxi1xi2*(b*ms - gamma*ms) + cxi2xi2*(b*ms - gamma*ms) + cxi4xi2*(beta + b*delta*ms - gamma*ms) - gamma*ms*ws + \
#        cxi3xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + \
#        cxi3xi2*(-a - gamma*mf + b*delta*ms - 2*gamma*ms - mu - gamma*wf - gamma*ws + b*(delta*mf + delta*ms + wf + ws)) + \
#        cxi1xi3*(gamma*mf + gamma*ms - b*wf + gamma*wf - 2*b*ws + 2*gamma*ws + 2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi2xi3*(a - beta - b*wf - 2*b*ws + gamma*ws - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi3xi3dt =  beta*mf + gamma*mf*ms + 2*gamma*ms**2 + 2*cxi1xi3*(b*ms - gamma*ms) + 2*cxi2xi3*(b*ms - gamma*ms) + \
#        2*cxi4xi3*(beta + b*delta*ms - gamma*ms) + ms*mu + gamma*ms*wf + gamma*ms*ws + \
#        2*cxi3xi3*(-a - gamma*mf + b*delta*ms - 2*gamma*ms - mu - gamma*wf - gamma*ws + b*(delta*mf + delta*ms + wf + ws)) + \
#        ms*(a - b*(delta*(mf + ms) + wf + ws))
# dcxi3xi4dt = -(beta*mf) + cxi1xi4*(b*ms - gamma*ms) + cxi2xi4*(b*ms - gamma*ms) + cxi4xi4*(beta + b*delta*ms - gamma*ms) + \
#        cxi1xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi3*(-(b*mf) - 2*b*ms + gamma*ms) + \
#        cxi3xi4*(-a - gamma*mf + b*delta*ms - 2*gamma*ms - mu - gamma*wf - gamma*ws + b*(delta*mf + delta*ms + wf + ws)) - \
#        ms*(2*a + gamma*(mf + 2*ms + wf + ws) - 2*b*(delta*(mf + ms) + wf + ws)) + \
#        cxi3xi3*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + gamma*ws + \
#           2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi4xi3*(a - beta - b*delta*mf - 2*b*delta*ms + gamma*ms - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi4xi1dt = cxi1xi1*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi1*(-(b*mf) - 2*b*ms + gamma*ms) - gamma*ms*ws + cxi2xi4*(beta + b*ws - gamma*ws) + \
#        cxi3xi4*(b*delta*ws - gamma*ws) + cxi4xi4*(b*delta*ws - gamma*ws) + \
#        cxi1xi4*(-a - gamma*mf - gamma*ms - mu - gamma*wf + b*ws - 2*gamma*ws + b*(delta*mf + delta*ms + wf + ws)) + \
#        cxi3xi1*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + gamma*ws + \
#           2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi4xi1*(a - beta - b*delta*mf - 2*b*delta*ms + gamma*ms - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi4xi2dt = cxi1xi2*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi2*(-(b*mf) - 2*b*ms + gamma*ms) + gamma*ms*ws + \
#        cxi3xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + \
#        cxi3xi2*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + gamma*ws + \
#           2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi1xi4*(gamma*mf + gamma*ms - b*wf + gamma*wf - 2*b*ws + 2*gamma*ws + 2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi4xi2*(a - beta - b*delta*mf - 2*b*delta*ms + gamma*ms - b*(delta*mf + delta*ms + wf + ws) - mu*xi) + \
#        cxi2xi4*(a - beta - b*wf - 2*b*ws + gamma*ws - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi4xi3dt = -(beta*mf) + cxi1xi4*(b*ms - gamma*ms) + cxi2xi4*(b*ms - gamma*ms) + cxi4xi4*(beta + b*delta*ms - gamma*ms) + \
#        cxi1xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi3*(-(b*mf) - 2*b*ms + gamma*ms) + \
#        cxi3xi4*(-a - gamma*mf + b*delta*ms - 2*gamma*ms - mu - gamma*wf - gamma*ws + b*(delta*mf + delta*ms + wf + ws)) - \
#        ms*(2*a + gamma*(mf + 2*ms + wf + ws) - 2*b*(delta*(mf + ms) + wf + ws)) + \
#        cxi3xi3*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + gamma*ws + \
#           2*(a - b*(delta*mf + delta*ms + wf + ws))) + \
#        cxi4xi3*(a - beta - b*delta*mf - 2*b*delta*ms + gamma*ms - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
# dcxi4xi4dt =  beta*mf + gamma*mf*ms + 2*gamma*ms**2 + 2*cxi1xi4*(-(b*mf) - 2*b*ms + gamma*ms) + 2*cxi2xi4*(-(b*mf) - 2*b*ms + gamma*ms) + \
#        gamma*ms*wf + gamma*ms*ws + mf*(a - b*(delta*(mf + ms) + wf + ws)) + 4*ms*(a - b*(delta*(mf + ms) + wf + ws)) + \
#        2*cxi3xi4*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + gamma*ws + \
#           2*(a - b*(delta*mf + delta*ms + wf + ws))) + mf*mu*xi + \
#        2*cxi4xi4*(a - beta - b*delta*mf - 2*b*delta*ms + gamma*ms - b*(delta*mf + delta*ms + wf + ws) - mu*xi)
	

# dwsdt
# dwfdt
# dmsdt
# dmfdt
# dcxi1xi1dt 
# dcxi1xi2dt 
# dcxi1xi3dt 
# dcxi1xi4dt 
# dcxi2xi1dt 
# dcxi2xi2dt 
# dcxi2xi3dt 
# dcxi2xi4dt 
# dcxi3xi1dt 
# dcxi3xi2dt 
# dcxi3xi3dt 
# dcxi3xi4dt 
# dcxi4xi1dt 
# dcxi4xi2dt 
# dcxi4xi3dt 
# dcxi4xi4dt 
                