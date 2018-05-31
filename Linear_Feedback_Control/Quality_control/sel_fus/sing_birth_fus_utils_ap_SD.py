
import numpy as np
import pdb
import re
from scipy.integrate import ode 

# Uses alternative parametrization of the linear feedback control



def two_dim_network(t, y, beta, gamma, kappa, b, xi, mu, delta, Q):
	ws, wf, ms, mf = y
	wsdot = -1.0 * gamma * ws * ws - gamma * ws * wf + beta * wf - (mu + b * (kappa - (ws + wf + delta * ms +  delta * mf))) * ws - mu * ws - gamma * mf * ws - gamma * ws * ms
	msdot = -1.0 * gamma * ms * ms - gamma * ms * mf + beta * mf - (mu + b * (kappa - (ws + wf + delta * ms +  delta * mf))) * ms - Q * mu * ms - gamma * wf * ms - gamma * ws * ms
	
	wfdot = 1.0 * gamma * ws * ws + gamma * ws * wf - beta * wf + (mu + b * (kappa - (ws + wf + delta * ms +  delta * mf))) * (2.0*ws + wf) - xi * mu * wf + gamma * mf * ws + gamma * ws * ms
	mfdot = 1.0 * gamma * ms * ms + gamma * ms * mf - beta * mf + (mu + b * (kappa - (ws + wf + delta * ms +  delta * mf))) * (2.0*ms + mf) - xi * mu * mf + gamma * wf * ms + gamma * ws * ms

	return [wsdot, wfdot, msdot, mfdot]

def make_trajectory(ws_init, wf_init, ms_init, mf_init, params, TMAX):
	beta, gamma, kappa, b, xi, mu, delta, Q = params


	y0, t0 = [ws_init, wf_init, ms_init, mf_init], 0.0
	r = ode(two_dim_network).set_integrator('dopri5', nsteps = 10*10)
	r.set_initial_value(y0, t0).set_f_params(beta, gamma, kappa, b, xi, mu, delta, Q)
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
	t_sol = np.array(t_sol); ws_sol = np.array(ws_sol); wf_sol = np.array(wf_sol); ms_sol = np.array(ms_sol); mf_sol = np.array(mf_sol)

	return t_sol, ws_sol,wf_sol,ms_sol,mf_sol

def make_trajectory_plot(axs, ws_init, wf_init, ms_init, mf_init, params, isfirst, TMAX):
	t_sol, ws_sol,wf_sol,ms_sol,mf_sol = make_trajectory(ws_init, wf_init, ms_init, mf_init, params, TMAX)

	ws_final = ws_sol[-1] 
	wf_final = wf_sol[-1]
	ms_final = ms_sol[-1]
	mf_final = mf_sol[-1]

	h_final = (mf_final + ms_final) / float(ws_final + wf_final + ms_final + mf_final)

	
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

	return h_final



def ss_sol_ext_mut(beta, gamma, kappa, b, xi, mu):
	# Steady state solution to network system with selective degradation of mutants where mutant copy
	# number is 0
	ws = (b*beta**2 + b**2*beta*kappa + b*beta*gamma*kappa + b**2*gamma*kappa**2 + 5*b*beta*mu - \
		beta*gamma*mu + b**2*kappa*mu + 3*b*gamma*kappa*mu + 2*b*mu**2 + 2*gamma*mu**2 - b*beta*mu*xi +\
		beta*gamma*mu*xi + b**2*kappa*mu*xi - b*gamma*kappa*mu*xi + 2*b*mu**2*xi - 2*gamma*mu**2*xi - \
		(beta + b*kappa + 2*mu)*np.sqrt(gamma**2*mu**2*(-1 + xi)**2 - 2*b*gamma*mu*(-1 + xi)*(-beta + gamma*kappa + mu + mu*xi)\
		+ b**2*(beta**2 + 2*beta*(gamma*kappa + 3*mu - mu*xi) + (gamma*kappa + mu + mu*xi)**2)))/(2.*(b - gamma)**2*mu*(-1 + xi))
	wf = -((gamma*mu - gamma*mu*xi + b*(-beta + gamma*kappa + mu + mu*xi) + np.sqrt(gamma**2*mu**2*(-1 + xi)**2 - \
		2*b*gamma*mu*(-1 + xi)*(-beta + gamma*kappa + mu + mu*xi) + b**2*(beta**2 + 2*beta*(gamma*kappa + 3*mu - mu*xi) + \
		(gamma*kappa + mu + mu*xi)**2)))*(-(b*(beta**2 + beta*gamma*kappa + 5*beta*mu + 3*gamma*kappa*mu + 2*mu**2)) - \
		gamma*(beta - 2*mu)*mu*(-1 + xi) + b*(beta + gamma*kappa - 2*mu)*mu*xi - b**2*kappa*(beta + gamma*kappa + mu + mu*xi) + \
		beta*np.sqrt(gamma**2*mu**2*(-1 + xi)**2 - 2*b*gamma*mu*(-1 + xi)*(-beta + gamma*kappa + mu + mu*xi) + b**2*(beta**2 + \
		2*beta*(gamma*kappa + 3*mu - mu*xi) + (gamma*kappa + mu + mu*xi)**2)) + b*kappa*np.sqrt(gamma**2*mu**2*(-1 + xi)**2 - \
		2*b*gamma*mu*(-1 + xi)*(-beta + gamma*kappa + mu + mu*xi) + b**2*(beta**2 + 2*beta*(gamma*kappa + 3*mu - mu*xi) + \
		(gamma*kappa + mu + mu*xi)**2)) + 2*mu*np.sqrt(gamma**2*mu**2*(-1 + xi)**2 - 2*b*gamma*mu*(-1 + xi)*(-beta + gamma*kappa + mu + mu*xi) + \
		b**2*(beta**2 + 2*beta*(gamma*kappa + 3*mu - mu*xi) + (gamma*kappa + mu + mu*xi)**2))))/(4.*b*beta*(b - gamma)**2*mu*(-1 + xi))
	return ws, wf

def jacobian_matrix(xi, gamma, beta, kappa, b, mu, delta, Q, ms, mf, ws, wf):
	J=np.array([[-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws,
	b*delta*ws - gamma*ws,
	beta + b*ws - gamma*ws,
	b*delta*ws - gamma*ws],
	[b*ms - gamma*ms,
	-(gamma*mf) + b*delta*ms - 2*gamma*ms - mu - mu*Q - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws,
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

def integrate_system(params, dt):
	xi, beta, gamma, kappa, b, mu, delta, Q, TFinal, ws_init, ms_init, wf_init, mf_init = params

	y0, t0 = np.hstack(([ws_init, wf_init, ms_init, mf_init], np.zeros(16))), 0.0
	r = ode(integrate_SSE_network_sing_birth_fus_SD).set_integrator('vode',method='bdf',nsteps=10**10) # Need to let each iteration take as long as it likes or an error can be thrown
	r.set_initial_value(y0, t0).set_f_params(xi, gamma, beta, kappa, b, mu, delta, Q)

	netsol = []
	t_sol = []
	netsol.append(y0)
	t_sol.append(t0)
	while r.successful() and r.t < TFinal:
		if r.t%100.0<1e-8:
			print(r.t)

		r.integrate(r.t+dt)
		t_sol.append(r.t)
		netsol.append(r.y)
		if r.successful() == False:
			pdb.set_trace()
			print('integration error')
			return -1, np.nan, np.nan

	t_sol = np.array(t_sol)
	netsol = np.array(netsol)

	return 0, t_sol, netsol

def integrate_SSE_network_sing_birth_fus_SD(t, y, xi, gamma, beta, kappa, b, mu, delta, Q):

	ws,wf,ms,mf,cxi1xi1,cxi1xi2,cxi1xi3,cxi1xi4,cxi2xi1,cxi2xi2,cxi2xi3,cxi2xi4,cxi3xi1,cxi3xi2,cxi3xi3,cxi3xi4,cxi4xi1,cxi4xi2,cxi4xi3,cxi4xi4 = y

	dwsdt = beta*wf - gamma*mf*ws - gamma*ms*ws - mu*ws - gamma*wf*ws - (mu + b*(kappa - delta*mf - delta*ms - wf - ws))*ws - gamma*ws**2
	dwfdt = -(beta*wf) + wf*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*mf*ws + gamma*ms*ws + gamma*wf*ws + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws))*ws + gamma*ws**2 - mu*wf*xi
	dmsdt = beta*mf - gamma*mf*ms - gamma*ms**2 - ms*mu*Q - gamma*ms*wf - ms*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - gamma*ms*ws
	dmfdt = -(beta*mf) + gamma*mf*ms + gamma*ms**2 + gamma*ms*wf + mf*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + 2*ms*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ms*ws - mf*mu*xi

	dcxi1xi1dt = beta*wf + 2*cxi1xi1*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + 2*cxi2xi1*(beta + b*ws - gamma*ws) + 2*cxi3xi1*(b*delta*ws - gamma*ws) + 2*cxi4xi1*(b*delta*ws - gamma*ws) + ws*(2*mu + b*(kappa - delta*(mf + ms) - wf - ws) + gamma*(mf + ms + wf + 2*ws))
	dcxi1xi2dt = -(beta*wf) + cxi1xi2*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi2xi2*(beta + b*ws - gamma*ws) + cxi3xi2*(b*delta*ws - gamma*ws) + cxi4xi2*(b*delta*ws - gamma*ws) + cxi3xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi1*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + ws*(-2*mu + 2*b*(-kappa + delta*(mf + ms) + wf + ws) - gamma*(mf + ms + wf + 2*ws)) + cxi2xi1*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi1xi3dt = cxi1xi1*(b*ms - gamma*ms) + cxi2xi1*(b*ms - gamma*ms) + cxi4xi1*(beta + b*delta*ms - gamma*ms) + gamma*ms*ws + cxi1xi3*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi3xi1*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - mu - mu*Q - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi2xi3*(beta + b*ws - gamma*ws) + cxi3xi3*(b*delta*ws - gamma*ws) + cxi4xi3*(b*delta*ws - gamma*ws)
	dcxi1xi4dt = cxi1xi1*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi1*(-(b*mf) - 2*b*ms + gamma*ms) - gamma*ms*ws + cxi1xi4*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi2xi4*(beta + b*ws - gamma*ws) + cxi3xi4*(b*delta*ws - gamma*ws) + cxi4xi4*(b*delta*ws - gamma*ws) + cxi3xi1*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + cxi4xi1*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)
	dcxi2xi1dt = -(beta*wf) + cxi1xi2*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi2xi2*(beta + b*ws - gamma*ws) + cxi3xi2*(b*delta*ws - gamma*ws) + cxi4xi2*(b*delta*ws - gamma*ws) + cxi3xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi1*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi1*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + ws*(-2*mu + 2*b*(-kappa + delta*(mf + ms) + wf + ws) - gamma*(mf + ms + wf + 2*ws)) + cxi2xi1*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi2xi2dt = beta*wf + mu*wf + gamma*mf*ws + gamma*ms*ws + 4*mu*ws + gamma*wf*ws + 2*gamma*ws**2 + b*(kappa - delta*(mf + ms) - wf - ws)*(wf + 4*ws) + 2*cxi3xi2*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + 2*cxi4xi2*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + 2*cxi1xi2*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + mu*wf*xi + 2*cxi2xi2*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi2xi3dt = cxi1xi2*(b*ms - gamma*ms) + cxi2xi2*(b*ms - gamma*ms) + cxi4xi2*(beta + b*delta*ms - gamma*ms) - gamma*ms*ws + cxi3xi2*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - mu - mu*Q - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi3xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi3*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + cxi2xi3*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi2xi4dt = cxi1xi2*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi2*(-(b*mf) - 2*b*ms + gamma*ms) + gamma*ms*ws + cxi3xi2*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + cxi3xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi4*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + cxi4xi2*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi) + cxi2xi4*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi3xi1dt = cxi1xi1*(b*ms - gamma*ms) + cxi2xi1*(b*ms - gamma*ms) + cxi4xi1*(beta + b*delta*ms - gamma*ms) + gamma*ms*ws + cxi1xi3*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi3xi1*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - mu - mu*Q - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi2xi3*(beta + b*ws - gamma*ws) + cxi3xi3*(b*delta*ws - gamma*ws) + cxi4xi3*(b*delta*ws - gamma*ws)
	dcxi3xi2dt = cxi1xi2*(b*ms - gamma*ms) + cxi2xi2*(b*ms - gamma*ms) + cxi4xi2*(beta + b*delta*ms - gamma*ms) - gamma*ms*ws + cxi3xi2*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - mu - mu*Q - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi3xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi3*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi3*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + cxi2xi3*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi3xi3dt = beta*mf + 2*cxi1xi3*(b*ms - gamma*ms) + 2*cxi2xi3*(b*ms - gamma*ms) + 2*cxi4xi3*(beta + b*delta*ms - gamma*ms) + 2*cxi3xi3*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - mu - mu*Q - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + ms*(mu*(1 + Q) + b*(kappa - delta*(mf + ms) - wf - ws) + gamma*(mf + 2*ms + wf + ws))
	dcxi3xi4dt = -(beta*mf) + cxi1xi4*(b*ms - gamma*ms) + cxi2xi4*(b*ms - gamma*ms) + cxi4xi4*(beta + b*delta*ms - gamma*ms) + cxi1xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi3xi4*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - mu - mu*Q - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi3xi3*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + ms*(-2*mu - gamma*(mf + 2*ms + wf + ws) + 2*b*(-kappa + delta*(mf + ms) + wf + ws)) + cxi4xi3*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)
	dcxi4xi1dt = cxi1xi1*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi1*(-(b*mf) - 2*b*ms + gamma*ms) - gamma*ms*ws + cxi1xi4*(-(gamma*mf) - gamma*ms - 2*mu - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) + b*ws - 2*gamma*ws) + cxi2xi4*(beta + b*ws - gamma*ws) + cxi3xi4*(b*delta*ws - gamma*ws) + cxi4xi4*(b*delta*ws - gamma*ws) + cxi3xi1*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + cxi4xi1*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)
	dcxi4xi2dt = cxi1xi2*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi2*(-(b*mf) - 2*b*ms + gamma*ms) + gamma*ms*ws + cxi3xi2*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + cxi3xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi4xi4*(-(b*delta*wf) - 2*b*delta*ws + gamma*ws) + cxi1xi4*(gamma*mf + gamma*ms - b*wf + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) - 2*b*ws + 2*gamma*ws) + cxi4xi2*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi) + cxi2xi4*(-beta + mu - b*wf + b*(kappa - delta*mf - delta*ms - wf - ws) - 2*b*ws + gamma*ws - mu*xi)
	dcxi4xi3dt = -(beta*mf) + cxi1xi4*(b*ms - gamma*ms) + cxi2xi4*(b*ms - gamma*ms) + cxi4xi4*(beta + b*delta*ms - gamma*ms) + cxi1xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi2xi3*(-(b*mf) - 2*b*ms + gamma*ms) + cxi3xi4*(-(gamma*mf) + b*delta*ms - 2*gamma*ms - mu - mu*Q - gamma*wf - b*(kappa - delta*mf - delta*ms - wf - ws) - gamma*ws) + cxi3xi3*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + ms*(-2*mu - gamma*(mf + 2*ms + wf + ws) + 2*b*(-kappa + delta*(mf + ms) + wf + ws)) + cxi4xi3*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)
	dcxi4xi4dt = beta*mf + gamma*mf*ms + 2*gamma*ms**2 + 2*cxi1xi4*(-(b*mf) - 2*b*ms + gamma*ms) + 2*cxi2xi4*(-(b*mf) - 2*b*ms + gamma*ms) + mf*mu + 4*ms*mu + gamma*ms*wf + gamma*ms*ws - b*(mf + 4*ms)*(-kappa + delta*(mf + ms) + wf + ws) + 2*cxi3xi4*(-(b*delta*mf) + gamma*mf - 2*b*delta*ms + 2*gamma*ms + gamma*wf + 2*(mu + b*(kappa - delta*mf - delta*ms - wf - ws)) + gamma*ws) + mf*mu*xi + 2*cxi4xi4*(-beta - b*delta*mf - 2*b*delta*ms + gamma*ms + mu + b*(kappa - delta*mf - delta*ms - wf - ws) - mu*xi)

	deriv = [dwsdt,dwfdt,dmsdt,dmfdt,
			 dcxi1xi1dt,dcxi1xi2dt,dcxi1xi3dt,dcxi1xi4dt,
			 dcxi2xi1dt,dcxi2xi2dt,dcxi2xi3dt,dcxi2xi4dt,
			 dcxi3xi1dt,dcxi3xi2dt,dcxi3xi3dt,dcxi3xi4dt,
			 dcxi4xi1dt,dcxi4xi2dt,dcxi4xi3dt,dcxi4xi4dt]
	return deriv