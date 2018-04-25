'''
Phase portrait of a simple network system

Author: Juvid Aryaman
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import ode 

from mitonetworks import utls

np.random.seed(1)

utls.reset_plots()



##############################################################
# Set parameters of model 
##############################################################

n_guess_per_species = 1000
beta = 33.12
gamma = beta / 2.0 /n_guess_per_species # match fusion & fission propensities

rho = 0.023 
eta = 1.1
alpha = 1.04

print(gamma, beta, rho, eta, alpha)

##############################################################
# System definition
##############################################################

# Find steady states given the parameters
sss = ((alpha - eta)*(-beta - rho + alpha*rho))/((-1 + eta)*gamma)
fss = -(((alpha - eta)**2*(-beta - rho + alpha*rho))/((-1 + alpha)*(-1 + eta)*gamma))

print(sss,fss)



def derivative(y, *args):
	s, f = y
	sdot = -1.0*gamma*s*s - gamma * f*s + beta * f +  alpha * rho * s - eta * rho * s 
	fdot = 1.0*gamma*s*s + gamma * f*s - beta * f +  alpha * rho * f - rho * f
	return [sdot, fdot]

def nullcline_sdot(s, *args):
	f = (gamma*s**2)/(beta + rho - alpha*rho - gamma*s)
	return f

def nullcline_fdot(s, *args):
	f = (s*(eta*rho - alpha*rho + gamma*s))/(beta - gamma*s)
	return f

def two_dim_network(t, y, gamma, beta, alpha, eta, rho):
	s = y[0]
	f = y[1]

	sdot = -1.0*gamma*s*s - gamma * f*s + beta * f +  alpha * rho * s - eta * rho * s 
	fdot = 1.0*gamma*s*s + gamma * f*s - beta * f +  alpha * rho * f - rho * f

	return [sdot, fdot]

#####################
# Stability analysis
#####################

# let sdot = F(s,f) and fdot = G(s,f)

dFds = -(gamma*fss) - eta*rho + alpha*rho - 2*gamma*sss
dFdf = beta - gamma*sss
dGds = gamma*fss + 2*gamma*sss
dGdf = -beta - rho + alpha*rho + gamma*sss

J = [[dFds, dFdf],[dGdf, dGds]]

J_eval, J_evec = np.linalg.eig(J)

print('Evals Jacobian:' )
print(J_eval)

#####################
# Phase portrait
#####################

s_lim = 2500
f_lim = 4000
s_space = np.linspace(0.0 , s_lim, 20)
f_space = np.linspace(0.0 , f_lim, 20)


# Create mesh grids
SM, FM = np.meshgrid(s_space, f_space)
SDOTM, FDOTM = np.zeros(SM.shape), np.zeros(FM.shape)

NI, NJ = SM.shape

for i in xrange(NI):
	for j in xrange(NJ):
		s = SM[i,j]
		f = FM[i,j]
		yprime = derivative([s,f],  gamma, beta, alpha, eta, rho)
		SDOTM[i,j] = yprime[0]
		FDOTM[i,j] = yprime[1]

fig, ax = plt.subplots(1,1)
ax.quiver(SM, FM, SDOTM, FDOTM, color = 'r')
ax.plot(sss, fss, 'sk', label = 'Steady State')
ax.plot(0, 0, 'xk')
nc_s_space = np.linspace(0, s_lim, 1000)
nc_s_space = np.linspace(0, f_lim, 1000)
#ax.plot(nc_s_space, nullcline_sdot(nc_s_space,  gamma, beta, alpha, eta, rho), '-.g', label = '$\dot{s} = 0$')
#ax.plot(nc_f_space, nullcline_fdot(nc_f_space,  gamma, beta, alpha, eta, rho), ':k', label = '$\dot{f} = 0$')

#################################
# Plot trajectories
#################################

for i in xrange(20):
	print(i)

	s_init, f_init = np.random.uniform(0,s_lim), np.random.uniform(0,f_lim)
	
	y0, t0 = [s_init, f_init], 0.0
	r = ode(two_dim_network)
	r.set_initial_value(y0, t0).set_f_params( gamma, beta, alpha, eta, rho)
	t1 = 1000
	dt = 1

	s_sol = []
	f_sol = []
	t_sol = []
	t_sol.append(t0)
	s_sol.append(y0[0])
	f_sol.append(y0[1])
	while r.successful() and r.t < t1:
		r.integrate(r.t+dt)
		t_sol.append(r.t)
		s_sol.append(r.y[0])
		f_sol.append(r.y[1])

	t_sol = np.array(t_sol)	
	s_sol = np.array(s_sol)	
	f_sol = np.array(f_sol)
	if i == 0:
		ax.plot(s_sol, f_sol, '-b', alpha = 0.2, label = 'Trajectory')
		ax.plot(s_sol[0], f_sol[0], '.k', label = 'Initial condition')
		ax.plot(s_sol[-1], f_sol[-1], 'xk', label = 'Final condition')
	else:
		ax.plot(s_sol, f_sol, '-b', alpha = 0.2)	
		ax.plot(s_sol, f_sol, '-b', alpha = 0.2)
		ax.plot(s_sol[0], f_sol[0], '.k')
		ax.plot(s_sol[-1], f_sol[-1], 'xk')

ax.set_xlim([-50.0, s_lim])
ax.set_ylim([-50.0, f_lim])
ax.set_xlabel('Singletons, $s$')
ax.set_ylabel('Fused, $f$')
ax.legend(loc = 'lower right', prop = {'size': 8})
plt.tight_layout()
plt.savefig('phase_portrait.svg')
plt.savefig('phase_portrait.png')


###########################################
# Difference in nullclines
###########################################

fig, ax = plt.subplots(1,1)

nc_s_space = np.linspace(sss - 100, sss + 100, 1000)
ax.plot(nc_s_space, nullcline_sdot(nc_s_space, gamma, beta, alpha, eta, rho) - nullcline_fdot(nc_s_space, gamma, beta, alpha, eta, rho), '-k')
ax.plot(nc_s_space, np.zeros(len(nc_s_space)), '--r')
ax.plot(sss, 0, '+k', label = '$s^*$')
ax.set_xlabel('Singletons, $s$')
ax.set_ylabel('Fused nullcline difference,\n $g(s) - l(s)$')

plt.tight_layout()
#plt.savefig('nullcline_difference.svg')
plt.savefig('nullcline_difference.png')
