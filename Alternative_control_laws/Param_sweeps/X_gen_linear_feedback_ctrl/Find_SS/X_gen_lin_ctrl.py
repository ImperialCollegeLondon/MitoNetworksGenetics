import numpy as np
import sys
import time
import pandas as pd
from scipy.integrate import ode 

start_time = time.time()
dcols=['t','ws','wf','ms','mf']

def two_dim_network(t, y, params):    
	xi = params['xi']
	gamma = params['gamma']
	beta = params['beta']
	kappa = params['kappa']
	b = params['b']
	mu = params['mu']
	deltas = params['deltas']

	ws, wf, ms, mf = y

	w = ws + wf
	m = ms + mf

	rep_rate = (mu+b*(kappa-(deltas[0]*ws+deltas[1]*wf+deltas[2]*ms+deltas[3]*mf)))

	wsdot = -1.0 * gamma * ws * ws - gamma * ws * wf + beta * wf - rep_rate * ws - mu * ws - gamma * mf * ws - gamma * ws * ms
	msdot = -1.0 * gamma * ms * ms - gamma * ms * mf + beta * mf - rep_rate * ms - mu * ms  - gamma * wf * ms - gamma * ws * ms

	wfdot = 1.0 * gamma * ws * ws + gamma * ws * wf - beta * wf + rep_rate * (2.0*ws + wf) - xi * mu * wf + gamma * mf * ws + gamma * ws * ms
	mfdot = 1.0 * gamma * ms * ms + gamma * ms * mf - beta * mf + rep_rate * (2.0*ms + mf) - xi * mu * mf + gamma * wf * ms + gamma * ws * ms

	return [wsdot, wfdot, msdot, mfdot]

def integrate_odes(t0, ws_init, wf_init, ms_init, mf_init, params, TMAX):
	y0 = [ws_init, wf_init, ms_init, mf_init]
	r = ode(two_dim_network).set_integrator("vode", nsteps=10**10)
	r.set_initial_value(y0, t0).set_f_params(params)
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
		return -1

	d = pd.DataFrame(data=zip(t_sol, ws_sol, wf_sol, ms_sol, mf_sol), columns = dcols)

	return d


grid_num = int(sys.argv[1])
h_target = float(sys.argv[2])
wall_time = float(sys.argv[3]) # in seconds


data = pd.DataFrame(columns = dcols)


#######################################
# Nominal params
#######################################

params = {
    'beta':33.12, 
    'gamma':0.03785142857142858,
    'b':1.2416523075924095e-05, 
    'kappa':11.662903457629223,
    'mu':0.023, 
    'xi':0.0, 
    'deltas':[0.8, 1.0, 0.2, 0.3],# ws, wf, ms, mf
}
gamma = params['gamma']
beta = params['beta']

rat_nominal = gamma/beta

#######################################
# Perturb params
#######################################

n_points = 11
lg_ls_low_churn = -2
lg_ls_high_churn = 2

lg_ls_low_rat = -2
lg_ls_high_rat = 2

mags = 10**np.linspace(lg_ls_low_churn,lg_ls_high_churn,n_points) # magnitudes of the fusion/fission rate
rats = 10**np.linspace(lg_ls_low_rat,lg_ls_high_rat,n_points) # ratios of the fusion/fission rate

i = grid_num/11
j = grid_num%11

mag = mags[i]
rat = rats[j]

newbeta = beta * mag
newgamma = newbeta * rat_nominal * rat

params['beta'] = newbeta
params['gamma'] = newgamma

##################################
# Initial conditions
##################################

w0 = 1000.0 * (1. - h_target)
m0 = 1000.0 * h_target

ws_old = 0.5 * w0
wf_old = 0.5 * w0
ms_old = 0.5 * m0
mf_old = 0.5 * m0

n_old = ws_old + wf_old + ms_old + mf_old
n_new = 0.0



t_old = 0.0
dt = 10.0 # steps of 10 days

#################################
# Integrate
#################################

t_new = 0; ws_new = 0; wf_new = 0; ms_new = 0; mf_new = 0

while True: # keep going till convergence or time is up
	
	update_d = integrate_odes(t_old, ws_old, wf_old, ms_old, mf_old, params, t_old+dt)

	if type(update_d) != int: # update without error
		data = data.append(update_d,ignore_index=True)
	else:
		print('Exit')
		break # there was an integration error

	t_new, ws_new, wf_new, ms_new, mf_new = update_d.iloc[-1]
	n_new = np.sum([ws_new, wf_new, ms_new, mf_new])
	print(n_old, n_new)


	is_converged = ((abs(ws_new - ws_old) < 0.5) and (abs(wf_new - wf_old) < 0.5) and (abs(ms_new - ms_old) < 0.5) and (abs(mf_new - mf_old) < 0.5))
	if ( is_converged or (time.time()-start_time > wall_time*0.9)):
		break

	# Update old -> new
	t_old, ws_old, wf_old, ms_old, mf_old = t_new, ws_new, wf_new, ms_new, mf_new	
	n_old = n_new

data.to_csv('ode_sol_{}.csv'.format(grid_num))