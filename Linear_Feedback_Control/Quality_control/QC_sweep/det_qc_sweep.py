# coding: utf-8
# QC Sweep


import numpy as np
import pdb
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker

import mitonetworks.det as mtd
import mitonetworks.stoch as mts

params = {
    'beta':33.12, 
    'gamma':0.03785142857142858,
    'b':1.2416523075924095e-05, 
    'kappa':11.662903457629223,
    'mu':0.023, 
    'xi':0.0, 
    'delta':1.0
}

gamma_nom = params['gamma']

epsilon_space = np.linspace(-5,-0.5,10)
Q_space = 1.0 + 10**epsilon_space
Q_space = 1./np.hstack(([1], Q_space)) # NB: det code uses opposite convention to avoid performing division


n_points = 21
lg_ls_low_gamma = -2
lg_ls_high_gamma = 2
gamma_space = 10**np.linspace(lg_ls_low_gamma,lg_ls_high_gamma,n_points) # gamma space
gamma_labels = np.log10(gamma_space)

deltah_arr = np.zeros((len(gamma_space),len(Q_space)))

make_data = True

T_measure = 1000.0


def get_delta_h(params):
    # Get ICs
    lc = mts.SubmitHPCSoluableFeedbackControl(ss_definition=mtd.E_linear_feedback_control_ss,nominal_params=params,
                                         hpc_workdir=None,c_filename=None,c_file_param_ordering_convention=None,hpc_jobname=None)
    ics, h, e = lc.find_ss_h_target(lc.nominal_params)
    if h is not np.nan:
        # Simulate ODEs
        lc = mtd.FeedbackControl(network_defn=mtd.E_linear_feedback_control, params=params, initial_state=ics, TMAX = T_measure)   
        initial_h = (ics[2]+ics[3])/float(np.sum(ics))
        lc.make_trajectory()
        final_state = lc.state
        final_h = (final_state[2]+final_state[3])/float(np.sum(final_state))
        delta_h = final_h - initial_h
        return delta_h
    else:
        raise Exception('Error')


if make_data:
    for i, gamma_change in enumerate(gamma_space):
        print(i)
        for j, newQ in enumerate(Q_space):
            # Update params
            newgamma = gamma_nom * gamma_change 
            params['gamma'] = newgamma
            params['Q_f'] = newQ                
            deltah_arr[i,j] = get_delta_h(params)            
        np.save('deltah_arr',deltah_arr)        
np.save('deltah_arr',deltah_arr)      