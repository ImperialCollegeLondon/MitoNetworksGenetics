
# # Process parallelised data and reduce to summary statistics
#
# Stochastic simulations are run in parallel, generating many files. To read in each file and generate copy number and heteroplasmy summary statistics. Mean and variance statistics will be computed online, so this code is scalable.
# This script is to be run for every parametrization


import numpy as np
import pandas as pd
import pdb 
import os
import time
import sys
import pdb

param_index = int(sys.argv[1])
nrowskip = int(sys.argv[2]) # lines of metadata in output
statfilename='online_stats_{:d}.csv'.format(param_index) # filename of output


t0 = time.time()

testing = False

if testing:
    max_iters = 4 # number of files expected
    data_dir = '.'
else:
    max_iters = 1000 # number of files expected
    data_dir = '.'

output_file_form = 'output_{0}_{1}.txt'


nrep_per_job = 100 # number of repeats in each job


## Stream data and compute summary statistics

data1 = pd.read_csv(data_dir+'/'+output_file_form.format(param_index,0), skiprows=nrowskip, delimiter=',')
t_vals = data1[data1['rep']==0]['t'].drop_duplicates().as_matrix()


# In[ ]:

# To be performed on every row of a repeat (so at constant time). This is parallelizable.
def online_mean_var(x): 
    #pdb.set_trace()
    if np.any(pd.isnull(x[['rep','ws','ms','wf','mf','n']])): # if any missing (excl h)
        return x # do nothing, i.e. ignore
    
    ns = x.ws + x.ms # number singletons
    nf = x.n - ns # number fused
    w = x.ws + x.wf # number of wild-types
    m = x.ms + x.mf # number of mutants
    if x.n > 0:
        fs = ns/float(x.n) # fraction singletons
        
        # CDF of heteroplasmy
        if x.h <= 0.1:
            x.cdf_h_0_1 += 1
        if x.h <= 0.2:
            x.cdf_h_0_2 += 1
        if x.h <= 0.3:
            x.cdf_h_0_3 += 1
        if x.h <= 0.4:
            x.cdf_h_0_4 += 1
        if x.h <= 0.5:
            x.cdf_h_0_5 += 1
        if x.h <= 0.6:
            x.cdf_h_0_6 += 1
        if x.h <= 0.7:
            x.cdf_h_0_7 += 1
        if x.h <= 0.8:
            x.cdf_h_0_8 += 1
        if x.h <= 0.9:
            x.cdf_h_0_9 += 1


        # CDF of fraction of singletons
        if fs <= 0.1:
            x.cdf_fs_0_1 += 1
        if fs <= 0.2:
            x.cdf_fs_0_2 += 1
        if fs <= 0.3:
            x.cdf_fs_0_3 += 1
        if fs <= 0.4:
            x.cdf_fs_0_4 += 1
        if fs <= 0.5:
            x.cdf_fs_0_5 += 1
        if fs <= 0.6:
            x.cdf_fs_0_6 += 1
        if fs <= 0.7:
            x.cdf_fs_0_7 += 1
        if fs <= 0.8:
            x.cdf_fs_0_8 += 1
        if fs <= 0.9:
            x.cdf_fs_0_9 += 1

        if x.h==0:
            x.p_h_fix_0 += 1
        if x.h == 1:
            x.p_h_fix_1 += 1
        if fs == 0:
            x.p_fs_fix_0 += 1
        if fs == 1:
            x.p_fs_fix_1 += 1

    if x.n==0:
        x.nextinct+=1
    
    # Update n
    x.counts += 1
    
    # Update delta
    x.delta_n = x.n - x.mean_n
    if x.n > 0:
        x.delta_h = x.h - x.mean_h
        x.delta_fs = fs - x.mean_fs
    x.delta_ns = ns - x.mean_ns
    x.delta_nf = nf - x.mean_nf
    x.delta_w = w - x.mean_w
    x.delta_m = m - x.mean_m
    
    # Update mean
    x.mean_n += x.delta_n/x.counts
    if x.n > 0:
        x.mean_h += x.delta_h/x.counts
        x.mean_fs += x.delta_fs/x.counts
    x.mean_nf += x.delta_nf/x.counts
    x.mean_ns += x.delta_ns/x.counts
    x.mean_w += x.delta_w/x.counts
    x.mean_m += x.delta_m/x.counts
    
    # Update delta2
    x.delta2_n = x.n - x.mean_n 
    if x.n > 0:
        x.delta2_h = x.h - x.mean_h
        x.delta2_fs = fs - x.mean_fs
    x.delta2_ns = ns - x.mean_ns 
    x.delta2_nf = nf - x.mean_nf
    x.delta2_w = w - x.mean_w
    x.delta2_m = m - x.mean_m
    
    # Update M2
    x.M2_n += x.delta_n*x.delta2_n
    if x.n > 0:
        x.M2_h += x.delta_h*x.delta2_h
        x.M2_fs += x.delta_fs*x.delta2_fs
    x.M2_ns += x.delta_ns*x.delta2_ns
    x.M2_nf += x.delta_nf*x.delta2_nf
    x.M2_w += x.delta_w*x.delta2_w
    x.M2_m += x.delta_m*x.delta2_m
    
    return x
        
        
        
# Initialise summary stats
online_cols =['t','counts','nextinct',
              'delta_n','mean_n','delta2_n','M2_n','var_n',
              'delta_h','mean_h','delta2_h','M2_h','var_h',
              'delta_w','mean_w','delta2_w','M2_w','var_w',
              'delta_m','mean_m','delta2_m','M2_m','var_m',
              'delta_ns','mean_ns','delta2_ns','M2_ns','var_ns',
              'delta_nf','mean_nf','delta2_nf','M2_nf','var_nf',
              'delta_fs','mean_fs','delta2_fs','M2_fs','var_fs',
              'cdf_h_0_1','cdf_h_0_2','cdf_h_0_3','cdf_h_0_4','cdf_h_0_5','cdf_h_0_6','cdf_h_0_7','cdf_h_0_8','cdf_h_0_9',
              'p_h_fix_0','p_h_fix_1',
              'cdf_fs_0_1','cdf_fs_0_2','cdf_fs_0_3','cdf_fs_0_4','cdf_fs_0_5','cdf_fs_0_6','cdf_fs_0_7','cdf_fs_0_8','cdf_fs_0_9',
              'p_fs_fix_0','p_fs_fix_1'
             ]




init_online_stats = np.zeros((len(t_vals),len(online_cols)))
init_online_stats[:,0] = t_vals

online_stats = pd.DataFrame(data=init_online_stats, columns=online_cols)




missing_data_counter = 0

for i in range(max_iters):
    if i % 100 == 0:
        # save
        online_stats.to_csv(statfilename,index=False)
        print(i)
    
    if os.path.isfile(data_dir+'/'+output_file_form.format(param_index,i)): # check file exists        
        data = pd.read_csv(data_dir+'/'+output_file_form.format(param_index,i), skiprows=nrowskip, delimiter=',')  

        for j in range(nrep_per_job):
            rep_data = data[data['rep']==j]
            if len(rep_data)>0:                
                temp_join = pd.merge(online_stats, rep_data, on = 't', how = 'left')
            
                temp_join= temp_join.apply(lambda x: online_mean_var(x), axis = 1) # to apply on rows
                # Update
                online_stats[online_cols] = temp_join[online_cols]
                
        
    else:
        missing_data_counter += 1

online_stats['var_n'] = online_stats['M2_n']/(online_stats['counts'] - 1.0)
online_stats['var_h'] = online_stats['M2_h']/(online_stats['counts'] - 1.0)
online_stats['var_fs'] = online_stats['M2_fs']/(online_stats['counts'] - 1.0)
online_stats['var_ns'] = online_stats['M2_ns']/(online_stats['counts'] - 1.0)
online_stats['var_nf'] = online_stats['M2_nf']/(online_stats['counts'] - 1.0)
online_stats['var_w'] = online_stats['M2_w']/(online_stats['counts'] - 1.0)
online_stats['var_m'] = online_stats['M2_m']/(online_stats['counts'] - 1.0)
   


# In[ ]:

online_stats.to_csv(statfilename,index=False)


# In[ ]:

print('Job completed successfully!')
print('Time taken (s):')


# In[ ]:

t1 = time.time()
total = t1-t0
print(total)
print('Missing data:')
print(missing_data_counter)

