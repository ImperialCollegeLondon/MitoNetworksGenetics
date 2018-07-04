# Parameter sweep of mitophagy rate on linear feedback control model

To make data:

1. Run LFC_mu_sweep.ipynb. Copy contents onto HPC. 
2. `$ bash master.sh` on HPC
3. When complete, `qsub process_sweep.pbs`
4. Copy data into the Data directory

To make plots run the Jupyter notebook in ./Analysis