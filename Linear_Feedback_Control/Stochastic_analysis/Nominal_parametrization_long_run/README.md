# Stochastic analysis of linear feedback control

To make data 

1. Copy `closed_loop_deg_single.c`, `make.sh`, `Process_stoch_sim.py` and the `.pbs` scripts onto the HPC 
3. On HPC, `$ bash make.sh`
4. When done, `qsub process_sweep.pbs`
5. When done, this will generate the contents of the `Data` folder. Copy from HPC into the `Data` folder

To generate plots, run the notebook `./Analysis/Analyse_LFC.ipynb`