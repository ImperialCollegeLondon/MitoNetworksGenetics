# Stochastic analysis of linear feedback control: delta sweep

This contains code designed to be run on the Imperial High Performance Computing facility around 09/17 (job quotas have since changed)

## Workflow

To make data, run steps 1-6. To use pre-computed data, go straight to step 6.

1. Run `./Make_sweep_vals/delta_SS_lines.ipynb`
2. Run `./Make_submission/Make_submission_param_sweep.ipynb`
2. Copy output files onto HPC
3. On HPC, compile C code as `gcc -o3 closed_loop_deg_single.c -lm -o closed_loop_deg_single.ce`, then `bash master.sh` will submit stochastic sims
4. When done, copy the file `./Post_processing/process_sweep.pbs` and `./Post_processing/Process_stoch_sim.py` onto HPC and then `qsub process_sweep.pbs`
5. When done, this will generate the contents of the `./Post_processing/Data` folder. Copy from cluster into `./Post_processing/Data` folder
6. To generate plots, run the notebook `./Post_processing/Analysis/Analyse_param_sweep.ipynb`