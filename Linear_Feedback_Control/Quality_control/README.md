# Stochastic analysis of quality control in a linear feedback control

Parameter sweeps for two models of quality control: selective degradation and selective fusion

To make data for any particular control law run 1-5. To just make plots, run 6.

1. Run the submission Jupyter notebook in the first level of the directory (`sel_deg` or `sel_fus`)
2. Copy output files onto HPC
3. On HPC, compile C code as `gcc -o3 *.c -lm -o *.ce` (where * is the name of the `.c` file in the first level directory), then `bash master.sh` will submit stochastic sims
4. When done, copy the file `process_sweep.pbs` and `Process_stoch_sim.py` onto HPC and then `qsub process_sweep.pbs`
5. When done, this will generate the contents of the `Data` folder. Copy from cluster into the `Data` folder
6. To generate plots, run the notebook in the `Analysis` subdirectory