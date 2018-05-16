# Stochastic analysis of alternative feedback control laws

Parameter sweeps for fusion and fission rates for various control laws. Control laws are pre-appended with a letter corresponding to Figure 1 A of [JohnstonJones2016AmJournHumGenet] for controls A-G. Controls X-Z are additional controls which were not considered in [JohnstonJones2016AmJournHumGenet].

To make data for any particular control law run 1-5. To just make plots, run 6.

1. Run the submission Jupyter notebook in the corresponding directory, e.g. `./A_relaxed_rep/Make_submit_sweep_rel_rep.ipynb`
2. Copy output files onto HPC
3. On HPC, compile C code as `gcc -o3 closed_loop_deg_single.c -lm -o closed_loop_deg_single.ce`, then `bash master.sh` will submit stochastic sims
4. When done, copy the file `process_sweep.pbs` and `Process_stoch_sim.py` onto HPC and then `qsub process_sweep.pbs`
5. When done, this will generate the contents of the `Data` folder. Copy from cluster into the `Data` folder
6. To generate plots, run the notebook `./Analyse_all/Analyse_all.ipynb`