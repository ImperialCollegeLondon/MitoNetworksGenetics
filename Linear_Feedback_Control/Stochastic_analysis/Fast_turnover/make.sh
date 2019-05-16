mkdir p0
gcc -Wall -o3 closed_loop_deg_single.c -lm -o closed_loop_deg_single.ce
qsub submit_stoch_sim.pbs
