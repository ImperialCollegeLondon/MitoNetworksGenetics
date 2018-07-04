gcc -o3 closed_loop_deg_single.c -lm -o closed_loop_deg_single.ce
echo "0"
mkdir -p $WORK/networks/sing_birth_fus/mu_sw/p0
qsub submission_0.pbs
sleep 0.5
echo "1"
mkdir -p $WORK/networks/sing_birth_fus/mu_sw/p1
qsub submission_1.pbs
sleep 0.5
echo "2"
mkdir -p $WORK/networks/sing_birth_fus/mu_sw/p2
qsub submission_2.pbs
sleep 0.5
echo "3"
mkdir -p $WORK/networks/sing_birth_fus/mu_sw/p3
qsub submission_3.pbs
sleep 0.5
