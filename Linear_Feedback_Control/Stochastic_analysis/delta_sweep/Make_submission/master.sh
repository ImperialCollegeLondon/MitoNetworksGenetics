echo "0"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p000
qsub submission_0.pbs
sleep 0.5
echo "1"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p001
qsub submission_1.pbs
sleep 0.5
echo "2"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p002
qsub submission_2.pbs
sleep 0.5
echo "3"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p003
qsub submission_3.pbs
sleep 0.5
echo "4"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p004
qsub submission_4.pbs
sleep 0.5
echo "5"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p005
qsub submission_5.pbs
sleep 0.5
echo "6"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p006
qsub submission_6.pbs
sleep 0.5
echo "7"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p007
qsub submission_7.pbs
sleep 0.5
echo "8"
mkdir -p /work/ja1109/networks/sing_birth_fus/param_sweep/delta_sweep/p008
qsub submission_8.pbs
sleep 0.5
