#GSL_HOME="/usr/include/gsl"
GSL_HOME="/home/juvid/gsl"
LD_LIBRARY_PATH="/home/juvid/gsl/lib"
export LD_LIBRARY_PATH
gcc -Wall -std=gnu99 -o infinite_alleles_mtDNA_N_sweep.ce  infinite_alleles_mtDNA_N_sweep.c -I$GSL_HOME/include -L$GSL_HOME/lib -lgsl -lgslcblas -lm
./infinite_alleles_mtDNA_N_sweep.ce 413032 0.023 5.6e-7 100
./infinite_alleles_mtDNA_N_sweep.ce 218360 0.023 5.6e-7 500
./infinite_alleles_mtDNA_N_sweep.ce 270161 0.023 5.6e-7 1000
./infinite_alleles_mtDNA_N_sweep.ce 79778 0.023 5.6e-7 1500
./infinite_alleles_mtDNA_N_sweep.ce 814137 0.023 5.6e-7 2000
