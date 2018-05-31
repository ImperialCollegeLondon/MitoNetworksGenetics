#GSL_HOME="/usr/include/gsl"
GSL_HOME="/home/juvid/gsl"
LD_LIBRARY_PATH="/home/juvid/gsl/lib"
export LD_LIBRARY_PATH
gcc -Wall -std=gnu99 -o infinite_alleles_mtDNA_fs_sweep.ce  infinite_alleles_mtDNA_fs_sweep.c -I$GSL_HOME/include -L$GSL_HOME/lib -lgsl -lgslcblas -lm
./infinite_alleles_mtDNA_fs_sweep.ce 885177 0.023 5.6e-7 0
./infinite_alleles_mtDNA_fs_sweep.ce 216620 0.023 5.6e-7 1
./infinite_alleles_mtDNA_fs_sweep.ce 988271 0.023 5.6e-7 2
./infinite_alleles_mtDNA_fs_sweep.ce 970642 0.023 5.6e-7 3