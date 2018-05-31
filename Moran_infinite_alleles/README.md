# Infinite alleles Moran model for mtDNA

An infinite alleles model with continuous time and Moran steps, where every Moran step has a binomially distributed number of mutation events.

To make data

`$ bash run.sh`

`GSL_HOME` and `LD_LIBRARY_PATH` in `run.sh` should be modified to point to the GSL installation directory on your machine. This is often `/usr/include/gsl` and `/usr/include/gsl/include` respectively, unless the installation was in your home directory. See https://www.gnu.org/software/gsl/ for help on installing GSL.

To make plots, run the Jupyter notebook.