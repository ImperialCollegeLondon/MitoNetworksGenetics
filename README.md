# MitoNetworksGenetics

Code associated with [INSERT TITLE] 

[TO DO: Assign figure numbers to readme]

## Installation

Clone the repository

```
$ git clone --recursive https://github.com/ImperialCollegeLondon/MitoNetworksGenetics.git
```

to make sure you get the mitonetworks repository too. Much of the python code depends on the mitonetworks package which contains a set of helper functions and classes. To install this package, perform the following commands

```
$ cd mitonetworks
$ python setup.py sdist 
$ pip install ./dist/mitonetworks-0.0.1.tar.gz
```

## Directories

### Alternative_control_laws

- `Param_sweeps` : Parameter sweeps for a set of replication/degradation rate functions, Figures S2 E-L

### mitonetworks

- A `python` package with various classes and helper functions for analysis and data processing

### Simple_Network_Process

- Deterministic treatment of a simple network process, Figure S1

### Linear_Feedback_Control

- ` Ablate_fusion_fission` : Figure 2D
- `Deterministic_phase_portrait` : Figures 2A-C
- `Quality_control/sel_fus/Analysis`: Figure 3A
- `Quality_control/QC_sweep`: Figure 3A inset
- `Quality_control/sel_deg/Analysis : Figure 3B
- `Stochastic_analysis/Nominal_parametrization` : Figure 2E-H
- `Stochastic_analysis/mu_sweep`: Figure 2I
- `Stochastic_analysis/kappa_sweep`: Figure 2J
- `Stochastic_analysis/delta_sweep`: Figure 2K
- `Stochastic_analysis/network_and_b_sweep/Post_processing/Analysis`: Figure 2L
- `Stochastic_analysis/network_and_b_sweep/Post_processing/Analysis` : Figure S2A,B
- `Stochastic_analysis/xi_sweep/Analysis` : Figure S2C
- `Stochastic_analysis/network_and_b_sweep/Post_processing/Analysis` : Figure S2D

### Moran_infinite_sites

 - Exploration of the infinite sites Moran model, Figure B1.

### Simple_moran

 - Exploration of a biallelic Moran model, Figure S3

### Simple_Network_Process

 - Exploration of a network process with constant rates, Figure XXX

## Proofs

To summarise the locations of the proofs in the paper [To do: complete table]

| Formula | Description | Directory|
|:-------:|:-----------:|:--------:|
|<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{V}(h))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{V}(h)" title="\mathbb{V}(h))" /></a> | LFC |	`Linear_Feedback_Control/Proof_Vh` |
