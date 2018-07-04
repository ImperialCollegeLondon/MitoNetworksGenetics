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

- `Param_sweeps` : Parameter sweeps for a set of replication/degradation rate functions, Figure XXX
- `Vh_proofs` : Heteroplasmy variance proofs for a set of feedback control laws

### mitonetworks

- A `python` package with various classes and helper functions for analysis and data processing

### Simple_Network_Process

- Deterministic treatment of a simple network process, Figure XXX

### Linear_Feedback_Control

- ` Ablate_fusion_fission` : Deterministic analysis reducing the fusion and fission rates, Figure XXX
- `Deterministic_phase_portrait` : Generates Figure XXX
- `Deterministic_Steady_State`: Mathematica notebook to find deterministic steady state
- `Nominal_parametrization` : Python code to find nominal parameterization 
- `Proof_Vh` : Analytic treatment of linear feedback control proving V(h) formula 
- `Quality_control` : Exploration of selective mitophagy and selective fusion, Figure XXX
- `Stochastic_analysis` : Python code to perform parameter sweeps, Figure XXX

 ### Moran_infinite_sites

 - Exploration of the infinite sites Moran model, Figure XXX

 ### Simple_moran

 - Exploration of a biallelic Moran model, Figure XXX

 ### Simple_Network_Process

 - Exploration of a network process with constant rates, Figure XXX

## Proofs

To summarise the locations of the proofs in the paper [To do: complete table]

| Formula | Description | Directory|
|:-------:|:-----------:|:--------:|
|<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{V}(h))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{V}(h)" title="\mathbb{V}(h))" /></a> | LFC |	`Linear_Feedback_Control/Proof_Vh` |