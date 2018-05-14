# MitoNetworksGenetics

Code associated with [INSERT TITLE] 

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

### mitonetworks

- A `python` package with various classes and helper functions for analysis and data processing

### Simple_Network_Process

- Deterministic treatment of a simple network process, Figure XXX

### Linear_Feedback_Control

- `Deterministic_phase_portrait` : Generates Figure XXX
- `Deterministic_Steady_State`: Mathematica notebook to find deterministic steady state
- `Nominal_parametrization` : Python code to find nominal parameterization 
- `Proof_Vh` : Analytic treatment of linear feedback control proving V(h) formula 
- `Stochastic_analysis` : Python code to perform parameter sweeps, Figure XXX

 deterministic steady state and a proof of a stochastic differential equation for heteroplasmy
