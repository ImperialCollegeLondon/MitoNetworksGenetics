# MitoNetworksGenetics

Code associated with the article "Mitochondrial network fragmentation modulates mutant mtDNA accumulation independently of absolute fission-fusion rates". See
https://www.biorxiv.org/content/early/2018/09/05/409128.article-metrics

## Installation

Clone the repository

```
$ git clone --recursive https://github.com/ImperialCollegeLondon/MitoNetworksGenetics.git
```

to make sure you get the mitonetworks repository too. Much of the python code depends on the `mitonetworks` package which contains a set of helper functions and classes. To install this package, perform the following commands

```
$ cd mitonetworks
$ python setup.py sdist 
$ pip install ./dist/mitonetworks-0.0.1.tar.gz
```

## Figures

### Alternative_control_laws

- `Param_sweeps` : Parameter sweeps for a set of replication/degradation rate functions, Figures S2I-P


### Linear_Feedback_Control

- ` Ablate_fusion_fission` : Figure S2B
- `Deterministic_phase_portrait` : Figures 2A,B and Figure S2A
- `Quality_control/sel_fus/Analysis`: Figure 3A
- `Quality_control/QC_sweep`: Figure S4
- `Quality_control/sel_deg/Analysis` : Figure 3B
- `Stochastic_analysis/delta_sweep/Post_processing/Analysis`: Figure 2G
- `Stochastic_analysis/kappa_sweep/Analysis`: Figure 2F
- `Stochastic_analysis/mu_sweep/Analysis`: Figure 2E
- `Stochastic_analysis/network_and_b_sweep/Post_processing/Analysis`: Figure S2D,E,F,H and Figure 2H
- `Stochastic_analysis/Nominal_parametrization/Analysis` : Figure 2C,D and Figure S2C
- `Stochastic_analysis/xi_sweep/Analysis` : Figure S2G

### Moran_infinite_sites

 - Exploration of the infinite sites Moran model, Figure B1.

### Simple_moran

 - Exploration of a biallelic Moran model, Figure S3

### Simple_Network_Process

- Deterministic treatment of a simple network process, Figure S1

## Proofs

The following table summarizes the locations of the Mathematica proofs in the paper. LFC = linear feedback control. 

| Formula | Description | Directory|
|:-------:|:-----------:|:--------:|
<a href="http://www.codecogs.com/eqnedit.php?latex=\mathbb{V}(h)=f_s&space;\frac{2\mu&space;t}{n}h(1-h)" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\mathbb{V}(h)=f_s&space;\frac{2\mu&space;t}{n}h(1-h)" title="\mathbb{V}(h)=f_s \frac{2\mu t}{n}h(1-h)" /></a> | LFC, Eq.(11) |	`Linear_Feedback_Control/Proof_Vh` |
<a href="http://www.codecogs.com/eqnedit.php?latex=\mathbb{V}(h)=f_s&space;\frac{2\mu&space;t}{n}h(1-h)" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\mathbb{V}(h)=f_s&space;\frac{2\mu&space;t}{n}h(1-h)" title="\mathbb{V}(h)=f_s \frac{2\mu t}{n}h(1-h)" /></a>  | Alternative_control_laws, Eq.(11) | `Alternative_control_laws/Vh_proofs`
<a href="http://www.codecogs.com/eqnedit.php?latex=\mathbb{V}(h)=&space;\frac{2\lambda&space;t}{n}h(1-h)" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\mathbb{V}(h)=&space;\frac{2\lambda&space;t}{n}h(1-h)" title="\mathbb{V}(h)= \frac{2\lambda t}{n}h(1-h)" /></a> | Alternative_control_laws, Eq.(86) | `Alternative_control_laws/Vh_proofs`
<a href="http://www.codecogs.com/eqnedit.php?latex=\mathbb{V}(h)=f_s&space;\frac{2\mu&space;t}{n}h(1-h)" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\mathbb{V}(h)=f_s&space;\frac{2\mu&space;t}{n}h(1-h)" title="\mathbb{V}(h)=f_s \frac{2\mu t}{n}h(1-h)" /></a> | <a href="http://www.codecogs.com/eqnedit.php?latex=X_S&space;\rightarrow&space;X_S&space;&plus;&space;X_S" target="_blank"><img src="http://latex.codecogs.com/gif.latex?X_S&space;\rightarrow&space;X_S&space;&plus;&space;X_S" title="X_S \rightarrow X_S + X_S" /></a> , Eq.(11)| `Alternative_model_vh_proof`
<a href="http://www.codecogs.com/eqnedit.php?latex=\dot{\mathbf{x}}=0" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\dot{\mathbf{x}}=0" title="\dot{\mathbf{x}}=0" /></a> | Deterministic steady state of LFC | `Linear_Feedback_Control/Deterministic_Steady_State`
