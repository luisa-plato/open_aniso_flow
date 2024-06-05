# Anisotropic flow simulation

This code simulates a coupled, anistropic Navier-Stokes-Nernst-Planck-Poisson system. 
To give it try set up a conda environment according to the yml-File and type

`python run_flow_simulation.py`

into the terminal.

## The mathematical model 

The coupled, anistropic Navier-Stokes-Nernst-Planck-Poisson system

```math
\begin{align*}
    \partial_t \boldsymbol{v} + (\boldsymbol{v} \cdot \nabla) \boldsymbol{v} - \nu \Delta \boldsymbol{v} + \nabla p &= B (c^+ - c^-) \nabla \psi \quad
    & &\text{ in } (0,T) \times \Omega\\
    \partial_t c^{\pm} + \nabla \cdot (c^{\pm} \boldsymbol{v}) - \nabla \cdot (\boldsymbol{\Lambda}(\boldsymbol{d}) (\mu \nabla c^{\pm} \pm F c^{\pm} \nabla \psi)) &= 0 \quad
    & &\text{ in } (0,T) \times \Omega\\
    \nabla \cdot (\boldsymbol{\mathcal{E}}(\boldsymbol{d}) \nabla \psi ) &= B( c^+ - c^-) \quad
    & &\text{ in } (0,T) \times \Omega
\end{align*}
```
is implemented, where $\boldsymbol{\Lambda}(\boldsymbol{d}) = \boldsymbol{I} + \lambda \boldsymbol{d} \otimes \boldsymbol{d}$, 
$\boldsymbol{\mathcal{E}} = \boldsymbol{I} + \varepsilon \boldsymbol{d} \otimes \boldsymbol{d}$ for $\lambda, \varepsilon > 0$,

- $\boldsymbol{v}$ denotes the fluid's velocity field and $p$ the pressure,
- $c^{\pm}$ are the densities of the positively and negatively charged particles, 
- $\psi$ is the electric potential,
- $\boldsymbol{d}\,$ is the so-called director which encodes the anisotropy.

The system is equipped with the boundary conditions
```math
\begin{align*}
    \boldsymbol{v} = 0, \quad
    (\boldsymbol{\Lambda}(\boldsymbol{d}) (\mu \nabla c^{\pm} \pm F c^{\pm} \nabla \psi)) \cdot \boldsymbol{n} = 0, \quad \text{ and } \quad
    (\boldsymbol{\mathcal{E}}(\boldsymbol{d}) \nabla \psi ) \cdot \boldsymbol{n} = 0 \quad \text{ on } [0,T] \times \partial \Omega.
\end{align*}
```
and solved through a linearization and fixed-point iteration in `run_experiment.py`.

All other appearing parameters are constants, which can be chosen in the simulation,
for example by typing

`python run_flow_simulation.py --F 10 --nu 2`.

The direcor $\boldsymbol{d}$ can be chosen from a number of predefined vecotr fields,
which are specified in `problems/director.py`, e.g. by

`python run_flow_simulation.py --director source_and_saddle`.

The generated data is stored in the data folder.


