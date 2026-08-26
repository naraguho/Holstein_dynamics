# Holstein Dynamics

Simulation code for nonequilibrium dynamics of the Holstein model.

## Implemented dynamics

- [`quasiclassical_limit/`](quasiclassical_limit/) contains the existing
  one-dimensional quasi-classical Holstein-model setup. The original Julia
  programs and their numerical parameters are preserved unchanged.
- [`ehrenfest_dynamics/`](ehrenfest_dynamics/) contains the GPU-accelerated
  Holstein–Ehrenfest simulation for interaction-quench dynamics.
- [`adiabatic_dynamics/`](adiabatic_dynamics/) contains a two-dimensional
  adiabatic thermal-quench simulation with a classical Langevin bath and its
  dimensionless derivation notes.

Additional dynamical limits can be added as sibling directories without
changing the quasi-classical implementation.

## AI assistance disclosure

This repository was organized and documented with assistance from OpenAI
Codex, an AI coding agent, using scientific code and research materials
provided by Ho Jang.

- The quasi-classical Julia programs were moved without changing their source.
- The Holstein–Ehrenfest program was uploaded without changing its source.
- The adiabatic program was reorganized and engineering-refactored by Codex
  from the supplied `data_classical_bath.jl`; its physical model and default
  simulation parameters were retained.
- Repository structure, README text, portability improvements, and validation
  checks include AI-generated contributions.

AI assistance is disclosed for transparency and does not replace independent
verification of the scientific implementation.

## Associated papers

### Holstein–Ehrenfest dynamics

**H. Jang and G.-W. Chern,
"Suppressed coarsening after an interaction quench in the Holstein chain."**

[arXiv:2602.05815](https://arxiv.org/abs/2602.05815)

### Quasi-classical limit

**H. Jang, Y. Yang, and G.-W. Chern,  
"Anomalous coarsening and nonlinear diffusion of kinks in a one-dimensional quasiclassical Holstein model,"  
Physical Review E (2026).**

DOI: 10.1103/gc41-d1yq
