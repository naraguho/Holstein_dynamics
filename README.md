# Holstein Dynamics

Simulation code for nonequilibrium dynamics of the Holstein model.

## Implemented dynamics

- [`quasiclassical_limit/`](quasiclassical_limit/) contains the existing
  one-dimensional quasi-classical Holstein-model setup. The original Julia
  programs and their numerical parameters are preserved unchanged.
- [`ehrenfest_dynamics/`](ehrenfest_dynamics/) contains the GPU-accelerated
  Holstein–Ehrenfest simulation for interaction-quench dynamics.

Additional dynamical limits can be added as sibling directories without
changing the quasi-classical implementation.

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
