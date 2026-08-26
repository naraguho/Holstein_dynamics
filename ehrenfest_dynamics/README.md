# Holstein–Ehrenfest dynamics

This directory contains GPU-accelerated simulation code for nonadiabatic
Ehrenfest dynamics of the semiclassical Holstein chain after an interaction
quench.

## Associated paper

**H. Jang and G.-W. Chern,
"Suppressed coarsening after an interaction quench in the Holstein chain."**

[arXiv:2602.05815](https://arxiv.org/abs/2602.05815)

## Program

- `evolve_100.jl` evolves the electronic one-particle density matrix and the
  classical lattice coordinates with fourth-order Runge–Kutta integration on
  a CUDA GPU.

The supplied setup uses a 1,024-site periodic chain, `dt = 0.01`,
`lambda = 0.6`, and `r = 0.3`. It uses `SLURM_JOB_ID` as the random seed when
available and otherwise defaults to seed 1. Output is currently written to the
`save_base` path defined near the end of the program; update that path for the
target computing environment before running.

## Requirements

- Julia
- CUDA.jl and a CUDA-capable GPU

Run with:

```bash
julia evolve_100.jl
```
