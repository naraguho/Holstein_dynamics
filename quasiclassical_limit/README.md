# Quasi-classical limit

This directory contains the original numerical setup for the one-dimensional
quasi-classical Holstein model. It generates nonequilibrium trajectories of the
coupled electronic and lattice system used to investigate charge-density-wave
formation, kink dynamics, and coarsening.

## Programs

- `quasiclassical_demonstration.jl`: adiabatic quasi-classical demonstration
  in which half-filled electronic occupations follow the lattice coordinates.
- `semi-classical.jl`: CPU implementation with instantaneous electronic
  diagonalization and Langevin lattice dynamics.
- `langevin_holstein_lineaer_cuda_v1.jl`: CUDA implementation of the coupled
  electron-lattice Langevin dynamics.
- `data_generation`: existing placeholder retained from the original setup.

The simulation programs and their numerical parameters have not been modified
during the repository reorganization.
