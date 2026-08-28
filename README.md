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

## Justification for the semiclassical lattice treatment

The simulations in this repository treat the phonon displacement and momentum
as classical lattice variables rather than as fully quantum-mechanical phonon
operators. Two physical considerations motivate this approximation in the
regimes studied here.

### 1. Phonons are slow degrees of freedom

When the characteristic phonon frequency is small compared with the electronic
energy scale,

$$
\frac{\hbar\Omega}{t_{\mathrm{nn}}} \ll 1,
$$

the lattice evolves slowly relative to the electrons. It can then be treated
as a collective coordinate moving on an electronic energy landscape. In the
adiabatic implementation, the electronic state is recomputed for each
instantaneous lattice configuration. In the Ehrenfest implementation, the
electrons retain their explicit time evolution while the lattice remains
classical.

### 2. Relative quantum fluctuations are small at strong coupling

The electron-phonon coupling produces a characteristic static distortion

$$
X_0 = \frac{g}{K},
$$

whereas the zero-point displacement of a harmonic oscillator is

$$
\Delta X = \sqrt{\frac{\hbar}{2m\Omega}},
\qquad
\Omega = \sqrt{\frac{K}{m}}.
$$

For fixed mass and stiffness, increasing $g$ increases $X_0$ without increasing
$\Delta X$. The relative fluctuation therefore becomes small when

$$
\frac{\Delta X}{X_0} \ll 1.
$$

This is the strong-coupling regime emphasized in much of this repository: the
coupling-induced distortion is large compared with the intrinsic quantum
spread of the oscillator, making a semiclassical lattice description more
natural. The approximation still neglects genuinely quantum-phonon effects
such as zero-point dynamics, tunneling between lattice configurations, and
electron-phonon entanglement; it is expected to be least reliable for fast,
light phonons or weak coupling.

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
