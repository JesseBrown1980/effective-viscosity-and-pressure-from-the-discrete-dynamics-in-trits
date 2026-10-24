# Independent Taylor–Green initialization and null-vortex control

**Date:** 2026-10-24
**Source audited:** `JesseBrown1980/Taylor_green_regularity-resolved-and-solved` commit `4b32817cf908e6f6106c3f50acfc28b82847e105`
**C++ source blob SHA:** `6e73b4931eb2e5113762b63b06f3bc02340cf47c`

## Independent test method

A separately implemented Python/NumPy D3Q27 experiment reconstructs:
- the 27 lattice velocities and all **288** collision LUTs;
- 23,328 local collision configurations, each checked for exact mass, momentum, energy preservation and strict involution;
- the orthogonal Taylor–Green initializer from the source commit;
- the same 32-gates-per-cell-per-tick schedule and periodic streaming;
- kinetic-energy, enstrophy, density and vorticity observers.

Two initial states, each on a 24³ lattice, are evolved for 48 transitions (49 recorded ticks):
1. **Taylor–Green:** orthogonal pair synthesis with `U0=0.20`.
2. **Null control:** every one of the 27 occupancies begins at `1` in every cell, so initial velocity and vorticity are exactly zero.

This is an independent reimplementation of the published equations and collision schedule, **not an execution of the original C++ executable or a physical Host8 disk experiment**.

## Previously raised initialization issue: corrected

Independent reconstruction reproduces the reported initial values:

| Property | Result |
|---|---:|
| Cell density | Exactly 27 at all 13,824 cells |
| Total momentum | Exactly (0, 0, 0) |
| Total doubled microscopic energy | Exactly 746,496 |
| Target velocity RMS | 0.100000 |
| Actual velocity RMS | 0.102798321 |
| Relative L² error | 0.279305191 |
| Cosine similarity | 0.962436969 |
| Initial kinetic energy density | 0.142661180 |
| Initial enstrophy | 0.001339989 |

The independent time evolution matches the published Taylor–Green measurements at ticks 1, 2, 3, 6 and 48 to their CSV reporting precision. Momentum and integer microscopic energy remain conserved for both experiments.

## New independent negative control

**All 288 of the current collision LUTs map the local uniform occupancy input `(1,1,1,1)` to a different output.** Thus the uniform, zero-flow initial microstate is not a stationary state of the collision rules.

| Tick | TGV enstrophy | Null enstrophy | TGV kinetic energy | Null kinetic energy |
|---:|---:|---:|---:|---:|
| 0 | 0.001340 | 0.000000 | 0.142661 | 0.000000 |
| 1 | 0.020046 | 0.015117 | 0.670623 | 0.531779 |
| 6 | 0.027536 | 0.017845 | 0.768078 | 0.649401 |
| 48 | 0.026968 | 0.015660 | 0.689470 | 0.642094 |

At tick 1 the Taylor–Green arm's density range is 14–38 particles/cell, while the zero-vortex arm's density range is 17–38. The initially uniform density does not remain uniform. The zero-vortex arm also develops nonzero vorticity (peak 0.39957 at tick 1).

The reported 20.55-fold TGV enstrophy amplification **is a real result of this finite lattice evolution**, but it cannot yet be assigned uniquely to the continuum vortex-stretching term: much of the growth also occurs when the initial macroscopic velocity is zero. At tick 6 the null arm's enstrophy is approximately 64.8% of the TGV arm's enstrophy. This percentage is descriptive, not a rigorous causal decomposition of nonlinear trajectories.

Moreover, macroscopic kinetic energy rises from 0.142661 to 0.689470 over the reported TGV interval. The claim of a net transfer **from macroscopic kinetic into thermal energy** is opposite to the net direction in that interval; the exact microscopic energy conservation remains valid.

## Scientific meaning

**Supported:** Exact finite-state collision invariants, orthogonal initializer fidelity, conservation during evolution, bounded fixed-grid observables, and reproducible finite-grid enstrophy growth.

**Not demonstrated by these data:** Incompressible Navier–Stokes evolution after tick 0, attribution of enstrophy growth specifically to vortex stretching, the continuum viscous kinetic-energy identity, a mesh-independent vorticity estimate, or solution of the Clay global-regularity problem.

Finite occupancy bounds vorticity at fixed spacing `Delta x`; this bound can scale as `1/Delta x` under refinement and is not itself a continuum regularity proof.

## Next validation

Add a compulsory `U0=0` null-vortex control to `reproduce.py`. Check density variability, discrete divergence, coherent Taylor–Green Fourier amplitude, enstrophy production and the actual stretching term for both arms. Investigate a collision operator with a stationary isotropic equilibrium, and compare the result to the current reversible lookup rules.

This is a standalone peer-review audit; the original results and mathematical documentation are left unchanged.
