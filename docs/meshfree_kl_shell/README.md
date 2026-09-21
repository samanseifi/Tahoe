# Meshfree RKPM Kirchhoff–Love shell (`meshfree_kl_shell`)

Implementation of the rotation-free, meshfree Kirchhoff–Love shell of Wang & Bazilevs,
*A general-purpose meshfree Kirchhoff–Love shell formulation*, Eng. Comput. 41 (2025) 1379–1410
(GitHub epic #59). Element `RKShellT` in `development/src/elements/meshfree_kl_shell/`.

## What it reuses and what is its own
- Reproducing-kernel shape functions and their first, second and third derivatives come from Tahoe's
  meshfree core: `MLSSolverT` with the cubic B-spline window (`CubicSplineWindowT`, the paper's
  kernel) or the legacy Gaussian window.
- Everything else is element-local: PCA local charts, KL kinematics (paper Eqs. 39–49,
  `KLShellKernels.h`), naturally stabilized nodal integration (Eq. 33), 3-point through-thickness
  Gauss integration, plane-stress J2 plasticity with linear/saturation hardening
  (`PlaneStressJ2.h`), Flanagan–Taylor co-rotation and thickness update (`FlanaganTaylor.h`),
  self-contact and rotational penalty coupling of patches. It is not built on `MeshFreeSupportT`
  or the EFG/RKPM solid elements.
- Essential BCs on the non-interpolatory field: `kinematic_BC ... collocation="1"` (direct nodal
  collocation, `CollocationKBCT`, #70) and the force controller `penalty_displacement_meshfree`
  (`MFPenaltyDisplacementT`): penalty enforcement of the physical displacement over the whole
  support with the physical reaction reported. Both talk to the element through
  `MeshFreeCollocationSupportT`.

## Deck parameters (element block `meshfree_kl_shell`)
`shell_thickness`, `Young_modulus`, `Poisson_ratio`, `density`; RK basis `support_factor`
(normalized support, paper 2.4–3.0), `completeness` (2 quadratic, 3 cubic), `kernel` (1 cubic
B-spline, 0 Gaussian); stabilization `stab_natural` (Eq. 33 Taylor, default 1), `stab_siggrad`
(stress-gradient, default 1), `stab_bending` (curvature, default 0, explicit-CFL sensitive),
`stab_scni`, `smoothed_gradient`; plasticity `yield_stress`, `hardening_modulus`,
`yield_saturation`, `saturation_rate`; kinematics `finite_strain`, `thickness_update`,
`corotational`; loads `load_x/y/z` (per area), `damping` (mass-proportional); self-contact
`contact_stiffness`, `contact_r_in/out`; patch coupling `penalty_coupling`, `couple_node_ID_a/b`;
diagnostics `monitor_node`, `monitor_stride`, `monitor_count`. The explicit force pass is
OpenMP-parallel (`OMP_NUM_THREADS`).

## Tests and benchmarks
- Unit tests: `tests/meshfree/test_KLShell*.cpp`, `test_RKShellConstitutive.cpp`, `test_Collocation`.
- `benchmark_XML/level.0/meshfree_kl_shell/`: elastic obstacle course (Scordelis–Lo, hemisphere,
  pinched cylinder, #66) and the essential-BC regressions `collocation_bc.xml` (#70),
  `penalty_displacement_bc.xml` (explicit) and `penalty_displacement_static.xml` (implicit, one
  Newton iteration) (#73).
- `benchmark_XML/level.2/meshfree_kl_shell/`: shortened, mass-scaled elasto-plastic paper cases
  `necking.xml` (§4.3) and `pinch_plastic.xml` (§4.4), ~45 s each single-threaded (#73). They are
  regressions of the plastic path; the validation against the paper is `fig18_validation/`.
- Run a level: `cd benchmark_XML/level.2 && printf "run.batch\nquit\n" | ../../build/bin/tahoe`
  then the same with `../../build/bin/compare`; or `./run_benchmarks.sh level.0 level.2` from the
  repository root.

## Validation status
| case (paper section) | result |
|---|---|
| Scordelis–Lo roof (4.2.1) | converges to −0.292 with the paper kernel; coarse 225-node deck is dilation-sensitive |
| pinched hemisphere (4.2.2) | 0.0909 vs 0.0924 at M40 (coarse meshes are far too stiff) |
| linear pinched cylinder (4.2.3) | 0.95 of 1.8248e-5 at 8,480 nodes, quadratic completeness |
| necking cylinder (4.3) | peak within 7%; post-peak localizes like the paper's unstabilized case |
| elasto-plastic pinched cylinder (4.4, Fig 18) | **matches** once the literature's quarter-model force (F/4) is accounted for; see `fig18_validation/README.md` |

Method notes: `METHODS_vs_PAPER.md`, `DEVIATIONS_FROM_PAPER.md`,
`IMPLEMENTATION_sigma_xi_and_multibody.md`. `history/` keeps the superseded Fig 18 investigation.
