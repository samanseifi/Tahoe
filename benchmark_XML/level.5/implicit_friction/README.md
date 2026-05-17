# Implicit Coulomb friction (#40)

Companion to the explicit friction benchmark in
[`../explicit_benchmark/vectorized_cubes_friction.xml`](../explicit_benchmark/vectorized_cubes_friction.xml)
(issue #26).  Same two-cube geometry, but driven by quasi-static
Newton-Raphson via `<updated_lagrangian>` instead of central-difference.

| File | Notes |
|------|-------|
| [`sliding_cubes.xml`](sliding_cubes.xml) | Top face compressed (`u_z = -0.005`) and dragged tangentially (`u_x = 0 → 0.01` over 100 quasi-static steps). |
| [`two_cubes.geom`](two_cubes.geom) | Copy of `explicit_benchmark/geometry/two_cubes_refined.geom`. |

## Status: converges

The implicit Coulomb path is complete:

- **Residual + slip history** (PR #41): `PenaltyContact3DT::RHSDriver` applies
  `f_t = -μ·|f_n|·Δu_t / sqrt(|Δu_t|² + ε²)` from per-striker slip history;
  `CloseStep` snapshots the contact-pair geometry; branch selection between
  explicit (#26) and implicit (#40) is automatic via `Field().Order()`.
- **Tangent** (this PR): `PenaltyContact3DT::LHSDriver` builds the friction
  contribution to K by a per-pair 12×12 forward-FD on the same residual,
  with the geometric chain rule (changes in normal and contact pressure
  under DOF perturbation) baked in.

Newton-Raphson converges quadratically:

```
Step: 1 of 100
   0: Relative error = 8.17e-01
   1: Relative error = 2.50e-05
   2: Relative error = 1.64e-08    ← converged

Step: 2 .. 100
   0: Relative error = ~3e-04
   1: Relative error = ~4e-09      ← ~1 iter/step
```

Whole 100-step run completes in ≈ 4 s (serial, profile_matrix solver).

## Running

```bash
cd benchmark_XML/level.5/implicit_friction
../../../build/bin/tahoe -f sliding_cubes.xml
```

## See also
- Issue #40 — implicit Coulomb friction tracking
- PR #41 — residual + slip-history scaffolding (merged)
- `tahoe/src/elements/contact/PenaltyContact3DT.cpp` — inline comment at
  the `LHSDriver` friction-tangent block documents the FD construction
  and the chain rule it captures (normal, contact pressure, slip)
