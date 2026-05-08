# Implicit Coulomb friction — work in progress (#40)

Companion to the explicit friction benchmark in
[`../explicit_benchmark/vectorized_cubes_friction.xml`](../explicit_benchmark/vectorized_cubes_friction.xml)
(issue #26).  Same two-cube geometry, but driven by quasi-static
Newton-Raphson via `<updated_lagrangian>` instead of central-difference.

| File | Notes |
|------|-------|
| [`sliding_cubes.xml`](sliding_cubes.xml) | Top face compressed (`u_z = -0.005`) and dragged tangentially (`u_x = 0 → 0.01` over 100 quasi-static steps). |
| [`two_cubes.geom`](two_cubes.geom) | Copy of `explicit_benchmark/geometry/two_cubes_refined.geom`. |

## Status: diverges as expected (#40 LHS-tangent gap)

PR #41 landed the **residual-only** implicit Coulomb branch in
`PenaltyContact3DT`.  `RHSDriver` correctly applies
`f_t = -μ·|f_n|·Δu_t / sqrt(|Δu_t|² + ε²)` from per-striker slip history;
`CloseStep` snapshots the contact-pair geometry; branch selection between
explicit (#26) and implicit (#40) is automatic via `Field().Order()`.

What's missing: the analytical friction tangent in `LHSDriver`,

```
∂f_t/∂Δu_t = s · (I − Δu_t ⊗ Δu_t / smooth²)
∂f_t/∂|f_n| coupling through the normal-spring linearisation
```

Without it, Newton has no way to predict how to balance the friction
perturbation it sees in the residual.  Even at small μ and tiny per-step
slip the iteration kicks out of the basin after iter 0:

```
Step: 1 of 100
   0: Relative error =    8.13e-01
   1: Relative error =    2.68e-01
   2: Relative error =    3.03e+03   ← diverges
```

This XML is kept as the natural integration test for the follow-up
that adds the LHS tangent — running it before/after will give immediate
feedback that the tangent is wired correctly.

## Running

```bash
cd benchmark_XML/level.5/implicit_friction
../../../build/bin/tahoe -f sliding_cubes.xml   # diverges today; will converge once LHS tangent lands
```

## See also
- Issue #40 — implicit Coulomb friction tracking
- PR #41 — residual + slip-history scaffolding
- `tahoe/src/elements/contact/PenaltyContact3DT.cpp` — inline comment in the friction block describes the missing tangent
