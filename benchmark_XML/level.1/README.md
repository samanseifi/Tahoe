# Level 1 — Constitutive model verification

Material-model tests under simple loading paths.  Each subdirectory holds
one or more *material.NN* (or *cohesive.NN*) cases that drive a single
element / patch with prescribed strain so the constitutive response can be
read off against an analytical reference.

Run as the other levels:
```bash
cd benchmark_XML/level.1
printf "run.batch\nquit\n" | ../../build/bin/tahoe
printf "run.batch\nquit\n" | ../../build/bin/compare
```

## Subdirectories

| Directory | Theme |
|-----------|-------|
| `material.solid/{1D,2D,3D}/material.NN/` | Bulk solid materials — linear elastic (`material.01`), elastoplastic (`material.02`), Drucker-Prager (`material.11`), hyperelastic Ogden / Mooney variants, etc.  Each subdir loads its material under a unit cell and checks stress vs analytical. |
| `material.surface/{2D,3D}/cohesive.NN/` | Cohesive-zone surface laws — `cohesive.0X` covers Xu-Needleman, linear, exponential, ductile traction-separation models. |
| `material.surface/{2D,3D}/SIMOD/` | SIMOD-format cohesive parameter conversion. |
| `stress.smoothing/` | Nodal stress recovery / averaging on a Q4 patch. |

## Status (May 2026)

100 / 103 PASS, 3 FAIL.  Remaining failures: `material.solid/{2D,3D}/material.11/mat.11.a.xml`
and `material.solid/2D/material.02/mat.2.a.xml` — unregistered materials in
the factory, tracked under #37.
