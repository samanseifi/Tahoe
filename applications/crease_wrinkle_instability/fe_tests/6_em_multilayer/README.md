# Case 6 — electromechanical multilayer
Electromechanical loading of a bonded multilayer DE (the partner to case 5).
No runs yet. Setup: bilayer mesh (`../../meshes/bilayer_*.geom`) + the two-stage
EM workflow of case 2 (static pre-strain → explicit voltage ramp) applied to the
laminate, with surface tension on the top face and the buried interface.
Use `common/run_sweep2.py` as the deck template (add the second block + interface
side-set) once the EM bilayer study is run.
