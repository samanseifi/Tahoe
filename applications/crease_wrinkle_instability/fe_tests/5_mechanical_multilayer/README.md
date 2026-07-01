# Case 5 — mechanical multilayer (bilayer)
Bonded two-layer film under lateral compression, interfacial surface energy off/on.
- `bilayer_intOFF.io1.exo`  γ_interface = 0.
- `bilayer_intON.io1.exo`   γ_interface > 0 (raises onset strain).
Feeds Figure `case3_fe` via `scripts/make_case3_validation.py`.
Mesh: `../../meshes/bilayer_*.geom`; generator `scripts/generate_bilayer_2D.py`.
