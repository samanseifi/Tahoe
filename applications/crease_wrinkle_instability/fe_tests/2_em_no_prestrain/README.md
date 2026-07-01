# Case 2 — electromechanical, no pre-strain (λ=1)
Refined 80×16 mesh. Voltage ramp Ẽ = Φ/H = 10·clip((t−40)/1200,0,1).
- `stage1_l100.*`  static restart (λ=1, zero grip).
- `stage2_l100_g{005,010,020,050,100,200}.*`  γ̄ = 0.5,1,2,5,10,20.
- `gamma_onset.json`  critical nominal field vs γ̄ (Ẽc√(ε/μ) = 2.30/2.59/2.92/3.36/3.66/4.00).
- ramp-rate (dynamic-effects) test: `s1dyn`, `s2dyn_{fast,base,slow}` — onset 2.91/2.60/2.39.
- `exploratory/`  earlier young-Laplace DE crease/wrinkle decks (staggered/monolithic).
Feeds Figure `lam1_compare` (FE crease vs linear-theory wrinkle).
