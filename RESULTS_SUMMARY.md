# Final Results Summary

## Central conclusion
The final synthetic suite supports the **quotient-invariant organization** of the revised BCQ mathematics, but it also falsifies the idea that the scalar generation condition `chi = 3` is a complete identity invariant.

## Headline results
- **Symplectic quotient:** 0/20,000 failures for both the scalar pairing and the marked relational charge fingerprint.
- **Independent toroidal replication:** 0/20,000 oriented-intersection failures.
- **Exhaustive collision search:** 117,649 integer charge vectors searched; **16,807 distinct vectors** satisfy `chi = 3`.
- **Reconstructibility ablation:** the full six-component marked relational signature separates all 16,807 collision states. Removing the one coordinate already fixed by `chi = 3` is redundant; removing any of the other five yields 2,401 unique signatures, i.e. mean seven-fold ambiguity.
- **HRR common-twist test:** 0/20,000 positive failures; one-sided negative controls had 0/20,000 escapes.
- **Gauge quotient:** maximum Casimir/spectral numerical error `1.4921e-13` while the median raw matrix displacement was `3.73385`; 10,000/10,000 trace-contamination controls detected.
- **Hypercharge basis:** 0/3,000 exact functional-covariance failures; coefficient coordinates changed in 2,992/3,000 transformations.
- **Discriminant stress:** condition number grew from `3.73` to about `3.0e7`; median coefficient reconstruction error grew from `1.25e-7` to `0.804`, while median observable error stayed around `1.1e-7`.
- **Unit quotient:** 29,179 admissible states; every raw dimensionful state changed under unit rescaling while the maximum normalized error in dimensionless variables was `4.13e-16`.
- **Ablation:** random charge mutations were detected by `chi` at 79.78% (95% Wilson CI 79.22-80.33%) versus 100% by the marked relational fingerprint. In 5,000 deliberately `chi`-preserving adversarial mutations, `chi` detected 0% and the relational fingerprint 100%.

## Revised invariant hierarchy
### Primitive identity layer (`I2`)
- brane charge/topological class `[Gamma_a]`
- background pairing structure `[X, Omega_X]`
- gauge representation/conjugacy class `[R]`
- hypercharge functional class `[Y]`
- spin/chirality class `[s, chi_4]`

### Derived layer (`I2`, derived)
- Euler/intersection pairing `chi_ab`
- generation multiplicity `N_gen = |chi_ab|`
- electric charge `Q_em = T3 + Y`

### Moduli layer (`I1`)
- dimensionless `g_s`
- `ell = L/sqrt(alpha')`
- `mu^2 = alpha' M^2`

### Removable coordinate layer (`I0`)
- raw Chan-Paton matrix entries
- hypercharge coefficient vector in a chosen basis
- raw dimensionful `L`, `M`, `alpha'`, and tension values as numerical coordinates

### Discriminant layer (`I3`)
- `Delta_Y = det(A^T A) = 0` as an explicit identifiability boundary
- `mu^2 = 0`, `ell = 0`, and `chi = 0` retained only as candidate critical loci pending model-specific physical interpretation

## Evidence boundary
These results test algebraic invariance, redundancy, reconstructibility, and failure boundaries in synthetic mathematical testbeds. They do **not** establish a physically valid BCQ compactification, D1-to-4D Standard Model embedding, anomaly/tadpole consistency, QCD confinement, or quark compositeness.
