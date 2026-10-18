# Threshold-consistent response

`threshold-kk-moment-budget-v2` repairs the dielectric construction, not the
projectile theory. It is a revised constrained optical-data model, not an
unchanged reproduction of ECN or an experimental validation. No ICRU data,
screening fit, or adjustable stopping-power multiplier enters this repair.

## Kramers--Kronig pair

Redistribution conserves the total imaginary valence spectrum above the common
onset g; below g it is zero. Retaining the untruncated real response left a
mismatched pair. We now subtract the KK transform of the removed interval:

    epsilon1_cut(E)-1 = epsilon1_raw(E)-1
        - (2/pi) PV integral_0^g w*epsilon2_raw(w)/(w^2-E^2) dw.

The finite subtraction is evaluated through rational integrals and logarithms,
with regularized Gaussian quadrature near repeated poles. For derivative Drude
bands, the susceptibility is `(1-gamma*d/dgamma)/D`, with
`D=E0^2-E^2-i*gamma*E`. The real-band decomposition precedes redistribution;
only its total is the KK partner of the total partitioned imaginary response.

The retained hard onset has a logarithmic real-part singularity. The exact
edge ELF limit is zero; the real value uses the nearest upper floating-point
energy. No smoothing parameter is introduced. Smoother physical onsets would
need constraints from optical data.

## Actual moment allocation

Let Rj(q),Ri(q) be the fractions of the unit excitation and outer-ionization
energy-weighted moments above g. Retain the prescribed excitation dispersion
and relative outer amplitudes, determining their common multiplier by

    A(q) = [1-FK(q)-sum_j fj(q)*Rj(q)] / sum_i fi(0)*Ri(q).

The retained valence conductivity moment is `1-FK(q)`. Its KK-consistent
inverse dielectric has the same high-frequency moment. Adding the specified
core ELF therefore closes the molecular f-sum. This replaces bookkeeping of
untruncated amplitudes; no final ELF or core spectrum is rescaled. The former
logarithmic post-inversion tail multiplier is removed because it invalidated
the preceding constraints. No eight-plus-two occupancy allocation is imposed.

Ions use the unscaled, q-dependent Heredia--Avalos continuum moment, including
for water. Electrons retain their fixed optical Drude core. The additive core
approximation `ELF_K ~= epsilon2_K` remains; the valence response in the RPWBA
density term is not an exact coupled core/valence dielectric. Electron core
dispersion, low-energy Born validity and projectile screening remain separate
model limitations. Moment closure does not establish experimental accuracy.

## Checks and regeneration

`tests/test_causal_dielectric.py` checks independent principal-value integrals,
full valence KK transforms, numerical conductivity and ELF moments at
q=0,0.5,1,2.5,8,30 a0^-1, threshold limits and scalar/batched/optical parity.
The corresponding ion-repository tests check ion kinematics and normalization;
this repository also tests an electron DAT round trip and cache rejection.

All old Born and dependent oscillator/polarization products require regeneration.
This commit does not regenerate production DATs or simulations. Retain accepted
old results as versioned comparison evidence; do not append new-response rows
to them. Recheck optical data, stopping and W-values before replacing transport
inputs. Water in the ion repository remains diagnostic-only pending that qualification.

The constraints are literature-backed; this particular repaired combination is
not attributed to a published model unchanged:

- Dingfelder (2014), Eqs.9--26 and Appendix C:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3692571/ .
- Emfietzoglou et al. (2009), sum rules and momentum dispersion:
  https://doi.org/10.1088/0031-9155/54/11/012 .
- Kyriakou et al. (2017), redistribution:
  https://doi.org/10.1667/RR14705.1 .
- Heredia--Avalos et al. (2005), retained ion core:
  https://doi.org/10.1103/PhysRevA.72.052902 .
