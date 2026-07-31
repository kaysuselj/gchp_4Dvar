# AGU Fall Meeting 2026 — Poster Abstract (DRAFT)

**Topic:** GCHP cubed-sphere transport adjoint for CO₂ flux inversion, emphasizing
adjoint validation across resolutions up to C180.

**Status:** Draft — 2026-07-31. Body ≈ 1,760 characters (AGU limit 2,000). ~234 words.

---

## Title (options)

- **A.** Building and Validating a Discrete Adjoint of the GCHP Cubed-Sphere
  Transport Model for CO₂ Flux Inversion Across Resolutions from C24 to C180
- **B.** Resolution-Dependent Validation of a GEOS-Chem High Performance Adjoint
  for 4D-Var Surface CO₂ Flux Estimation

## Authors

K. Suselj, [co-authors], [affiliations]

## Suggested section

Atmospheric Sciences (A) or Global Environmental Change (GC) — carbon-cycle /
biogeosciences inverse-modeling sessions.

---

## Abstract body

Top-down estimation of surface CO₂ fluxes from satellite column observations
(e.g., OCO-2) requires the adjoint of an atmospheric transport model to compute
the gradient of the misfit with respect to millions of flux parameters. We
present the development and validation of a discrete adjoint of the GEOS-Chem
High Performance (GCHP) model, which uses the FV3 finite-volume cubed-sphere
dynamical core, and its application to variational (4D-Var) CO₂ surface-flux
inversion.

A central challenge is guaranteeing adjoint correctness as horizontal resolution
increases: the transport adjoint must remain the exact transpose of the
tangent-linear advection operator on the cubed sphere, including at panel edges
and in the presence of flux limiters and spatially varying winds. We validate
the adjoint using offline adjoint-identity (dot-product) tests and
finite-difference gradient checks, applied systematically across resolutions
from C24 to C180. These tests isolate errors to individual operators and reveal
discretization-dependent inaccuracies in the reversed-wind transport
adjoint—exact for uniform flow but sensitive to wind gradients—which we diagnose
and correct.

We embed the validated adjoint in an Observing System Simulation Experiment
(OSSE) in which fixed synthetic OCO-2 retrievals, generated from ORCHIDEE-ECCO2
fluxes, constrain a monthly terrestrial-respiration flux control vector; adjoint
memory at high resolution is managed by checkpointing. We report the accuracy
and computational scaling of the adjoint and
gradient at each resolution, and discuss implications for high-resolution
regional-to-global carbon-flux inversions. This work establishes a verified,
resolution-scalable adjoint foundation for GCHP-based CO₂ data assimilation.

---

## Notes / to finalize before submission

- **Resolution list:** currently "C24 to C180." If C48/C90 intermediate points
  will be shown, name them for concreteness.
- **Wind-gradient adjoint bug:** framed here as "diagnose and correct." If the
  correction is not in hand by December, soften to "diagnose and characterize."
- **Placeholders:** co-authors, affiliations, AGU section, plain-language
  significance statement (AGU requests one separately).
- **~280 characters of headroom** remain under the 2,000-char limit for one more
  concrete sentence (e.g., a preliminary result or the intermediate resolutions).
