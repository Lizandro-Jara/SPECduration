# Changelog

## [1.2.0] – 2026-05-25

### Added
- Robust derivative-based spectral window detection using percentile stabilization
- Operational alpha strategy with:
  
  $$
  \alpha = 0.30
  $$

  as the default configuration
- Restricted adaptive alpha search:
  
  $$
  0.20 \leq \alpha \leq 0.50
  $$

  for complex spectral cases
- Automatic physical validation of detected windows and fitted solutions
- New output metadata:
  - `selected_alpha`
  - `alpha_strategy`
  - `window_method`
  - `valid_solution`

### Changed
- Replaced the original derivative threshold based on the maximum derivative by a robust percentile-based formulation:

  $$
  \left| \frac{dR}{dT} \right|
  \leq
  \alpha \cdot P95\left(\left| \frac{dR}{dT} \right|\right)
  $$

- Removed the global alpha optimization strategy based on maximum $R^2$
- Replaced unrestricted alpha scans with operational and restricted-search strategies
- Formalized the envelope fallback as an official stabilization stage for complex DRS geometries
- Improved consistency and physical coherence of automatically detected spectral windows
- Increased robustness against:
  - micro-window overfitting
  - isolated derivative spikes
  - unstable local spectral regions

### Improved
- More stable duration estimation across near-fault and impulsive records
- Better handling of multiscale and geometrically complex DRS cases
- Reduced sensitivity to local spectral irregularities
- Improved reproducibility and methodological traceability
- Cleaner and more interpretable automatic outputs

---

## [1.1.0] – 2026-04-14

### Changed
- Removed the constant-offset formulation from the official workflow
- Simplified sinusoidal fitting to the physically interpretable model without constant offset
- Improved consistency between software outputs, README, and manuscript

### Improved
- Cleaner and more consistent output files for reproducible research
- Clearer traceability in final results, including window information and fitting diagnostics

---

## [1.0.0] – 2026-02-05

### Added
- Initial release of drs-duration
- Command-line interface (CLI) for spectral estimation of seismic duration
- High-resolution Displacement Response Spectrum (DRS) computation using the Newmark–β method
- Automatic detection of the predominant spectral window
- Alpha-scan procedure for derivative threshold exploration
- Nonlinear sinusoidal fitting
- Generation of reproducible CSV tables and diagnostic figures