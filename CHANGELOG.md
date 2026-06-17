# Changelog

## [1.2.0] – 2026-05-25

### Added

- Final curvature-based spectral window detection using the second derivative of the displacement response spectrum.
- Automatic identification of the dominant negative-curvature region around the peak period \(T_p\).
- Detection of the fitting limits \(T_l\) and \(T_r\) from the internal boundaries of the \(S''(T)<0\) region.
- Default spectral resolution:

  $$
  \Delta T = 0.01 \ \text{s}
  $$

- New diagnostic output for second-derivative-based window detection:
  - `second_derivative_window.png`

- Updated output structure:
  - `drs_full.csv`
  - `second_derivative.csv`
  - `results.csv`
  - `results_diagnostics.csv`
  - `drs_full.png`
  - `drs_window.png`
  - `second_derivative_window.png`

- Updated result metadata:
  - `Tp`
  - `Tl`
  - `Tr`
  - `n`
  - `td`
  - `R2`
  - `status`

### Changed

- Renamed the software from `drs-duration` to `SPECduration`.
- Updated the command-line interface from:

  ```bash
  drs-duration --input <file> --out <folder>
  ```

  to:

  ```bash
  specduration --input <file> --out <folder>
  ```

- Replaced the previous alpha-based window detection workflow with a second-derivative curvature-based detector.
- Replaced the previous notation \(R(T)\) with \(S(T)\) for the displacement response spectrum.
- Updated the harmonic fitting workflow to operate within the automatically detected curvature-based window \([T_l,T_r]\).
- Updated the README, citation metadata, package metadata, and manuscript terminology to match the final SPECduration workflow.
- Adopted \(\Delta T=0.01\) s as the default operational spectral resolution.

### Removed

- Operational alpha strategy.
- Restricted adaptive alpha search.
- Envelope-based fallback stabilization.
- Alpha-scan outputs:
  - `alpha_scan.csv`
  - `r2_vs_alpha.png`
- Alpha-related metadata:
  - `selected_alpha`
  - `alpha_strategy`
  - `window_method`

### Improved

- More consistent detection of the dominant DRS lobe.
- Better preservation of local spectral geometry for window detection.
- Reduced dependence on empirical threshold parameters.
- Improved methodological clarity and reproducibility.
- Cleaner and more traceable diagnostic outputs.

---

## [1.1.0] – 2026-04-14

### Added

- Robust derivative-based spectral window detection using percentile stabilization.
- Operational alpha strategy with:

  $$
  \alpha = 0.30
  $$

  as the default configuration.

- Restricted adaptive alpha search:

  $$
  0.20 \leq \alpha \leq 0.50
  $$

  for complex spectral cases.

- Auxiliary envelope-based fallback stabilization for complex spectral geometries.
- Automatic physical validation of detected windows and fitted solutions.
- New output metadata:
  - `selected_alpha`
  - `alpha_strategy`
  - `window_method`
  - `valid_solution`

### Changed

- Combined the previous 1.1.0 and intermediate development updates into a single historical release entry.
- Removed the constant-offset formulation from the official workflow.
- Simplified the harmonic fitting model to the formulation without constant offset.
- Replaced the original derivative threshold based on the maximum derivative by a robust percentile-based formulation:

  $$
  \left| \frac{dR}{dT} \right|
  \leq
  \alpha \cdot P95\left(\left| \frac{dR}{dT} \right|\right)
  $$

- Removed the global alpha optimization strategy based on maximum \(R^2\).
- Replaced unrestricted alpha scans with operational and restricted-search strategies.
- Improved consistency between software outputs, README, and manuscript.

### Improved

- More stable duration estimation across near-fault and impulsive records.
- Better handling of multiscale and geometrically complex DRS cases.
- Reduced sensitivity to local spectral irregularities.
- Cleaner output files for reproducible research.
- Clearer traceability in final results, including window information and fitting diagnostics.

---

## [1.0.0] – 2026-02-05

### Added

- Initial release of `drs-duration`.
- Command-line interface for spectral estimation of seismic duration.
- High-resolution Displacement Response Spectrum computation using the Newmark–β method.
- Automatic detection of the predominant spectral window.
- Alpha-scan procedure for derivative threshold exploration.
- Nonlinear harmonic fitting.
- Generation of reproducible CSV tables and diagnostic figures.