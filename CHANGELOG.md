# Changelog

## [1.2.0] - 2026-05-25

### Added

- Final spectral window detection based on the second derivative of the displacement response spectrum, S''(T).
- Automatic detection of the predominant spectral window around the peak period, T_peak.
- Harmonic fitting within the detected window limits, T_left and T_right.
- Default spectral resolution set to dT = 0.01 s.
- Diagnostic figure for the second-derivative window detection:
  - `second_derivative_window.png`
- Main Excel output for user-friendly result visualization:
  - `results.xlsx`
- Additional CSV outputs for reproducibility:
  - `summary.csv`
  - `results_diagnostics.csv`
  - `drs_full.csv`
  - `second_derivative.csv`

### Changed

- Renamed the software from `drs-duration` to `SPECduration`.
- Updated the command-line interface from `drs-duration` to `specduration`.
- Replaced the previous alpha-based window detection workflow with the final second-derivative curvature-based detector.
- Replaced the previous notation R(T) with S(T) for the displacement response spectrum.
- Updated the harmonic fitting workflow to operate only inside the automatically detected spectral window.
- Updated the output structure to separate the main user results from diagnostic and reproducibility files.
- Updated the README, citation metadata, package metadata, and manuscript terminology to match the final SPECduration workflow.

### Removed

- Operational alpha strategy from the official workflow.
- Restricted adaptive alpha search from the official workflow.
- Envelope-based fallback stabilization from the official workflow.
- Alpha-scan output files:
  - `alpha_scan.csv`
  - `r2_vs_alpha.png`
- Alpha-related result fields from the official output structure:
  - `selected_alpha`
  - `alpha_strategy`

### Improved

- Clearer and more reproducible spectral window detection.
- Reduced dependence on empirical threshold parameters.
- Cleaner output files for users and researchers.
- Improved separation between main results and diagnostic information.
- More consistent terminology between the software, README, citation file, and manuscript.

---

## [1.1.0] - 2026-04-14

### Added

- Derivative-based spectral window detection using percentile stabilization.
- Operational alpha strategy with alpha = 0.30 as the default configuration.
- Restricted adaptive alpha search for complex spectral cases.
- Auxiliary envelope-based fallback stabilization.
- Automatic validation of detected windows and fitted solutions.
- Additional metadata for alpha-based diagnostics.

### Changed

- Removed the constant-offset formulation from the official harmonic fitting model.
- Simplified the harmonic fitting model to the formulation without constant offset.
- Replaced unrestricted alpha scans with operational and restricted-search strategies.
- Improved consistency between software outputs, README, and manuscript.

### Improved

- More stable duration estimation across near-fault and impulsive records.
- Better handling of multiscale and geometrically complex DRS cases.
- Reduced sensitivity to local spectral irregularities.
- Clearer traceability of fitting diagnostics.

---

## [1.0.0] - 2026-02-05

### Added

- Initial release of `drs-duration`.
- Command-line interface for spectral estimation of seismic duration.
- High-resolution displacement response spectrum computation using the Newmark-beta method.
- Automatic detection of the predominant spectral window.
- Alpha-scan procedure for derivative threshold exploration.
- Nonlinear harmonic fitting.
- Generation of CSV tables and diagnostic figures.