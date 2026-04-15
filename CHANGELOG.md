# Changelog

## [1.1.0] – 2026-04-14
### Changed
- Removed the constant-offset formulation from the official workflow
- Simplified sinusoidal fitting to the physically interpretable model without constant offset
- Improved consistency between software outputs, README, and manuscript

### Improved
- Cleaner and more consistent output files for reproducible research
- Clearer traceability in final results, including window information and fitting diagnostics

## [1.0.0] – 2026-02-05
### Added
- Initial release of **drs-duration**
- Command-line interface (CLI) for spectral estimation of seismic duration
- High-resolution Displacement Response Spectrum (DRS) computation using the Newmark–β method
- Automatic detection of the predominant spectral window
- α-scan procedure for optimal derivative threshold selection
- Nonlinear sinusoidal fitting
- Generation of reproducible CSV tables and diagnostic figures