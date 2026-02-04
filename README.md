# drs-duration

**drs-duration** is a Python package for estimating **seismic duration** from the
harmonic structure of the **Displacement Response Spectrum (DRS)**.

The method is based on a physically grounded observation: high-resolution DRS
curves derived from impulsive and near-fault ground motions exhibit **smooth and
harmonic-like segments across the spectrum**, reflecting the free oscillatory
response of SDOF systems subjected to impulsive excitation.

Although harmonic patterns may appear along different portions of the DRS, this
software focuses the analysis on the **predominant spectral region**, located
around the global DRS peak. This region concentrates the largest spectral demand
and is therefore the most relevant for engineering interpretation and seismic
performance assessment.

Instead of relying on time-domain energy thresholds, **drs-duration** estimates a
**spectral-based duration parameter** by automatically identifying the
predominant harmonic window of the DRS and performing a nonlinear sinusoidal fit
within that region.

The entire procedure is **fully automatic, objective, and reproducible**, and can
be executed through a single command-line instruction.

---

## Scientific background (brief)

The proposed approach exploits the fact that the residual response of an SDOF
system subjected to impulsive excitation inherently contains terms of the form:

\[R(T) = A \left| \sin\left(\frac{\pi t_d}{T} + \varphi \right) \right|\]

where \( t_d \) represents a characteristic duration parameter and \( T \) is the
structural period. This mathematical structure explains the appearance of
harmonic-like geometries in displacement response spectra, particularly in
near-fault records dominated by strong velocity pulses.

To account for small vertical biases that may arise from numerical discretization
or spectral processing, an extended model including a constant offset is also
considered:

\[R(T) = A \left| \sin\left(\frac{\pi t_d}{T} + \varphi \right) \right| + C\]

The constant \( C \) has **no physical meaning** and acts solely as a numerical
corrector.
---

## Features

- Computation of **Displacement Response Spectrum (DRS)** using Newmark-β integration  
- Automatic detection of the **predominant harmonic window** around the DRS peak  
- Derivative-based criterion with robust fallback strategies  
- Systematic **α-scan** to select the optimal derivative threshold  
- Nonlinear **sinusoidal duration fitting**, with and without constant offset  
- Generation of **publication-ready figures** (high resolution, consistent style)  
- Fully automated **command-line interface (CLI)** for batch processing  
- Reproducible outputs (CSV tables and PNG figures)

---

## Installation

### Requirements

- Python ≥ 3.9

### Install from source (recommended for research)

Clone the repository and install in editable mode:

```bash
pip install -e .

## USAGE

Run the software from the command line:

drs-duration --input examples/San_Fernando_en_metros.txt --out outputs

## COMMAND-LINE OPTIONS

--input : Path to the ground-motion record (.txt)

--out : Output directory

--no-plots : Disable figure generation (tables only)

## INPUT FORMAT

The input file must be a plain text file with two columns:

1. Time (or time step) in seconds
2. Ground acceleration (e.g., in m/s²)

The file may include a header line. Units must be consistent.

Example:

t(s)    ag(m/s^2)
0.00    -0.00123
0.02     0.00345
0.04    -0.00210
...

## OUTPUT STRUCTURE

For each executed record, the software creates a dedicated subfolder inside the
specified output directory containing:

- drs_full.csv – Full displacement response spectrum

- drs_full.png – Full DRS figure

- drs_zoom.png – Zoomed DRS around the predominant spectral window

- alpha_scan.csv – α-scan results

- r2_vs_alpha.png – R² versus α plot

- sine_fit.png – Final sinusoidal fits within the detected window

- results.csv – Summary of final parameters


## INTERPRETATION OF RESULTS

The estimated duration parameter 𝑡_𝑑 represents a spectral-based measure
of effective seismic duration, derived from the harmonic structure of the DRS.
Unlike traditional time-domain duration metrics, this parameter is less
sensitive to weak energy tails, late secondary peaks, and noise contamination.

The constant-offset model is provided for numerical robustness only. The
parameter 𝐶 has no physical interpretation.

## REPRODUCIBILITY

All results are fully reproducible. Given the same input record and command-line
options, the software will generate identical outputs.

## LICENSE

This project is released under the MIT License.

## CITATION

If you use this software in academic work, please cite it using the information
provided in the CITATION.cff file included in this repository.