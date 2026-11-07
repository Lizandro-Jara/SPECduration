
# SPECduration

SPECduration is an open-source Python package for estimating seismic duration from the displacement response spectrum (DRS) of processed acceleration records.

The software automatically identifies the predominant sinusoidal peak of the DRS and estimates the seismic duration $T_s$ by means of harmonic fitting:

$$
DRS(T)=P\sin\left(\frac{\pi T_s}{T}-c\right)
$$

where $P$ is the fixed maximum DRS value, $T$ is the natural period, and $c$ is the phase constant.

A full and clear explanation of SPECduration (in fewer than 3,500 words) is provided in the accompanying SoftwareX manuscript. The theoretical basis of the spectral-duration method is described in the [original research paper](https://doi.org/10.1177/10775463251360873).

## Features

- Spectral seismic-duration estimation in the natural-period domain using the predominant sinusoidal shape of the DRS.
- Automatic determination of the predominant fitting window using the second spectral derivative, replacing manual and visual selection.
- Nonlinear harmonic fitting with fixed maximum DRS amplitude to estimate the seismic duration $T_s$.
- Evaluation of fitting quality through the coefficient of determination $R^2$.
- Automatic analysis from a processed acceleration record, without requiring user-defined duration thresholds or manual window selection.

## Installation

Clone the repository and install the package:

```bash
git clone https://github.com/Lizandro-Jara/SPECduration.git
cd SPECduration
python -m pip install .
```

Alternatively, download the repository as a ZIP file, extract it, and run the same installation command from the project directory.

A virtual environment is recommended.

## Usage

Run SPECduration from the command line using a processed earthquake acceleration record:

```bash
specduration --input "examples/LOMA PRIETA.txt" --out "outputs"
```

Alternatively, execute the Python module directly:

```bash
python -m specduration.cli --input "examples/LOMA PRIETA.txt" --out "outputs"
```

The input must be a two-column `.txt` file containing time in seconds and ground acceleration in m/s², preceded by a header line.

Example:

```text
time(s)  acceleration(m/s^2)
0.00     -0.00123
0.02      0.00345
0.04     -0.00210
```

The default analysis uses a damping ratio of 5% and a spectral period step of 0.01 s for DRS computation.

For additional command-line options, run:

```bash
specduration --help
```

## Outputs

SPECduration generates the following files for each processed record:

**Main report**

- `results.xlsx` — Excel report containing the estimated seismic duration $T_s$ and coefficient of determination $R^2$.

**CSV files**

- `summary.csv`
- `results_diagnostics.csv`
- `drs_full.csv`
- `second_derivative.csv`

**Diagnostic figures**

- `drs_full.png`
- `drs_window.png`
- `second_derivative_window.png`

### Fig. 2. Displacement response spectrum for Loma Prieta at Corralitos

![Displacement response spectrum for Loma Prieta at Corralitos](docs/figures/drs_full.png)

### Fig. 3. Predominant DRS peak (blue) and harmonic fit (orange) over computed window for Corralitos record

![Predominant DRS peak and harmonic fit](docs/figures/drs_window.png)

## Citation and License

If you use SPECduration in academic work or elsewhere, please cite the software using the information provided in [`CITATION.cff`](CITATION.cff).

SPECduration is distributed under the [MIT License](LICENSE).