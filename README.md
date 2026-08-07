# PaleoTimeseries-toolkit
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.11+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/arjun-vh/PaleoTimeseries-toolkit/HEAD?urlpath=lab)

A reference-grade, pedagogical, and fully reproducible toolkit for analyzing multi-proxy paleoclimate time series. 
This repository provides quality Python workflows designed specifically for unevenly spaced geological time series. It bridges the gap between raw multi-proxy datasets, multivariate dimensionality reduction, and advanced unevenly-spaced spectral analysis.

---
## 🌟 Key Features
1. **Python REDFIT Engine (`redfit38e`)**
   * A line-by-line Python port of the classic Fortran `redfit38e.f90` algorithm (Schulz & Mudelsee, 2002).
   * Complete mathematical parity, including exact replication of the LCG random number generator (`ran1`) and Box-Muller transformation (`gasdev`) for red-noise AR(1) surrogate generation.
   * Resolves numerical edge cases and adds support for $\chi^2$-based confidence intervals ($80\%$, $90\%$, $95\%$, $99\%$).
2. **Weighted Wavelet Z-Transform (WWZ)**
   * Time-frequency analysis utilizing the WWZ method (Foster, 1996) optimized for unevenly spaced series via `pyleoclim`.
   * Integrates cross-method validation, overlaying significant REDFIT peaks directly on the WWZ scalogram.
3. **Multi-Proxy Dimensionality Reduction (PCA)**
   * Pipeline for handling multivariate proxy datasets (e.g., microfossil abundances, geochemical markers).
   * Implements exploratory data analysis, correlation profiling, and Principal Component Analysis (PCA) to extract independent climatic modes.
4. **Reproducibility & Pedagogy**
   * Comprehensive markdown cells explaining the mathematical physics, parameter ranges, and potential geological pitfalls.
   * Reproducibility logs capturing seeds, package versions, and configuration matrices.
---
## 📂 Repository Structure
```directory
paleo-timeseries-toolkit/
├── README.md                            # You are here
├── LICENSE                              # MIT License
├── environment.yml                      # Conda environment definition
├── requirements.txt                     # Pip dependencies
│
├── notebooks/
│   ├── REDFIT_WWZ_Reference.ipynb       # Spectral & Wavelet Analysis Notebook
│   └── TimeSeries_PCA_Reference.ipynb   # Multi-Proxy PCA Notebook (Companion)
│
├── redfit_engine/
│   ├── __init__.py
│   └── redfit38e_python.py              # Standalone, importable REDFIT engine
│
├── data/
│   ├── README.md                        # Data documentation & licenses
│   └── example_15ka.xlsx                # Sample multi-proxy dataset
│
└── tests/
    ├── test_engine_fortran_parity.py    # Numerical unit tests matching Fortran output
    └── test_reproducibility.py          # Seed verification tests
```
---
## 🚀 Quick Start
> **No installation at all?** Click the Binder badge at the top to run the notebooks live in your browser — no local setup required.
---
## 🛠️ Installation Guide
### Prerequisites
| Requirement | Minimum Version | Notes |
|---|---|---|
| Python | 3.9 | 3.11 recommended |
| Conda / Mamba | Any recent | Strongly preferred over bare pip |
| Git | Any | For cloning the repository |
| Jupyter | 6.x or 7.x | Included in `environment.yml` |
* Fitted autocorrelation ($\rho$) and characteristic relaxation time ($\tau$).
* Output frequencies, spectral amplitudes, and theoretical $\chi^2$ significance thresholds.
---
## 📝 Citation
If you use this toolkit in your research, please cite the repository and the foundational papers:
```bibtex
@software{paleo_timeseries_toolkit,
  author       = {Arjun, V. H.},
  title        = {paleo-timeseries-toolkit: Reference-grade Python workflows for paleoclimate time series},
  month        = aug,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.placeholder},
  url          = {https://github.com/yourusername/paleo-timeseries-toolkit}
}
@article{schulz2002redfit,
  title={REDFIT: estimating red-noise spectra of unevenly spaced paleoclimatic time series},
  author={Schulz, Michael and Mudelsee, Manfred},
  journal={Computers \& Geosciences},
  volume={28},
  number={3},
  pages={421--426},
  year={2002},
  publisher={Elsevier}
}
@article{foster1996wavelets,
  title={Wavelets for period analysis of unevenly sampled time series},
  author={Foster, Grant},
  journal={The Astronomical Journal},
  volume={112},
  pages={1709},
  year={1996}
}
```
---
## 📄 License
This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
