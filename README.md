# PaleoTimeseries-toolkit
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.placeholder-blue)](https://zenodo.org/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/)
A reference-grade, pedagogical, and fully reproducible toolkit for analyzing multi-proxy paleoclimate time series. 
This repository provides audited, production-quality Python workflows designed specifically for unevenly spaced geological time series. It bridges the gap between raw multi-proxy datasets, multivariate dimensionality reduction, and advanced unevenly-spaced spectral analysis.
---
## 🌟 Key Features
1. **Audited REDFIT Engine (`redfit38e`)**
   * A line-by-line, verified Python port of the classic Fortran `redfit38e.f90` algorithm (Schulz & Mudelsee, 2002).
   * Audited to guarantee mathematical parity, including exact replication of the LCG random number generator (`ran1`) and Box-Muller transformation (`gasdev`) for red-noise AR(1) surrogate generation.
   * Resolves numerical edge cases and adds out-of-the-box support for $\chi^2$-based confidence intervals ($80\%$, $90\%$, $95\%$, $99\%$).
2. **Weighted Wavelet Z-Transform (WWZ)**
   * Time-frequency analysis utilizing the WWZ method (Foster, 1996) optimized for unevenly spaced series via `pyleoclim`.
   * Integrates cross-method validation, overlaying significant REDFIT peaks directly on the WWZ scalogram.
3. **Multi-Proxy Dimensionality Reduction (PCA)**
   * Dedicated pipeline for handling multivariate proxy datasets (e.g., microfossil abundances, geochemical markers).
   * Implements robust exploratory data analysis, correlation profiling, and Principal Component Analysis (PCA) to extract independent climatic modes.
4. **Rigorous Reproducibility & Pedagogy**
   * Over 6,000 words of comprehensive markdown cells explaining the mathematical physics, parameter ranges, and potential geological pitfalls.
   * Auto-generated reproducibility logs capturing seeds, package versions, and configuration matrices.
---
## 📂 Repository Structure
```directory
paleo-timeseries-toolkit/
├── README.md                            # You are here
├── LICENSE                              # MIT License
