# FRED-MD Data Cleaning Pipeline

An R pipeline for cleaning and processing the [FRED-MD](https://research.stlouisfed.org/econ/mccracken/fred-databases/) macroeconomic database. Implements McCracken & Ng (2016, JBES) transformation codes and Bai & Ng (2002, ECTA) factor extraction via a compiled C++ (RcppArmadillo) backend.

## Overview

FRED-MD is a large monthly macroeconomic database maintained by the Federal Reserve Bank of St. Louis. This pipeline reads the raw CSV, applies the 7-code transformation scheme to each series, removes outliers, and extracts common factors using an EM algorithm with PCA. The output is a balanced panel (`fred_md_balanced`) ready for downstream modeling.

## Requirements

- R ≥ 4.0
- A C++ compiler (Xcode Command Line Tools on macOS, `build-essential` on Linux)

Install R packages once:

```r
install.packages(c("tidyverse", "lubridate", "janitor", "Rcpp", "RcppArmadillo", "testthat"))
```

## Usage

Place the raw FRED-MD CSV in `Data/` (e.g. `Data/FRED_MD_062024.csv`), then run:

```bash
Rscript load.R
```

`src/fred_factors.cpp` is compiled automatically on the first run (~10–30 s). Subsequent runs reuse the cached shared library.

### Parameters

All tunable parameters are at the top of `load.R`:

| Variable | Default | Description |
|----------|---------|-------------|
| `demean` | `2` | Normalization: `1` = demean only, `2` = demean + standardize, `3` = recursive, `0` = none |
| `jj` | `2` | Information criterion: `1` = PC\_p1, `2` = PC\_p2, `3` = PC\_p3 |
| `kmax` | `8` | Max factors to consider; set to `99` to force exactly 8 factors |

### Output

The pipeline produces `fred_md_balanced`: a tibble with a `date` column followed by all imputed series (snake\_case names), with any remaining `NA` rows dropped.

## Running Tests

```bash
Rscript -e "testthat::test_file('tests/test_fred_factors.R')"
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 23 ]`

## Architecture

```
load.R                      # Orchestrator
R/prepare-missing.R         # Applies 7-code transformations per series
R/remove-outliers.R         # Flags |x − median| > 10×IQR as NA
src/fred_factors.cpp        # Compiled C++ (RcppArmadillo): EM + PCA kernels
tests/test_fred_factors.R   # C++ vs. R reference regression tests
Data/FRED_MD_062024.csv     # Raw input (127 monthly series, 1959–2024)
```

### Pipeline Steps

1. **`prepare_missing()`** — reads transformation codes from row 1 of the CSV, applies `transxf()` to each series (levels, differences, log-differences, etc.)
2. **`remove_outliers()`** — replaces values more than 10×IQR from the median with `NA`
3. **`factors_em_cpp()`** — compiled EM algorithm: imputes missing values, demeans/standardizes, runs PCA, iterates until convergence (tolerance < 1e-6)

### C++ Kernels (`src/fred_factors.cpp`)

| Function | Description |
|----------|-------------|
| `factors_em_cpp(x, kmax, jj, DEMEAN)` | Main entry point; hot T×N inner loop; `kmax=99` forces 8 factors |
| `baing_cpp(X, kmax, jj)` | Bai & Ng (2002) information-criterion factor count selection |
| `pc2_cpp(X, nfac)` | PCA via economy SVD of the gram matrix X'X |
| `transform_data_cpp(x2, DEMEAN)` | Normalization; returns T0×N `mut`/`std` matrices |
| `minindc_cpp(X)` | Column-wise argmin (1-based row index) |

The `factors_em_cpp` function expects a **pre-stripped numeric matrix** (no date column). `load.R` strips the date before the call and re-attaches it to the output tibble.

## References

- McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research. *Journal of Business & Economic Statistics*, 34(4), 574–589.
- Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. *Econometrica*, 70(1), 191–221.

