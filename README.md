# tidyXCMS

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/github/actions/workflow/status/stanstrup/tidyXCMS/R-CMD-check.yaml?label=R-CMD-check)](https://github.com/stanstrup/tidyXCMS/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://img.shields.io/codecov/c/github/stanstrup/tidyXCMS/main)](https://codecov.io/gh/stanstrup/tidyXCMS)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![GitHub release](https://img.shields.io/github/release/stanstrup/tidyXCMS.svg)](https://github.com/stanstrup/tidyXCMS/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Tidy interface for XCMS metabolomics data processing

## Overview

tidyXCMS provides a modern, tidy interface for working with XCMS metabolomics data, making it easier to integrate XCMS workflows with the tidyverse ecosystem.

## Installation

```r
# Install from GitHub
# devtools::install_github("stanstrup/tidyXCMS")
```

## Usage

```r
library(tidyXCMS)
library(xcms)

# Load preprocessed XCMS data
xdata <- loadXcmsData("xmse")

# Create long-format peak table
peak_table <- tidy_peaklist(xdata)

# Result: one row per feature per sample
# Includes peak intensities, m/z, RT, and sample metadata
head(peak_table)
```

Optional annotations:

- **CAMERA**: Group features by retention time, correlation across samples or EIC similarity. Add isotope, adduct, and pseudo-spectra annotations
- **MsFeatures**: Group features by retention time and correlation across samples.

See the [vignette](https://stanstrup.github.io/tidyXCMS/articles/long-format-peaklist.html) for detailed examples.

## Development

This package uses:

- **Semantic versioning** via semantic-release
- **Automated documentation** via pkgdown
- **Continuous integration** via GitHub Actions

## License

MIT License
