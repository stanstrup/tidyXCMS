# tidyXCMS

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
peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata)

# Result: one row per feature per sample
# Includes peak intensities, m/z, RT, and sample metadata
head(peak_table)
```

Optional annotations:
- **CAMERA**: Add isotope, adduct, and pseudospectrum annotations
- **MsFeatures**: Group features by retention time, correlation, or EIC similarity

See the [vignette](https://stanstrup.github.io/tidyXCMS/articles/long-format-peaklist.html) for detailed examples.

## Development

This package uses:

- **Semantic versioning** via semantic-release
- **Automated documentation** via pkgdown
- **Continuous integration** via GitHub Actions

## License

MIT License
