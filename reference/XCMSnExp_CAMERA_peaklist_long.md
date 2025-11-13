# Create Long-Format Peak Table with Optional CAMERA Annotations

Generates a comprehensive long-format peak table where each row
represents one feature for one sample. Optionally includes CAMERA
annotations for isotopes, adducts, and pseudospectrum groups if
provided.

## Usage

``` r
XCMSnExp_CAMERA_peaklist_long(XCMSnExp, xsAnnotate = NULL)
```

## Arguments

- XCMSnExp:

  An [xcms::XCMSnExp](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
  or
  [xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html)
  object containing peak detection and feature grouping results.

- xsAnnotate:

  Optional. An `xsAnnotate` object from the CAMERA package containing
  peak annotations (isotopes, adducts, pseudospectrum groups). If NULL
  (default), annotation columns will not be included in the output.

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
in long format with one row per feature per sample. The tibble contains:

- Feature-level columns:

  f_mzmed, f_mzmin, f_mzmax, f_rtmed, f_rtmin, f_rtmax - feature summary
  statistics

- CAMERA annotations:

  isotopes, adduct, pcgroup (pseudospectrum group) - only present when
  xsAnnotate is provided

- MsFeatures grouping:

  feature_group - feature group ID from
  [`MsFeatures::groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html),
  only present if groupFeatures was applied

- Peak-level columns:

  mz, rt, into, intb, maxo, sn - individual peak measurements

- Sample information:

  filepath, filename, fromFile - sample identifiers

- Additional columns:

  Any columns from pData(XCMSnExp) or sampleData(XcmsExperiment)

## Details

This function integrates data from multiple sources:

- Chromatographic peaks from
  [`xcms::chromPeaks()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)

- Feature definitions from
  [`xcms::featureDefinitions()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)

- Feature intensities from
  [`xcms::featureValues()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-peak-grouping-results.html)

- CAMERA annotations from the xsAnnotate object (if provided)

- Sample metadata from
  [`Biobase::pData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html)

The function supports both
[xcms::XCMSnExp](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html) and
[xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html)
objects, providing flexibility for modern XCMS workflows.

The function performs several data processing steps:

1.  Extracts peaks and adds file path information

2.  Combines feature definitions with CAMERA annotations

3.  Unnests peak indices to create feature-peak relationships

4.  Joins feature intensities to match peaks with features

5.  Filters to ensure peak intensities match feature intensities

6.  Handles duplicate peaks by taking the first occurrence

7.  Completes the data frame to include all feature-sample combinations

8.  Adds sample metadata from pData

Missing values (NA) indicate that a feature was not detected in a
particular sample.

## Note

- XCMS sometimes creates duplicate peaks; this function keeps only the
  first occurrence.

- The function uses "maxint" method for feature values, meaning it takes
  the maximum intensity peak for each feature in each sample.

- CAMERA annotations are optional. If xsAnnotate is NULL, the isotopes,
  adduct, and pcgroup columns will not be present in the output.

- MsFeatures grouping is automatically detected. If
  [`MsFeatures::groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)
  was applied to the object, the feature_group column will be included
  in the output.

- For XcmsExperiment objects, sample metadata is accessed via
  sampleData() instead of pData().

## Examples

``` r
if (FALSE) { # \dontrun{
library(xcms)
library(BiocParallel)

# Load example data
faahko_sub <- loadXcmsData("faahko_sub")

# Peak detection
cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())

# Peak grouping
pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                        minFraction = 0.5)
xdata <- groupChromPeaks(xdata, param = pdp)

# Example 1: Without CAMERA annotations
peak_table_no_camera <- XCMSnExp_CAMERA_peaklist_long(xdata)
head(peak_table_no_camera)

# Example 2: With CAMERA annotations
library(CAMERA)
library(commonMZ)

# Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
xset <- as(xdata, "xcmsSet")

# CAMERA annotation
xs <- xsAnnotate(xset, polarity = "positive")
xs <- groupFWHM(xs, perfwhm = 0.1, intval = "into", sigma = 6)
xs <- findIsotopes(xs, ppm = 10, mzabs = 0.01, intval = "into")
xs <- groupCorr(xs, calcIso = FALSE, calcCiS = TRUE, calcCaS = FALSE,
                cor_eic_th = 0.7, pval = 1E-6)

# Get adduct/fragment rules from commonMZ
rules_pos <- commonMZ::MZ_CAMERA(mode = "pos", warn_clash = TRUE, clash_ppm = 5)
rules_pos <- as.data.frame(rules_pos)

# Find adducts using the rules
xs <- findAdducts(xs, ppm = 10, mzabs = 0.01, multiplier = 4,
                  polarity = "positive", rules = rules_pos)

# Create long-format peak table (pass XCMSnExp object, not xcmsSet)
peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata, xs)
head(peak_table)

# Explore the data
library(dplyr)

# Count features per sample
peak_table %>%
  group_by(filename) %>%
  summarize(n_features = sum(!is.na(into)))

# View adduct annotations
peak_table %>%
  filter(!is.na(adduct)) %>%
  select(feature_id, f_mzmed, f_rtmed, adduct, pcgroup)

# Example 3: With XcmsExperiment object
library(MsExperiment)
library(msdata)

# Load data as XcmsExperiment
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
xdata_exp <- readMsExperiment(spectraFiles = fls[1:2], BPPARAM = SerialParam())

# Peak detection and grouping
cwp <- CentWaveParam(peakwidth = c(5, 20), noise = 100)
xdata_exp <- findChromPeaks(xdata_exp, param = cwp, BPPARAM = SerialParam())
pdp <- PeakDensityParam(sampleGroups = rep(1, 2), minFraction = 0.5)
xdata_exp <- groupChromPeaks(xdata_exp, param = pdp)

# Create peak table without CAMERA
peak_table_exp <- XCMSnExp_CAMERA_peaklist_long(xdata_exp)
head(peak_table_exp)
} # }
```
