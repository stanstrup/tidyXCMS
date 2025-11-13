# Create Tidy Long-Format Peak Table with Optional Annotations

Generates a comprehensive long-format peak table where each row
represents one feature for one sample. Optionally includes CAMERA
annotations for isotopes, adducts, and pseudospectrum groups if
provided.

## Usage

``` r
tidy_peaklist(x, xsAnnotate = NULL)
```

## Arguments

- x:

  An [xcms::XCMSnExp](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
  or
  [xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html)
  object containing peak detection results. Feature grouping (via
  [`xcms::groupChromPeaks()`](https://rdrr.io/pkg/xcms/man/groupChromPeaks.html))
  is optional but recommended.

- xsAnnotate:

  Optional. An `xsAnnotate` object from the CAMERA package containing
  peak annotations (isotopes, adducts, pseudospectrum groups). If NULL
  (default), annotation columns will not be included in the output.

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
in long format with one row per feature per sample (or one row per peak
per sample if features are not defined). The tibble contains:

- Feature-level columns:

  f_mzmed, f_mzmin, f_mzmax, f_rtmed, f_rtmin, f_rtmax - feature summary
  statistics (only if features are defined)

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

  Any columns from pData(x) or sampleData(x)

## Details

This function integrates data from multiple sources:

- Chromatographic peaks from
  [`xcms::chromPeaks()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)

- Feature definitions from
  [`xcms::featureDefinitions()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)
  (if available)

- Feature intensities from
  [`xcms::featureValues()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-peak-grouping-results.html)
  (if features defined)

- CAMERA annotations from the xsAnnotate object (if provided)

- Sample metadata from
  [`Biobase::pData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html)
  or
  [`MsExperiment::sampleData()`](https://rdrr.io/pkg/MsExperiment/man/MsExperiment.html)

The function supports both
[xcms::XCMSnExp](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html) and
[xcms::XcmsExperiment](https://rdrr.io/pkg/xcms/man/XcmsExperiment.html)
objects, providing flexibility for modern XCMS workflows.

**Feature grouping is optional**: If
[`xcms::groupChromPeaks()`](https://rdrr.io/pkg/xcms/man/groupChromPeaks.html)
has not been run, the function will return peak-level data only, without
feature-level aggregation.

**Intensity matching**: When multiple peaks from the same sample are
assigned to a single feature, XCMS has a parameter (method = "maxint")
that selects which peak intensity to use. This function filters peaks to
match the intensity values that XCMS selected, ensuring consistency with
XCMS's feature intensity matrix.

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

- For large datasets (many features × many samples), the output can be
  memory-intensive. A warning is issued if the expected output exceeds
  10 million rows.

## Examples

``` r
# \donttest{
library(xcms)
#> Loading required package: BiocParallel
#> 
#> This is xcms version 4.8.0 
library(BiocParallel)

# Load example data
faahko_sub <- loadXcmsData("faahko_sub")

# Peak detection
cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())
#> Removed previously identified chromatographic peaks.
#> Detecting mass traces at 25 ppm ... 
#> OK
#> Detecting chromatographic peaks in 2075 regions of interest ...
#>  OK: 561 found.
#> Detecting mass traces at 25 ppm ... 
#> OK
#> Detecting chromatographic peaks in 2308 regions of interest ...
#>  OK: 741 found.
#> Detecting mass traces at 25 ppm ... 
#> OK
#> Detecting chromatographic peaks in 2081 regions of interest ...
#>  OK: 448 found.

# Peak grouping
pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                        minFraction = 0.5)
xdata <- groupChromPeaks(xdata, param = pdp)

# Example 1: Without CAMERA annotations
peak_table_no_camera <- tidy_peaklist(xdata)
head(peak_table_no_camera)
#> # A tibble: 6 × 24
#>   sampleNames fromFile feature_id f_mzmed f_mzmin f_mzmax f_rtmed f_rtmin
#>   <chr>          <dbl>      <int>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1 ko15.CDF           1          1    205     205     205    2788.   2785.
#> 2 ko15.CDF           1          2    206     206     206    2785.   2783.
#> 3 ko15.CDF           1          3    207.    207.    207.   2713.   2711.
#> 4 ko15.CDF           1          4    233     233     233    3025.   3024.
#> 5 ko15.CDF           1          5    233.    233.    233.   3016.   3015.
#> 6 ko15.CDF           1          6    240.    240.    240.   3677.   3675.
#> # ℹ 16 more variables: f_rtmax <dbl>, ms_level <int>, filepath <chr>,
#> #   filename <chr>, peakidx <dbl>, mz <dbl>, mzmin <dbl>, mzmax <dbl>,
#> #   rt <dbl>, rtmin <dbl>, rtmax <dbl>, into <dbl>, intb <dbl>, maxo <dbl>,
#> #   sn <dbl>, is_filled <lgl>

# Example 2: With CAMERA annotations
library(CAMERA)
library(commonMZ)

# Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
xset <- as(xdata, "xcmsSet")
#> Note: you might want to set/adjust the 'sampclass' of the returned xcmSet object before proceeding with the analysis.

# CAMERA annotation
xs <- xsAnnotate(xset, polarity = "positive")
xs <- groupFWHM(xs, perfwhm = 0.1, intval = "into", sigma = 6)
#> Start grouping after retention time.
#> Created 227 pseudospectra.
xs <- findIsotopes(xs, ppm = 10, mzabs = 0.01, intval = "into")
#> Generating peak matrix!
#> Run isotope peak annotation
#>  % finished: 10  20  30  40  50  60  70  80  90  100  
#> Found isotopes: 23 
xs <- groupCorr(xs, calcIso = FALSE, calcCiS = TRUE, calcCaS = FALSE,
                cor_eic_th = 0.7, pval = 1E-6)
#> Start grouping after correlation.
#> Generating EIC's .. 
#> Warning: Found NA peaks in selected sample.
#> 
#> Calculating peak correlations in 227 Groups... 
#>  % finished: 10  20  30  40  50  60  70  80  90  100  
#> 
#> Calculating graph cross linking in 227 Groups... 
#>  % finished: 10  20  30  40  50  60  70  80  90  100  
#> New number of ps-groups:  280 
#> xsAnnotate has now 280 groups, instead of 227 

# Get adduct/fragment rules from commonMZ
rules_pos <- commonMZ::MZ_CAMERA(mode = "pos", warn_clash = TRUE, clash_ppm = 5)
#> Warning: The following adducts/fragments seem to collide. 
#> # A tibble: 2 × 2
#>   first       second                             
#>   <chr>       <chr>                              
#> 1 [M+H-NH3]+  [M+NH4]+                           
#> 2 [M+H-C3H4]+ [M+H+(CH3)2CO-H2O]+ (acetone cond.)
#> 
#> 
#> Consider removing one of them. Example: 
#>  rules=rules[            !grepl("[M+NH4]+",rules[,"name"],fixed=TRUE)         ,]
rules_pos <- as.data.frame(rules_pos)

# Find adducts using the rules
xs <- findAdducts(xs, ppm = 10, mzabs = 0.01, multiplier = 4,
                  polarity = "positive", rules = rules_pos)
#> Generating peak matrix for peak annotation!
#> Polarity is set in xsAnnotate: positive 
#> Found and use user-defined ruleset!
#> Calculating possible adducts in 280 Groups... 
#>  % finished: 10  20  30  40  50  60  70  80  90  100  

# Create long-format peak table
peak_table <- tidy_peaklist(xdata, xs)
head(peak_table)
#> # A tibble: 6 × 27
#>   sampleNames fromFile feature_id f_mzmed f_mzmin f_mzmax f_rtmed f_rtmin
#>   <chr>          <dbl>      <int>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1 ko15.CDF           1          1    205     205     205    2788.   2785.
#> 2 ko15.CDF           1          2    206     206     206    2785.   2783.
#> 3 ko15.CDF           1          3    207.    207.    207.   2713.   2711.
#> 4 ko15.CDF           1          4    233     233     233    3025.   3024.
#> 5 ko15.CDF           1          5    233.    233.    233.   3016.   3015.
#> 6 ko15.CDF           1          6    240.    240.    240.   3677.   3675.
#> # ℹ 19 more variables: f_rtmax <dbl>, ms_level <int>, isotopes <chr>,
#> #   adduct <chr>, pcgroup <int>, filepath <chr>, filename <chr>, peakidx <dbl>,
#> #   mz <dbl>, mzmin <dbl>, mzmax <dbl>, rt <dbl>, rtmin <dbl>, rtmax <dbl>,
#> #   into <dbl>, intb <dbl>, maxo <dbl>, sn <dbl>, is_filled <lgl>

# Explore the data
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:MSnbase':
#> 
#>     combine
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from 'package:generics':
#> 
#>     explain
#> The following objects are masked from 'package:xcms':
#> 
#>     collect, groups
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# Count features per sample
peak_table %>%
  group_by(filename) %>%
  summarize(n_features = sum(!is.na(into)))
#> # A tibble: 3 × 2
#>   filename n_features
#>   <chr>         <int>
#> 1 ko15.CDF        292
#> 2 ko16.CDF        318
#> 3 ko18.CDF        303

# View adduct annotations
peak_table %>%
  filter(!is.na(adduct)) %>%
  select(feature_id, f_mzmed, f_rtmed, adduct, pcgroup)
#> # A tibble: 1,044 × 5
#>    feature_id f_mzmed f_rtmed adduct                                     pcgroup
#>         <int>   <dbl>   <dbl> <chr>                                        <int>
#>  1          1    205    2788. "[M+H-C3H9N-C2H4O2]+ 323.091"                   10
#>  2          2    206    2785. "[M+H-NH3-C3H4]+ 262.041 [M+H-C4H6O]+ 275…      36
#>  3          3    207.   2713. ""                                             164
#>  4          4    233    3025. ""                                             181
#>  5          5    233.   3016. ""                                             208
#>  6          6    240.   3677. ""                                             228
#>  7          7    241.   3670. ""                                              97
#>  8          8    242.   3674. ""                                             202
#>  9          9    244.   2829. ""                                             146
#> 10         10    249.   3664. "[M+H-C4H6-H2O]+ 320.158 [M+H-C4H6-NH3-H2…     111
#> # ℹ 1,034 more rows

# Example 3: With XcmsExperiment object
library(MsExperiment)
library(msdata)
#> 
#> Attaching package: 'msdata'
#> The following object is masked from 'package:dplyr':
#> 
#>     ident

# Load data as XcmsExperiment
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
xdata_exp <- readMsExperiment(spectraFiles = fls[1:2], BPPARAM = SerialParam())

# Peak detection and grouping
cwp <- CentWaveParam(peakwidth = c(5, 20), noise = 100)
xdata_exp <- findChromPeaks(xdata_exp, param = cwp, BPPARAM = SerialParam())
pdp <- PeakDensityParam(sampleGroups = rep(1, 2), minFraction = 0.5)
xdata_exp <- groupChromPeaks(xdata_exp, param = pdp)

# Create peak table without CAMERA
peak_table_exp <- tidy_peaklist(xdata_exp)
head(peak_table_exp)
#> # A tibble: 6 × 25
#>   sample_index spectraOrigin fromFile feature_id f_mzmed f_mzmin f_mzmax f_rtmed
#>          <int> <chr>            <dbl>      <int>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1            1 /home/runner…        1          1    105.    105.    105.    63.1
#> 2            1 /home/runner…        1          2    105.    105.    105.   158. 
#> 3            1 /home/runner…        1          3    105.    105.    105.   203. 
#> 4            1 /home/runner…        1          4    106.    106.    106.   181. 
#> 5            1 /home/runner…        1          5    106.    106.    106.   195. 
#> 6            1 /home/runner…        1          6    106.    106.    106.    30.8
#> # ℹ 17 more variables: f_rtmin <dbl>, f_rtmax <dbl>, ms_level <int>,
#> #   filepath <chr>, filename <chr>, peakidx <dbl>, mz <dbl>, mzmin <dbl>,
#> #   mzmax <dbl>, rt <dbl>, rtmin <dbl>, rtmax <dbl>, into <dbl>, intb <dbl>,
#> #   maxo <dbl>, sn <dbl>, is_filled <lgl>
# }
```
