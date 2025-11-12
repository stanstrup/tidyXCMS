# Creating Long-Format Peak Tables with Optional CAMERA Annotations

## Introduction

The `tidyXCMS` package provides functions to work with XCMS metabolomics
data in a tidy, long-format structure. This vignette demonstrates how to
use the
[`XCMSnExp_CAMERA_peaklist_long()`](https://stanstrup.github.io/tidyXCMS/reference/XCMSnExp_CAMERA_peaklist_long.md)
function to create a comprehensive peak table that integrates:

- Peak detection results from XCMS
- Feature grouping (correspondence) information
- CAMERA annotations (isotopes, adducts, pseudospectrum groups) -
  **optional**
- Sample metadata

The resulting long-format table has one row per feature per sample,
making it ideal for downstream analysis with tidyverse tools like
`dplyr` and `ggplot2`.

**Note:** CAMERA annotations are optional. You can create a peak table
without CAMERA annotations by simply omitting the `xsAnnotate`
parameter.

## Setup

``` r
library(tidyXCMS)
library(xcms)
library(CAMERA)
library(BiocParallel)
library(dplyr)
library(ggplot2)
```

## Example Workflow

### Load Example Data

We’ll use the `faahko_sub` dataset from XCMS, which contains LC-MS data
from wild-type and knockout samples.

``` r
# Load example data from XCMS
faahko_sub <- loadXcmsData("faahko_sub")

# Check the files
fileNames(faahko_sub)
#> [1] "/home/runner/work/_temp/Library/faahKO/cdf/KO/ko15.CDF"
#> [2] "/home/runner/work/_temp/Library/faahKO/cdf/KO/ko16.CDF"
#> [3] "/home/runner/work/_temp/Library/faahKO/cdf/KO/ko18.CDF"
```

### Add Sample Metadata

Let’s add some sample information that will be included in our final
table.

``` r
library(Biobase)

pd <- data.frame(
  sample_name = basename(fileNames(faahko_sub)),
  sample_group = c("WT", "KO", "WT")[1:length(fileNames(faahko_sub))],
  injection_order = 1:length(fileNames(faahko_sub)),
  row.names = basename(fileNames(faahko_sub))
)

pData(faahko_sub) <- pd
```

### Peak Detection

Perform chromatographic peak detection using the CentWave algorithm.

``` r
# Define peak detection parameters
cwp <- CentWaveParam(
  peakwidth = c(20, 80),
  noise = 5000,
  ppm = 25
)

# Find peaks (using SerialParam for reproducibility)
xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())

# Check results
cat("Detected", nrow(chromPeaks(xdata)), "peaks across",
    length(fileNames(xdata)), "samples\n")
#> Detected 1750 peaks across 3 samples
```

### Feature Grouping (Correspondence)

Group peaks across samples into features (aligned peaks).

``` r
# Define grouping parameters
pdp <- PeakDensityParam(
  sampleGroups = pData(xdata)$sample_group,
  minFraction = 0.5,
  bw = 30
)

# Group peaks
xdata <- groupChromPeaks(xdata, param = pdp)

# Check results
cat("Grouped peaks into", nrow(featureDefinitions(xdata)), "features\n")
#> Grouped peaks into 842 features
```

### Create Long-Format Peak Table (Without CAMERA)

You can create a long-format peak table directly from the XCMS results
without using CAMERA annotations. This is useful when you want to
quickly explore your data or when CAMERA annotations are not needed for
your analysis.

``` r
# Create peak table without CAMERA
peak_table_no_camera <- XCMSnExp_CAMERA_peaklist_long(xdata)

# Check the structure
dim(peak_table_no_camera)
#> [1] 2526   26

# View first few rows
head(peak_table_no_camera)
#> # A tibble: 6 × 26
#>   sample_name sample_group injection_order fromFile feature_id f_mzmed f_mzmin
#>   <chr>       <chr>                  <int>    <dbl>      <int>   <dbl>   <dbl>
#> 1 ko15.CDF    WT                         1        1          1    200.    200.
#> 2 ko15.CDF    WT                         1        1          2    201.    201.
#> 3 ko15.CDF    WT                         1        1          3    205     205 
#> 4 ko15.CDF    WT                         1        1          4    206.    206.
#> 5 ko15.CDF    WT                         1        1          5    206     206 
#> 6 ko15.CDF    WT                         1        1          6    207.    207.
#> # ℹ 19 more variables: f_mzmax <dbl>, f_rtmed <dbl>, f_rtmin <dbl>,
#> #   f_rtmax <dbl>, ms_level <int>, filepath <chr>, filename <chr>,
#> #   peakidx <dbl>, mz <dbl>, mzmin <dbl>, mzmax <dbl>, rt <dbl>, rtmin <dbl>,
#> #   rtmax <dbl>, into <dbl>, intb <dbl>, maxo <dbl>, sn <dbl>, is_filled <lgl>

# Check column names - note that isotopes, adduct, and pcgroup are NOT present
colnames(peak_table_no_camera)
#>  [1] "sample_name"     "sample_group"    "injection_order" "fromFile"       
#>  [5] "feature_id"      "f_mzmed"         "f_mzmin"         "f_mzmax"        
#>  [9] "f_rtmed"         "f_rtmin"         "f_rtmax"         "ms_level"       
#> [13] "filepath"        "filename"        "peakidx"         "mz"             
#> [17] "mzmin"           "mzmax"           "rt"              "rtmin"          
#> [21] "rtmax"           "into"            "intb"            "maxo"           
#> [25] "sn"              "is_filled"
```

The resulting table contains all feature and peak information. The
CAMERA annotation columns (`isotopes`, `adduct`, `pcgroup`) are **not
included** when no CAMERA annotations are provided, keeping the output
cleaner and more memory-efficient.

### CAMERA Annotation (Optional)

Use CAMERA to annotate isotopes, adducts, and group features into
pseudospectra.

**Note:** CAMERA requires an `xcmsSet` object, so we first convert the
`XCMSnExp` object using `as("xcmsSet")`.

#### Creating a CAMERA object

In this first step we create a CAMERA object from the XCMS object. We
specify the polarity to help CAMERA determine appropriate adducts and
fragments.

``` r
# Convert to xcmsSet for CAMERA (required for compatibility)
xset <- as(xdata, "xcmsSet")

# Create xsAnnotate object with polarity
xs <- xsAnnotate(xset, polarity = "positive")
```

#### Grouping coeluting peaks

The first step to grouping features is to group co-eluting peaks. This
is a naïve approach that we will refine later.

``` r
# Group peaks by retention time
xs <- groupFWHM(xs, perfwhm = 0.1, intval = "into", sigma = 6)
#> Start grouping after retention time.
#> Created 465 pseudospectra.
```

#### Grouping based on correlation

Now we group the features based on correlations. This looks at each
group from the previous step and splits them into separate groups for
peaks that correlate with each other.

**Note:** For this small example dataset, we use a simplified approach.
In real workflows with more samples, you would use stricter parameters
and enable correlation across samples (`calcCaS = TRUE`).

``` r
# Group by correlation
# For small datasets, we skip correlation analysis to avoid errors
# In production with adequate samples, use:
# xs <- groupCorr(xs, calcCaS = TRUE, cor_eic_th = 0.7, pval = 0.05)

# For this example, we skip the correlation step and proceed directly
# xs <- groupCorr(xs, calcIso = FALSE, calcCiS = FALSE, calcCaS = FALSE)
```

#### Isotope annotation

This annotates peaks that are possible isotopes based on m/z difference
and intensity patterns.

``` r
# Find isotopes
xs <- findIsotopes(xs, ppm = 10, mzabs = 0.01, intval = "into", maxcharge = 2)
#> Generating peak matrix!
#> Run isotope peak annotation
#>  % finished: 10  20  30  40  50  60  70  80  90  100  
#> Found isotopes: 39
```

#### Annotation of adducts and fragments

Now we try to annotate adducts and fragments based on expected mass
differences.

``` r
# Find adducts
xs <- findAdducts(xs, ppm = 10, mzabs = 0.01, multiplier = 4, polarity = "positive")
#> Generating peak matrix for peak annotation!
#> Polarity is set in xsAnnotate: positive 
#> Ruleset could not read from object! Recalculate
#> 
#> Calculating possible adducts in 465 Groups... 
#>  % finished: 10  20  30  40  50  60  70  80  90  100
```

### Create Long-Format Peak Table

Now we can create our comprehensive long-format peak table! Note that we
use the original `xdata` (XCMSnExp) object along with the CAMERA
annotations from `xs`.

``` r
peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata, xs)

# Check the structure
dim(peak_table)
#> [1] 2526   29
colnames(peak_table)
#>  [1] "sample_name"     "sample_group"    "injection_order" "fromFile"       
#>  [5] "feature_id"      "f_mzmed"         "f_mzmin"         "f_mzmax"        
#>  [9] "f_rtmed"         "f_rtmin"         "f_rtmax"         "isotopes"       
#> [13] "adduct"          "pcgroup"         "ms_level"        "filepath"       
#> [17] "filename"        "peakidx"         "mz"              "mzmin"          
#> [21] "mzmax"           "rt"              "rtmin"           "rtmax"          
#> [25] "into"            "intb"            "maxo"            "sn"             
#> [29] "is_filled"
```

## Understanding the Output

The resulting tibble contains comprehensive information about each
feature in each sample. Below is a complete description of all columns:

### Column Reference Table

The table has been organized into logical groups for easier
understanding:

| Column                                                               | Content                                                                                           |
|----------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| **Sample information**                                               |                                                                                                   |
| `filepath`                                                           | Path to the raw data file                                                                         |
| `filename`                                                           | The filename without path                                                                         |
| `fromFile`                                                           | The file number (the order files were supplied in)                                                |
| *Plus any columns from `pData(XCMSnExp)`*                            | Sample metadata columns you added                                                                 |
| **Feature identifiers**                                              |                                                                                                   |
| `feature_id`                                                         | The index of the feature after grouping across samples                                            |
| `peakidx`                                                            | The index of the peak before grouping across samples                                              |
| **Feature-level m/z statistics (across all samples)**                |                                                                                                   |
| `f_mzmed`                                                            | The median m/z found for that feature **across samples**                                          |
| `f_mzmin`                                                            | The minimum m/z found for that feature across samples                                             |
| `f_mzmax`                                                            | The maximum m/z found for that feature across samples                                             |
| **Feature-level retention time statistics (across all samples)**     |                                                                                                   |
| `f_rtmed`                                                            | The median retention time found for that feature across samples                                   |
| `f_rtmin`                                                            | The minimum retention time found for that feature across samples                                  |
| `f_rtmax`                                                            | The maximum retention time found for that feature across samples                                  |
| **CAMERA annotations**                                               |                                                                                                   |
| `isotopes`                                                           | The isotope annotation from CAMERA (e.g., “\[M\]+”, “\[M+1\]+”, “\[M+2\]+”)                       |
| `adduct`                                                             | The adduct annotation from CAMERA (e.g., “\[M+H\]+”, “\[M+Na\]+”, “\[M+NH4\]+”)                   |
| `pcgroup`                                                            | The feature grouping index from CAMERA - features with same ID likely come from the same compound |
| **Peak-level m/z measurements (in that specific sample)**            |                                                                                                   |
| `mz`                                                                 | The median m/z found for that feature in **that sample**                                          |
| `mzmin`                                                              | The minimum m/z found for that feature in that sample                                             |
| `mzmax`                                                              | The maximum m/z found for that feature in that sample                                             |
| **Peak-level retention time measurements (in that specific sample)** |                                                                                                   |
| `rt`                                                                 | The retention time found for that feature in that sample                                          |
| `rtmin`                                                              | The minimum retention time found for that feature in that sample                                  |
| `rtmax`                                                              | The maximum retention time found for that feature in that sample                                  |
| **Peak intensity measurements**                                      |                                                                                                   |
| `into`                                                               | The area under the peak                                                                           |
| `intb`                                                               | The area under the peak after baseline removal                                                    |
| `maxo`                                                               | The maximum intensity (i.e., height) of the peak                                                  |
| `sn`                                                                 | The signal to noise ratio of that peak                                                            |
| **Gaussian peak fitting parameters**                                 |                                                                                                   |
| `egauss`                                                             | RMSE of Gaussian fit                                                                              |
| `mu`                                                                 | Gaussian parameter mu (center of the Gaussian; unit is scan number)                               |
| `sigma`                                                              | Gaussian parameter sigma                                                                          |
| `h`                                                                  | Gaussian parameter h (height of the Gaussian peak)                                                |
| **CentWave algorithm parameters**                                    |                                                                                                   |
| `f`                                                                  | Region number of m/z ROI where the peak was localized                                             |
| `dppm`                                                               | m/z deviation of mass trace across scans in ppm                                                   |
| `scale`                                                              | Scale on which the peak was localized                                                             |
| `scpos`                                                              | Center of peak position found by wavelet analysis                                                 |
| `scmin`                                                              | Left peak limit found by wavelet analysis (scan number)                                           |
| `scmax`                                                              | Right peak limit found by wavelet analysis (scan number)                                          |
| **Additional columns**                                               |                                                                                                   |
| `ms_level`                                                           | MS level (e.g., MS1, MS2)                                                                         |
| `is_filled`                                                          | Was the intensity found by gap filling (TRUE) or peak picking (FALSE)                             |

**Important notes:**

- **Each row represents one feature in one sample**
- If a feature was not detected in a sample, peak-level columns (mz, rt,
  into, etc.) will be `NA`
- Feature-level statistics (f_mzmed, f_rtmed, etc.) are always present
  and represent values across all samples
- CAMERA annotations apply at the feature level and are the same across
  all samples for a given feature
- The distinction between “feature-level” and “peak-level” is important:
  - **Feature-level** (prefix `f_`): Statistics aggregated across all
    samples
  - **Peak-level** (no prefix): Values for the specific peak in that
    specific sample

### Feature-Level Information

These columns describe each feature (grouped peak):

``` r
peak_table %>%
  select(feature_id, f_mzmed, f_rtmed) %>%
  distinct() %>%
  head()
#> # A tibble: 6 × 3
#>   feature_id f_mzmed f_rtmed
#>        <int>   <dbl>   <dbl>
#> 1          1    200.   2877.
#> 2          2    201.   3141.
#> 3          3    205    2788.
#> 4          4    206.   2796.
#> 5          5    206    2785.
#> 6          6    207.   2713.
```

### CAMERA Annotations

- `isotopes`: Isotope annotation (e.g., “\[M\]+”, “\[M+1\]+”)
- `adduct`: Adduct annotation (e.g., “\[M+H\]+”, “\[M+Na\]+”)
- `pcgroup`: Pseudospectrum correlation group ID

``` r
# View features with adduct annotations
peak_table %>%
  filter(!is.na(adduct)) %>%
  select(feature_id, f_mzmed, f_rtmed, isotopes, adduct, pcgroup) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  head()
#> # A tibble: 6 × 6
#>   feature_id f_mzmed f_rtmed isotopes adduct pcgroup
#>        <int>   <dbl>   <dbl> <chr>    <chr>    <int>
#> 1          1    200.   2877. ""       ""         228
#> 2          2    201.   3141. ""       ""         227
#> 3          3    205    2788. ""       ""          10
#> 4          4    206.   2796. ""       ""         457
#> 5          5    206    2785. ""       ""          41
#> 6          6    207.   2713. ""       ""         245
```

### Peak-Level Information

For each feature in each sample:

- `mz`, `rt`: Detected peak position
- `into`: Integrated peak intensity
- `intb`: Baseline-corrected intensity
- `maxo`: Maximum intensity
- `sn`: Signal-to-noise ratio

``` r
# View detected peaks (non-NA intensities)
peak_table %>%
  filter(!is.na(into)) %>%
  select(feature_id, filename, mz, rt, into, sn) %>%
  head()
#> # A tibble: 6 × 6
#>   feature_id filename    mz    rt     into    sn
#>        <int> <chr>    <dbl> <dbl>    <dbl> <dbl>
#> 1          3 ko15.CDF  205  2785. 1924712.    64
#> 2          4 ko15.CDF  206. 2796.   57359.  8177
#> 3          5 ko15.CDF  206  2783.  213659.    14
#> 4          6 ko15.CDF  207. 2711.  349011.    17
#> 5          7 ko15.CDF  221. 3079.  163255.  7841
#> 6          9 ko15.CDF  231  3079.  328034. 10059
```

### Sample Information

- `filename`, `filepath`: Sample file information
- `fromFile`: Sample index
- Plus any columns from
  [`pData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html) (sample
  metadata)

``` r
peak_table %>%
  select(feature_id, filename, sample_group, injection_order) %>%
  head()
#> # A tibble: 6 × 4
#>   feature_id filename sample_group injection_order
#>        <int> <chr>    <chr>                  <int>
#> 1          1 ko15.CDF WT                         1
#> 2          2 ko15.CDF WT                         1
#> 3          3 ko15.CDF WT                         1
#> 4          4 ko15.CDF WT                         1
#> 5          5 ko15.CDF WT                         1
#> 6          6 ko15.CDF WT                         1
```

### Missing Values

Features not detected in a sample have `NA` for peak-level columns:

``` r
# Count detected features per sample
peak_table %>%
  group_by(filename) %>%
  summarise(
    n_features_detected = sum(!is.na(into)),
    n_features_total = n()
  )
#> # A tibble: 3 × 3
#>   filename n_features_detected n_features_total
#>   <chr>                  <int>            <int>
#> 1 ko15.CDF                 436              842
#> 2 ko16.CDF                 606              842
#> 3 ko18.CDF                 365              842
```

## Downstream Analysis Examples

### Visualize Feature Intensities

``` r
# Plot intensity distribution by sample group
peak_table %>%
  filter(!is.na(into)) %>%
  ggplot(aes(x = sample_group, y = log10(into), fill = sample_group)) +
  geom_boxplot() +
  labs(
    title = "Feature Intensity Distribution",
    x = "Sample Group",
    y = "log10(Intensity)"
  ) +
  theme_minimal()
```

![Boxplot showing the distribution of log10-transformed feature
intensities across sample groups. The plot displays the median,
quartiles, and outliers for each
group.](long-format-peaklist_files/figure-html/plot_intensities-1.png)

### Identify Features with Adduct Annotations

``` r
# Count features by adduct type
peak_table %>%
  filter(!is.na(adduct)) %>%
  distinct(feature_id, adduct) %>%
  count(adduct, sort = TRUE)
#> # A tibble: 75 × 2
#>    adduct                       n
#>    <chr>                    <int>
#>  1 ""                         768
#>  2 "[2M+Na+2K]3+ 590.193"       1
#>  3 "[2M+Na+3K-H]3+ 589.198"     1
#>  4 "[3M+2Na]2+ 210.138"         1
#>  5 "[3M+H+Na]2+ 384.195"        1
#>  6 "[M+2H-CH4]2+ 524.213"       1
#>  7 "[M+2Na-H]+ 384.234"         1
#>  8 "[M+H+NH3]+ 304.113"         1
#>  9 "[M+H-C2H4]+ 328.015"        1
#> 10 "[M+H-C2H4]+ 519.288"        1
#> # ℹ 65 more rows
```

### Find Features Present in All Samples

``` r
# Features detected in all samples
complete_features <- peak_table %>%
  group_by(feature_id) %>%
  summarise(
    n_detected = sum(!is.na(into)),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  filter(n_detected == n_samples)

cat("Found", nrow(complete_features), "features detected in all samples\n")
#> Found 217 features detected in all samples
```

### Coefficient of Variation Analysis

``` r
# Calculate CV for each feature
feature_cv <- peak_table %>%
  filter(!is.na(into)) %>%
  group_by(feature_id, f_mzmed, f_rtmed) %>%
  summarise(
    mean_intensity = mean(into),
    sd_intensity = sd(into),
    cv = sd_intensity / mean_intensity * 100,
    .groups = "drop"
  ) %>%
  arrange(cv)

# Show most stable features
head(feature_cv)
#> # A tibble: 6 × 6
#>   feature_id f_mzmed f_rtmed mean_intensity sd_intensity    cv
#>        <int>   <dbl>   <dbl>          <dbl>        <dbl> <dbl>
#> 1        498    486.   2798.        227942.         246. 0.108
#> 2        454    466.   3193.        173672.         607. 0.350
#> 3        245    369.   3770.        996200.        6299. 0.632
#> 4        286    386.   3446.        205233.        1332. 0.649
#> 5        235    366.   2767.       1096511.       14736. 1.34 
#> 6        108    311.   3070.        277669.        4432. 1.60
```

### Extract Specific Features for Further Analysis

``` r
# Get a specific feature across all samples
feature_123 <- peak_table %>%
  filter(feature_id == 123) %>%
  select(feature_id, filename, sample_group, mz, rt, into)

head(feature_123)
#> # A tibble: 3 × 6
#>   feature_id filename sample_group    mz    rt    into
#>        <int> <chr>    <chr>        <dbl> <dbl>   <dbl>
#> 1        123 ko15.CDF WT             NA    NA      NA 
#> 2        123 ko16.CDF KO            316. 3450. 167161.
#> 3        123 ko18.CDF WT             NA    NA      NA
```

## Session Info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] MSnbase_2.36.0      ProtGenerics_1.42.0 S4Vectors_0.48.0   
#>  [4] mzR_2.44.0          Rcpp_1.1.0          ggplot2_4.0.0      
#>  [7] dplyr_1.1.4         CAMERA_1.66.0       Biobase_2.70.0     
#> [10] BiocGenerics_0.56.0 generics_0.1.4      xcms_4.8.0         
#> [13] BiocParallel_1.44.0 tidyXCMS_0.99.2    
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.2.3                   RBGL_1.86.0                
#>   [3] gridExtra_2.3               rlang_1.1.6                
#>   [5] magrittr_2.0.4              clue_0.3-66                
#>   [7] MassSpecWavelet_1.76.0      matrixStats_1.5.0          
#>   [9] compiler_4.5.2              systemfonts_1.3.1          
#>  [11] vctrs_0.6.5                 reshape2_1.4.4             
#>  [13] stringr_1.6.0               pkgconfig_2.0.3            
#>  [15] MetaboCoreUtils_1.18.0      crayon_1.5.3               
#>  [17] fastmap_1.2.0               backports_1.5.0            
#>  [19] XVector_0.50.0              labeling_0.4.3             
#>  [21] utf8_1.2.6                  rmarkdown_2.30             
#>  [23] graph_1.88.0                preprocessCore_1.72.0      
#>  [25] ragg_1.5.0                  purrr_1.2.0                
#>  [27] xfun_0.54                   MultiAssayExperiment_1.36.0
#>  [29] cachem_1.1.0                jsonlite_2.0.0             
#>  [31] progress_1.2.3              DelayedArray_0.36.0        
#>  [33] prettyunits_1.2.0           parallel_4.5.2             
#>  [35] cluster_2.1.8.1             R6_2.6.1                   
#>  [37] bslib_0.9.0                 stringi_1.8.7              
#>  [39] RColorBrewer_1.1-3          limma_3.66.0               
#>  [41] rpart_4.1.24                GenomicRanges_1.62.0       
#>  [43] jquerylib_0.1.4             Seqinfo_1.0.0              
#>  [45] SummarizedExperiment_1.40.0 iterators_1.0.14           
#>  [47] knitr_1.50                  base64enc_0.1-3            
#>  [49] IRanges_2.44.0              BiocBaseUtils_1.12.0       
#>  [51] nnet_7.3-20                 Matrix_1.7-4               
#>  [53] igraph_2.2.1                tidyselect_1.2.1           
#>  [55] rstudioapi_0.17.1           abind_1.4-8                
#>  [57] yaml_2.3.10                 doParallel_1.0.17          
#>  [59] codetools_0.2-20            affy_1.88.0                
#>  [61] lattice_0.22-7              tibble_3.3.0               
#>  [63] plyr_1.8.9                  withr_3.0.2                
#>  [65] S7_0.2.0                    evaluate_1.0.5             
#>  [67] foreign_0.8-90              desc_1.4.3                 
#>  [69] Spectra_1.20.0              pillar_1.11.1              
#>  [71] affyio_1.80.0               BiocManager_1.30.26        
#>  [73] MatrixGenerics_1.22.0       checkmate_2.3.3            
#>  [75] foreach_1.5.2               MALDIquant_1.22.3          
#>  [77] ncdf4_1.24                  hms_1.1.4                  
#>  [79] scales_1.4.0                MsExperiment_1.12.0        
#>  [81] glue_1.8.0                  Hmisc_5.2-4                
#>  [83] MsFeatures_1.18.0           lazyeval_0.2.2             
#>  [85] tools_4.5.2                 mzID_1.48.0                
#>  [87] data.table_1.17.8           QFeatures_1.20.0           
#>  [89] vsn_3.78.0                  fs_1.6.6                   
#>  [91] XML_3.99-0.20               grid_4.5.2                 
#>  [93] impute_1.84.0               tidyr_1.3.1                
#>  [95] colorspace_2.1-2            MsCoreUtils_1.21.0         
#>  [97] PSMatch_1.14.0              htmlTable_2.4.3            
#>  [99] Formula_1.2-5               cli_3.6.5                  
#> [101] textshaping_1.0.4           S4Arrays_1.10.0            
#> [103] AnnotationFilter_1.34.0     pcaMethods_2.2.0           
#> [105] gtable_0.3.6                sass_0.4.10                
#> [107] digest_0.6.37               SparseArray_1.10.1         
#> [109] htmlwidgets_1.6.4           farver_2.1.2               
#> [111] htmltools_0.5.8.1           pkgdown_2.2.0              
#> [113] lifecycle_1.0.4             statmod_1.5.1              
#> [115] MASS_7.3-65
```
