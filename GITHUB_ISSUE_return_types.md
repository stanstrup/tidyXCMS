# GitHub Issue: Document explicit column types in tidy_peaklist return value

**Title:** Document explicit column types in tidy_peaklist return value

**Labels:** documentation, enhancement

**Body:**

## Issue Description

The `tidy_peaklist()` function returns a tibble but does not explicitly document the data types of each column in the return value documentation.

## Problem

Users don't know if columns are character, numeric, integer, or factor. This creates uncertainty and potential compatibility issues in downstream code if column types change between versions.

## Proposed Solution

Update the roxygen2 documentation for `tidy_peaklist()` to explicitly specify column types in the @return section:

```r
#' @return A [tibble::tibble] with the following columns:
#'   \describe{
#'     \item{feature_id}{integer: Feature identifier (1 to n_features)}
#'     \item{f_mzmed}{numeric: Median m/z of feature (Da)}
#'     \item{f_rtmed}{numeric: Median retention time of feature (seconds)}
#'     \item{f_mzmin}{numeric: Minimum m/z of feature (Da)}
#'     \item{f_mzmax}{numeric: Maximum m/z of feature (Da)}
#'     \item{f_rtmin}{numeric: Minimum retention time (seconds)}
#'     \item{f_rtmax}{numeric: Maximum retention time (seconds)}
#'     \item{isotopes}{character: Isotope annotation (e.g., "[M]+", "[M+1]+"), only if xsAnnotate provided}
#'     \item{adduct}{character: Adduct annotation (e.g., "[M+H]+"), only if xsAnnotate provided}
#'     \item{pcgroup}{integer: CAMERA pseudospectrum group ID, only if xsAnnotate provided}
#'     \item{feature_group}{character: MsFeatures group ID (e.g., "FG.001"), only if groupFeatures applied}
#'     \item{mz}{numeric: Peak m/z in this sample (Da)}
#'     \item{rt}{numeric: Peak retention time in this sample (seconds)}
#'     \item{into}{numeric: Integrated peak intensity (area under curve)}
#'     \item{intb}{numeric: Baseline-corrected peak intensity}
#'     \item{maxo}{numeric: Maximum peak intensity (height)}
#'     \item{sn}{numeric: Signal-to-noise ratio}
#'     \item{filepath}{character: Full path to raw data file}
#'     \item{filename}{character: Basename of raw data file}
#'     \item{fromFile}{integer: File index in input object}
#'     \item{...}{Additional columns from pData()/sampleData() with various types}
#'   }
```

## Benefits

1. **Type stability**: Documents expected types, making breaking changes more visible
2. **User confidence**: Clear expectations for downstream code
3. **Better IDE support**: Tools can provide better autocomplete/hints
4. **Prevents bugs**: Users know what type conversions are needed

## Priority

**Low-Medium** - Improves documentation quality but doesn't affect functionality.

## Related Files

- `R/tidy_peaklist.R` (main function documentation)

---

**To create this issue on GitHub:**
```bash
cd /home/user/tidyXCMS
gh issue create --title "Document explicit column types in tidy_peaklist return value" \
  --body-file GITHUB_ISSUE_return_types.md \
  --label "documentation,enhancement"
```

Or create it manually at: https://github.com/stanstrup/tidyXCMS/issues/new
