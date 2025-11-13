#' Create Tidy Long-Format Peak Table with Optional Annotations
#'
#' Generates a comprehensive long-format peak table where each row represents
#' one feature for one sample. Optionally includes CAMERA annotations for isotopes,
#' adducts, and pseudospectrum groups if provided.
#'
#' @param x An [xcms::XCMSnExp-class] or [xcms::XcmsExperiment-class]
#'   object containing peak detection results. Feature grouping (via
#'   [xcms::groupChromPeaks()]) is optional but recommended.
#' @param xsAnnotate Optional. An `xsAnnotate` object from the CAMERA package
#'   containing peak annotations (isotopes, adducts, pseudospectrum groups).
#'   If NULL (default), annotation columns will not be included in the output.
#'
#' @return A [tibble::tibble] in long format with one row per feature per sample
#'   (or one row per peak per sample if features are not defined).
#'   The tibble contains the following columns with their data types:
#'   \describe{
#'     \item{feature_id}{integer: Feature identifier (1, 2, 3, ...). Only present
#'       when features are defined via [xcms::groupChromPeaks()].}
#'     \item{f_mzmed, f_mzmin, f_mzmax}{numeric: Feature-level m/z statistics
#'       (median, minimum, maximum). Only present when features are defined.}
#'     \item{f_rtmed, f_rtmin, f_rtmax}{numeric: Feature-level retention time
#'       statistics in seconds (median, minimum, maximum). Only present when
#'       features are defined.}
#'     \item{ms_level}{integer: MS level (typically 1 for MS1 features). Only
#'       present when features are defined.}
#'     \item{isotopes}{character: Isotope annotations from CAMERA (e.g., "\code{"[M]+"}",
#'       "\code{"[M+1]+"}"). Only present when xsAnnotate is provided.}
#'     \item{adduct}{character: Adduct annotations from CAMERA (e.g., "\code{"[M+H]+"}",
#'       "\code{"[M+Na]+"}"). Only present when xsAnnotate is provided.}
#'     \item{pcgroup}{integer: Pseudospectrum group ID from CAMERA. Features with
#'       the same pcgroup are believed to originate from the same compound. Only
#'       present when xsAnnotate is provided.}
#'     \item{feature_group}{character: Feature group identifier from
#'       [MsFeatures::groupFeatures()] (e.g., "FG.001", "FG.002"). Only present
#'       if groupFeatures was applied.}
#'     \item{mz, mzmin, mzmax}{numeric: Peak m/z values (centroid, minimum, maximum).}
#'     \item{rt, rtmin, rtmax}{numeric: Peak retention time in seconds (apex, minimum,
#'       maximum).}
#'     \item{into}{numeric: Integrated peak intensity (integrated area of original
#'       (raw) peak).}
#'     \item{intb}{numeric: Baseline-corrected integrated peak intensity.}
#'     \item{maxo}{numeric: Maximum peak intensity (apex).}
#'     \item{sn}{numeric: Signal-to-noise ratio.}
#'     \item{is_filled}{logical: Whether the peak was gap-filled (TRUE) or originally
#'       detected (FALSE).}
#'     \item{filepath}{character: Full file path to the raw data file.}
#'     \item{filename}{character: Basename of the raw data file.}
#'     \item{fromFile}{numeric: Sample index (1-based) corresponding to the order
#'       in [xcms::fileNames,MsExperiment-method].}
#'     \item{peakidx}{numeric: Internal peak index from [xcms::chromPeaks()].}
#'     \item{sample_index}{integer: Sample index from sample metadata.}
#'     \item{spectraOrigin}{character: Original spectra file path.}
#'     \item{Additional columns}{Variable types: Any columns from [Biobase::pData()]
#'       or [MsExperiment::sampleData()] are included with their original data types.}
#'   }
#'
#' @details
#' This function integrates data from multiple sources:
#' \itemize{
#'   \item Chromatographic peaks from [xcms::chromPeaks()]
#'   \item Feature definitions from [xcms::featureDefinitions()] (if available)
#'   \item Feature intensities from [xcms::featureValues()] (if features defined)
#'   \item CAMERA annotations from the xsAnnotate object (if provided)
#'   \item Sample metadata from [Biobase::pData()] or [MsExperiment::sampleData()]
#' }
#'
#' The function supports both [xcms::XCMSnExp-class] and [xcms::XcmsExperiment-class]
#' objects, providing flexibility for modern XCMS workflows.
#'
#' **Feature grouping is optional**: If [xcms::groupChromPeaks()] has not been run,
#' the function will return peak-level data only, without feature-level aggregation.
#'
#' **Intensity matching**: When multiple peaks from the same sample are assigned to
#' a single feature, XCMS has a parameter (method = "maxint") that selects which
#' peak intensity to use. This function filters peaks to match the intensity values
#' that XCMS selected, ensuring consistency with XCMS's feature intensity matrix.
#'
#' Missing values (NA) indicate that a feature was not detected in a particular
#' sample.
#'
#' @note
#' \itemize{
#'   \item XCMS sometimes creates duplicate peaks; this function keeps only the
#'     first occurrence.
#'   \item The function uses "maxint" method for feature values, meaning it
#'     takes the maximum intensity peak for each feature in each sample.
#'   \item CAMERA annotations are optional. If xsAnnotate is NULL, the isotopes,
#'     adduct, and pcgroup columns will not be present in the output.
#'   \item MsFeatures grouping is automatically detected. If [MsFeatures::groupFeatures()]
#'     was applied to the object, the feature_group column will be included in the output.
#'   \item For XcmsExperiment objects, sample metadata is accessed via sampleData()
#'     instead of pData().
#'   \item For large datasets (many features x many samples), the output can be
#'     memory-intensive. A warning is issued if the expected output exceeds 10 million rows.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(xcms)
#' library(BiocParallel)
#'
#' # Load example data
#' faahko_sub <- loadXcmsData("faahko_sub")
#'
#' # Peak detection
#' cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
#' xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())
#'
#' # Peak grouping
#' pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
#'                         minFraction = 0.5)
#' xdata <- groupChromPeaks(xdata, param = pdp)
#'
#' # Example 1: Without CAMERA annotations
#' peak_table_no_camera <- tidy_peaklist(xdata)
#' head(peak_table_no_camera)
#'
#' # Example 2: With CAMERA annotations
#' library(CAMERA)
#' library(commonMZ)
#'
#' # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
#' xset <- as(xdata, "xcmsSet")
#'
#' # CAMERA annotation
#' xs <- xsAnnotate(xset, polarity = "positive")
#' xs <- groupFWHM(xs, perfwhm = 0.1, intval = "into", sigma = 6)
#' xs <- findIsotopes(xs, ppm = 10, mzabs = 0.01, intval = "into")
#' xs <- groupCorr(xs, calcIso = FALSE, calcCiS = TRUE, calcCaS = FALSE,
#'                 cor_eic_th = 0.7, pval = 1E-6)
#'
#' # Get adduct/fragment rules from commonMZ
#' rules_pos <- commonMZ::MZ_CAMERA(mode = "pos", warn_clash = TRUE, clash_ppm = 5)
#' rules_pos <- as.data.frame(rules_pos)
#'
#' # Find adducts using the rules
#' xs <- findAdducts(xs, ppm = 10, mzabs = 0.01, multiplier = 4,
#'                   polarity = "positive", rules = rules_pos)
#'
#' # Create long-format peak table
#' peak_table <- tidy_peaklist(xdata, xs)
#' head(peak_table)
#'
#' # Explore the data
#' library(dplyr)
#'
#' # Count features per sample
#' peak_table %>%
#'   group_by(filename) %>%
#'   summarize(n_features = sum(!is.na(into)))
#'
#' # View adduct annotations
#' peak_table %>%
#'   filter(!is.na(adduct)) %>%
#'   select(feature_id, f_mzmed, f_rtmed, adduct, pcgroup)
#'
#' # Example 3: With XcmsExperiment object
#' library(MsExperiment)
#' library(msdata)
#'
#' # Load data as XcmsExperiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' xdata_exp <- readMsExperiment(spectraFiles = fls[1:2], BPPARAM = SerialParam())
#'
#' # Peak detection and grouping
#' cwp <- CentWaveParam(peakwidth = c(5, 20), noise = 100)
#' xdata_exp <- findChromPeaks(xdata_exp, param = cwp, BPPARAM = SerialParam())
#' pdp <- PeakDensityParam(sampleGroups = rep(1, 2), minFraction = 0.5)
#' xdata_exp <- groupChromPeaks(xdata_exp, param = pdp)
#'
#' # Create peak table without CAMERA
#' peak_table_exp <- tidy_peaklist(xdata_exp)
#' head(peak_table_exp)
#' }
#'
#' @importFrom xcms chromPeaks chromPeakData featureDefinitions featureValues fileNames
#' @importFrom MsExperiment sampleData
#' @importFrom dplyr %>% mutate rename left_join right_join filter group_by
#'   ungroup slice select rename_with any_of bind_cols as_tibble across n
#' @importFrom tidyr unnest gather complete nesting
#' @importFrom tibble as_tibble
#' @importFrom Biobase pData
#' @importFrom rlang syms
tidy_peaklist <- function(x, xsAnnotate = NULL) {
  # Input validation
  .validate_inputs(x, xsAnnotate)

  # Extract peaks with file information
  temp_peaks <- .extract_peaks_with_metadata(x)

  # Check if features are defined
  has_features <- nrow(featureDefinitions(x)) > 0

  if (!has_features) {
    # If no features defined, return peak-level data only
    warning("No features defined (groupChromPeaks not run). Returning peak-level data only.")
    return(.add_sample_metadata(temp_peaks, x))
  }

  # Extract and annotate features
  temp_features <- .extract_and_annotate_features(x, xsAnnotate)

  # Extract feature intensities
  temp_peaks_tab <- .extract_feature_intensities(x)

  # Join features with peaks and filter
  out <- .join_and_filter_features_peaks(temp_features, temp_peaks, temp_peaks_tab)

  # Check for potential memory issues
  .check_memory_usage(out)

  # Complete to include all feature-sample combinations
  out <- .complete_feature_sample_combinations(out, temp_features, xsAnnotate, x)

  # Add sample metadata
  out <- .add_sample_metadata(out, x)

  return(out)
}


# Helper functions --------------------------------------------------------


#' Validate input parameters
#' @noRd
.validate_inputs <- function(x, xsAnnotate) {
  # Validate object class
  if (!inherits(x, c("XCMSnExp", "XcmsExperiment"))) {
    stop("'x' must be an XCMSnExp or XcmsExperiment object")
  }

  # Check if peaks exist
  if (nrow(chromPeaks(x)) == 0) {
    stop("No chromatographic peaks found. Run findChromPeaks() first.")
  }

  # Validate xsAnnotate if provided
  if (!is.null(xsAnnotate)) {
    if (!inherits(xsAnnotate, "xsAnnotate")) {
      stop("'xsAnnotate' must be an xsAnnotate object from CAMERA package")
    }

    # Check dimensions match (only if features are defined)
    if (nrow(featureDefinitions(x)) > 0) {
      camera_nfeatures <- nrow(CAMERA::getPeaklist(xsAnnotate))
      xcms_nfeatures <- nrow(featureDefinitions(x))
      if (camera_nfeatures != xcms_nfeatures) {
        stop(sprintf(
          "Feature count mismatch: XCMS has %d features but CAMERA has %d",
          xcms_nfeatures, camera_nfeatures
        ))
      }
    }
  }
}


#' Extract peaks with file path metadata
#' @noRd
.extract_peaks_with_metadata <- function(x) {
  peaks <- .chromPeaks_classic(x)
  n_peaks <- nrow(peaks)

  peaks %>%
    mutate(peakidx = seq_len(n_peaks)) %>%
    rename(fromFile = sample) %>%
    mutate(filepath = fileNames(x)[fromFile],
           filename = basename(filepath))
}


#' Extract and annotate features with CAMERA if provided
#' @noRd
.extract_and_annotate_features <- function(x, xsAnnotate) {
  temp_features <- featureDefinitions(x) %>%
    as_tibble()

  # Add CAMERA annotations if provided
  if (!is.null(xsAnnotate)) {
    temp_features <- .add_camera_annotations(temp_features, xsAnnotate)
  }

  # Add feature IDs and unnest peak indices
  n_features <- nrow(temp_features)
  temp_features <- temp_features %>%
    mutate(feature_id = seq_len(n_features)) %>%
    unnest(peakidx)

  return(temp_features)
}


#' Add CAMERA annotations to features
#' @noRd
.add_camera_annotations <- function(features, xsAnnotate) {
  camera_annot <- as_tibble(CAMERA::getPeaklist(xsAnnotate))

  # Extract CAMERA annotation columns
  camera_cols <- intersect(c("isotopes", "adduct", "pcgroup"), colnames(camera_annot))

  if (length(camera_cols) > 0) {
    features <- features %>%
      bind_cols(camera_annot[, camera_cols]) %>%
      mutate(across(any_of("pcgroup"), as.integer))
  }

  return(features)
}


#' Extract feature intensity values
#' @noRd
.extract_feature_intensities <- function(x) {
  n_features <- nrow(featureDefinitions(x))

  featureValues(x, value = "into", method = "maxint") %>%
    as_tibble() %>%
    mutate(feature_id = seq_len(n_features)) %>%
    gather("filename", "into_f", -feature_id)
}


#' Join features with peaks and filter to match XCMS intensity selection
#' @noRd
.join_and_filter_features_peaks <- function(temp_features, temp_peaks, temp_peaks_tab) {
  out <- temp_features %>%
    select(mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, peakidx, feature_id,
           any_of(c("isotopes", "adduct", "pcgroup", "feature_group"))) %>%
    rename_with(~paste0("f_", .), .cols = any_of(c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax"))) %>%
    left_join(temp_peaks, by = "peakidx") %>%
    left_join(temp_peaks_tab, by = c("feature_id", "filename"))

  # Filter to match peak intensities with feature intensities
  # This ensures we use the same peaks that XCMS selected when method="maxint"
  # Multiple peaks from the same sample can be assigned to one feature, and
  # this filter keeps only the peak(s) whose intensity matches what XCMS reported
  out <- out %>%
    filter(abs(into - into_f) < .Machine$double.eps^0.5) %>%
    select(-into_f) %>%
    group_by(feature_id, filename) %>%
    slice(1) %>%  # Take first if duplicates exist
    ungroup()

  return(out)
}


#' Get feature columns for nesting based on what's available
#' @noRd
.get_feature_columns <- function(data, has_camera, has_features_group) {
  # Base feature columns (should always exist if features are defined)
  base_cols <- c("feature_id", "f_mzmed", "f_mzmin", "f_mzmax",
                 "f_rtmed", "f_rtmin", "f_rtmax", "ms_level")

  # Add CAMERA columns if present
  if (has_camera) {
    base_cols <- c(base_cols, "isotopes", "adduct", "pcgroup")
  }

  # Add feature_group if present
  if (has_features_group) {
    base_cols <- c(base_cols, "feature_group")
  }

  # Return only columns that actually exist in the data
  intersect(base_cols, colnames(data))
}


#' Check for potential memory issues with large datasets
#' @noRd
.check_memory_usage <- function(data) {
  n_features <- length(unique(data$feature_id))
  n_samples <- length(unique(data$filename))
  expected_rows <- n_features * n_samples

  # Warn if expected output will be very large (> 10 million rows)
  if (expected_rows > 1e7) {
    warning(sprintf(
      "Creating large data frame with %s rows (%d features x %d samples). Consider filtering features first.",
      format(expected_rows, big.mark = ","),
      n_features,
      n_samples
    ))
  }
}


#' Complete data to include all feature-sample combinations
#' @noRd
.complete_feature_sample_combinations <- function(out, temp_features, xsAnnotate, x) {
  # Determine which optional columns are present
  has_camera <- !is.null(xsAnnotate)
  has_feature_group <- "feature_group" %in% colnames(temp_features)

  # Get the columns to use for nesting
  feature_cols <- .get_feature_columns(out, has_camera, has_feature_group)

  # Complete the data frame
  out <- out %>%
    complete(nesting(!!!syms(feature_cols)),
             nesting(filepath, filename, fromFile))

  return(out)
}


#' Add sample metadata to the output
#' @noRd
.add_sample_metadata <- function(data, x) {
  # Handle both XCMSnExp and XcmsExperiment
  if (inherits(x, "XcmsExperiment")) {
    # For XcmsExperiment, use sampleData()
    sample_metadata <- as.data.frame(sampleData(x)) %>%
      mutate(fromFile = seq_len(n())) %>%
      as_tibble()
  } else {
    # For XCMSnExp, use pData()
    sample_metadata <- pData(x) %>%
      mutate(fromFile = seq_len(n())) %>%
      as_tibble()
  }

  # Join sample metadata
  sample_metadata %>%
    right_join(data, by = "fromFile", multiple = "all")
}
