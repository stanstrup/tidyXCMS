#' Create Long-Format Peak Table with Optional CAMERA Annotations
#'
#' Generates a comprehensive long-format peak table where each row represents
#' one feature for one sample. Optionally includes CAMERA annotations for isotopes,
#' adducts, and pseudospectrum groups if provided.
#'
#' @param XCMSnExp An [xcms::XCMSnExp-class] or [xcms::XcmsExperiment-class]
#'   object containing peak detection and feature grouping results.
#' @param xsAnnotate Optional. An `xsAnnotate` object from the CAMERA package
#'   containing peak annotations (isotopes, adducts, pseudospectrum groups).
#'   If NULL (default), annotation columns will not be included in the output.
#'
#' @return A [tibble::tibble] in long format with one row per feature per sample.
#'   The tibble contains:
#'   \describe{
#'     \item{Feature-level columns}{f_mzmed, f_mzmin, f_mzmax, f_rtmed, f_rtmin,
#'       f_rtmax - feature summary statistics}
#'     \item{CAMERA annotations}{isotopes, adduct, pcgroup (pseudospectrum group) -
#'       only present when xsAnnotate is provided}
#'     \item{Peak-level columns}{mz, rt, into, intb, maxo, sn - individual peak
#'       measurements}
#'     \item{Sample information}{filepath, filename, fromFile - sample identifiers}
#'     \item{Additional columns}{Any columns from pData(XCMSnExp) or sampleData(XcmsExperiment)}
#'   }
#'
#' @details
#' This function integrates data from multiple sources:
#' \itemize{
#'   \item Chromatographic peaks from [xcms::chromPeaks()]
#'   \item Feature definitions from [xcms::featureDefinitions()]
#'   \item Feature intensities from [xcms::featureValues()]
#'   \item CAMERA annotations from the xsAnnotate object (if provided)
#'   \item Sample metadata from [Biobase::pData()]
#' }
#'
#' The function supports both [xcms::XCMSnExp-class] and [xcms::XcmsExperiment-class]
#' objects, providing flexibility for modern XCMS workflows.
#'
#' The function performs several data processing steps:
#' \enumerate{
#'   \item Extracts peaks and adds file path information
#'   \item Combines feature definitions with CAMERA annotations
#'   \item Unnests peak indices to create feature-peak relationships
#'   \item Joins feature intensities to match peaks with features
#'   \item Filters to ensure peak intensities match feature intensities
#'   \item Handles duplicate peaks by taking the first occurrence
#'   \item Completes the data frame to include all feature-sample combinations
#'   \item Adds sample metadata from pData
#' }
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
#'   \item For XcmsExperiment objects, sample metadata is accessed via sampleData()
#'     instead of pData().
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' peak_table_no_camera <- XCMSnExp_CAMERA_peaklist_long(xdata)
#' head(peak_table_no_camera)
#'
#' # Example 2: With CAMERA annotations
#' library(CAMERA)
#'
#' # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
#' xset <- as(xdata, "xcmsSet")
#'
#' # CAMERA annotation
#' xs <- xsAnnotate(xset)
#' xs <- groupFWHM(xs)
#' xs <- findIsotopes(xs)
#' xs <- groupCorr(xs)
#' xs <- findAdducts(xs)
#'
#' # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
#' peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata, xs)
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
#' peak_table_exp <- XCMSnExp_CAMERA_peaklist_long(xdata_exp)
#' head(peak_table_exp)
#' }
#'
#' @importFrom xcms chromPeaks chromPeakData featureDefinitions featureValues fileNames
#' @importFrom MsExperiment sampleData
#' @importFrom dplyr %>% mutate rename left_join right_join filter group_by
#'   ungroup slice select rename_with any_of bind_cols as_tibble across
#' @importFrom tidyr unnest gather complete nesting
#' @importFrom tibble as_tibble
#' @importFrom Biobase pData
XCMSnExp_CAMERA_peaklist_long <- function(XCMSnExp, xsAnnotate = NULL) {
  # Extract peaks with file information
  temp_peaks <- .chromPeaks_classic(XCMSnExp) %>%
    mutate(peakidx = 1:nrow(.)) %>%
    rename(fromFile = sample) %>%
    mutate(filepath = fileNames(XCMSnExp)[fromFile])

  # Extract features and optionally combine with CAMERA annotations
  temp_features <- featureDefinitions(XCMSnExp) %>%
    as_tibble()

  # Add CAMERA annotations if provided
  if (!is.null(xsAnnotate)) {
    camera_annot <- as_tibble(CAMERA::getPeaklist(xsAnnotate))
    # Extract CAMERA annotation columns
    camera_cols <- intersect(c("isotopes", "adduct", "pcgroup"), colnames(camera_annot))
    if (length(camera_cols) > 0) {
      temp_features <- temp_features %>%
        bind_cols(camera_annot[, camera_cols]) %>%
        mutate(across(any_of("pcgroup"), as.integer))
    }
  }

  temp_features <- temp_features %>%
    mutate(feature_id = 1:nrow(.)) %>%
    unnest(peakidx)

  # Extract feature values (intensities)
  temp_peaks_tab <- featureValues(XCMSnExp, value = "into", method = "maxint") %>%
    as_tibble() %>%
    mutate(feature_id = 1:nrow(.)) %>%
    gather("filename", "into_f", -feature_id)

  # Join features with peaks
  out <- temp_features %>%
    select(mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, peakidx, feature_id, any_of(c("isotopes", "adduct", "pcgroup"))) %>%
    rename_with(~paste0("f_", .), .cols = any_of(c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax"))) %>%
    left_join(temp_peaks, by = "peakidx") %>%
    mutate(filename = basename(filepath)) %>%
    left_join(temp_peaks_tab, by = c("feature_id", "filename"))

  # Filter to match peak intensities with feature intensities
  out <- out %>%
    filter(into == into_f) %>%
    select(-into_f) %>%
    group_by(feature_id, filename) %>%
    slice(1) %>%  # Take first if duplicates exist
    ungroup()

  # Complete to include all feature-sample combinations
  if (!is.null(xsAnnotate)) {
    out <- out %>%
      complete(nesting(feature_id, f_mzmed, f_mzmin, f_mzmax, f_rtmed, f_rtmin, f_rtmax,
                       isotopes, adduct, pcgroup, ms_level),
               nesting(filepath, filename, fromFile))
  } else {
    out <- out %>%
      complete(nesting(feature_id, f_mzmed, f_mzmin, f_mzmax, f_rtmed, f_rtmin, f_rtmax, ms_level),
               nesting(filepath, filename, fromFile))
  }

  # Add sample metadata - handle both XCMSnExp and XcmsExperiment
  if (inherits(XCMSnExp, "XcmsExperiment")) {
    # For XcmsExperiment, use sampleData()
    sample_metadata <- as.data.frame(sampleData(XCMSnExp)) %>%
      mutate(fromFile = 1:n()) %>%
      as_tibble()
  } else {
    # For XCMSnExp, use pData()
    sample_metadata <- pData(XCMSnExp) %>%
      mutate(fromFile = 1:n()) %>%
      as_tibble()
  }

  out <- sample_metadata %>%
    right_join(out, by = "fromFile", multiple = "all")

  return(out)
}
