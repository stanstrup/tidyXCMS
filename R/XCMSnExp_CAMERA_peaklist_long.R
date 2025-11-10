#' Create Long-Format Peak Table with CAMERA Annotations
#'
#' Generates a comprehensive long-format peak table where each row represents
#' one feature for one sample, including CAMERA annotations for isotopes,
#' adducts, and pseudospectrum groups.
#'
#' @param XCMSnExp An [xcms::XCMSnExp-class] object containing peak detection
#'   and feature grouping results.
#' @param xsAnnotate An `xsAnnotate` object from the CAMERA package containing
#'   peak annotations (isotopes, adducts, pseudospectrum groups).
#'
#' @return A [tibble::tibble] in long format with one row per feature per sample.
#'   The tibble contains:
#'   \describe{
#'     \item{Feature-level columns}{f_mzmed, f_mzmin, f_mzmax, f_rtmed, f_rtmin,
#'       f_rtmax - feature summary statistics}
#'     \item{CAMERA annotations}{isotopes, adduct, pcgroup (pseudospectrum group)}
#'     \item{Peak-level columns}{mz, rt, into, intb, maxo, sn - individual peak
#'       measurements}
#'     \item{Sample information}{filepath, filename, fromFile - sample identifiers}
#'     \item{Additional columns}{Any columns from pData(XCMSnExp)}
#'   }
#'
#' @details
#' This function integrates data from multiple sources:
#' \itemize{
#'   \item Chromatographic peaks from [xcms::chromPeaks()]
#'   \item Feature definitions from [xcms::featureDefinitions()]
#'   \item Feature intensities from [xcms::featureValues()]
#'   \item CAMERA annotations from the xsAnnotate object
#'   \item Sample metadata from [Biobase::pData()]
#' }
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
#'   \item CAMERA must be run prior to using this function to obtain the
#'     xsAnnotate object.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(xcms)
#' library(CAMERA)
#' library(BiocParallel)
#' library(faahKO)
#'
#' # Load example data
#' data("faahko_sub", package = "faahKO")
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
#' # CAMERA annotation
#' xs <- xsAnnotate(xdata)
#' xs <- groupFWHM(xs)
#' xs <- findIsotopes(xs)
#' xs <- groupCorr(xs)
#' xs <- findAdducts(xs)
#'
#' # Create long-format peak table
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
#' }
#'
#' @importFrom xcms chromPeaks chromPeakData featureDefinitions featureValues fileNames
#' @importFrom dplyr %>% mutate rename left_join right_join filter group_by
#'   ungroup slice select rename_at vars bind_cols as_tibble
#' @importFrom tidyr unnest gather complete nesting
#' @importFrom tibble as_tibble
#' @importFrom purrr one_of
#' @importFrom Biobase pData
XCMSnExp_CAMERA_peaklist_long <- function(XCMSnExp, xsAnnotate) {
  # Extract peaks with file information
  temp_peaks <- chromPeaks_classic(XCMSnExp) %>%
    mutate(peakidx = 1:nrow(.)) %>%
    rename(fromFile = sample) %>%
    mutate(filepath = fileNames(XCMSnExp)[fromFile])

  # Extract features and combine with CAMERA annotations
  temp_features <- featureDefinitions(XCMSnExp) %>%
    as_tibble() %>%
    bind_cols(as_tibble(CAMERA::getPeaklist(xsAnnotate))[, c("isotopes", "adduct", "pcgroup")]) %>%
    mutate(pcgroup = as.integer(pcgroup)) %>%
    mutate(feature_id = 1:nrow(.)) %>%
    unnest(peakidx)

  # Extract feature values (intensities)
  temp_peaks_tab <- featureValues(XCMSnExp, value = "into", method = "maxint") %>%
    as_tibble() %>%
    mutate(feature_id = 1:nrow(.)) %>%
    gather("filename", "into_f", -feature_id)

  # Join features with peaks
  out <- temp_features %>%
    select(mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, peakidx, feature_id, isotopes, adduct, pcgroup) %>%
    rename_at(vars(one_of("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax")),
              list(~paste0("f_", .))) %>%
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
  out <- out %>%
    complete(nesting(feature_id, f_mzmed, f_mzmin, f_mzmax, f_rtmed, f_rtmin, f_rtmax,
                     isotopes, adduct, pcgroup, ms_level),
             nesting(filepath, filename, fromFile))

  # Add sample metadata
  out <- pData(XCMSnExp) %>%
    mutate(fromFile = 1:n()) %>%
    right_join(out, by = "fromFile", multiple = "all")

  return(out)
}
