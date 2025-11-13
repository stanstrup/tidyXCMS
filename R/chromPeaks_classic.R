#' Extract Chromatographic Peaks as Tibble (Internal)
#'
#' Extracts chromatographic peaks from an XCMSnExp or XcmsExperiment object
#' and returns them as a tibble, combining peak matrix data with peak data.
#'
#' @param x An [xcms::XCMSnExp-class] or [xcms::XcmsExperiment-class] object
#'   containing chromatographic peak data.
#'
#' @return A [tibble::tibble] with one row per chromatographic peak. Contains
#'   all columns from [xcms::chromPeaks()] and [xcms::chromPeakData()].
#'
#' @details
#' This function provides a convenient way to extract peak information in a
#' tidy format. It combines the peak matrix (mz, rt, intensity values, etc.)
#' with additional peak annotations stored in chromPeakData.
#'
#' **Historical note**: In older versions of xcms, there was a single function
#' that returned all peak information together. In modern xcms, peak data is
#' split between chromPeaks() (the numeric matrix with m/z, RT, intensities)
#' and chromPeakData() (a DataFrame with additional annotations). This helper
#' function merges both sources to recreate the classic unified format as a
#' tibble.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom xcms chromPeaks chromPeakData
#' @importFrom dplyr %>% bind_cols
#' @importFrom tibble as_tibble
.chromPeaks_classic <- function(x) {
  chromPeaks(x) %>%
    as.data.frame() %>%
    bind_cols(as.data.frame(chromPeakData(x))) %>%
    as_tibble()
}
