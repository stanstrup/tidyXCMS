# Helper file to create test data for testthat tests
# This file is automatically loaded before tests run

library(xcms)
library(BiocParallel)

# Load pre-processed test data objects
# These objects already have peak detection and grouping completed

# xdata: XCMSnExp object with peak detection and grouping
# This is more efficient than manual peak picking in tests
xdata <- loadXcmsData("xdata")

# xmse: XcmsExperiment object with peak detection and grouping
# Use this for testing XcmsExperiment compatibility
xmse <- loadXcmsData("xmse")


# CAMERA needs the legacy xcmsSet class. The conversion is not available on
# every platform/version combination (it fails with infinite recursion on
# R-devel with Bioconductor devel), so CAMERA tests are skipped when it fails.
camera_works <- requireNamespace("CAMERA", quietly = TRUE) &&
  tryCatch({
    methods::as(xdata, "xcmsSet")
    TRUE
  }, error = function(e) FALSE)

skip_if_no_camera <- function() {
  testthat::skip_if_not(camera_works,
                        "CAMERA or the conversion to xcmsSet is not available")
}
