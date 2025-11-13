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
