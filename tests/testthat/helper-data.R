# Helper file to create test data for testthat tests
# This file is automatically loaded before tests run

library(xcms)
library(BiocParallel)

# Load faahko_sub test data
# This loads the pre-processed data object from xcms package
# and manually fixes file paths to avoid loadXcmsData() issues during R CMD check
create_faahko_sub <- function() {
  # Load the data object from xcms package
  # This is a pre-processed XCMSnExp object with detected peaks
  e <- new.env()
  data("faahko_sub", envir = e, package = "xcms")
  obj <- get("faahko_sub", e)

  # Get the actual file paths from faahKO package
  cdf_path <- system.file("cdf/KO", package = "faahKO")

  # Check if package is available
  if (cdf_path == "") {
    stop("Package 'faahKO' not available. Please install using BiocManager::install('faahKO')")
  }

  # Get the files that match the basenames in the object
  file_basenames <- basename(fileNames(obj))
  files <- file.path(cdf_path, file_basenames)

  # Verify files exist
  if (!all(file.exists(files))) {
    stop("Some CDF files not found: ", paste(files[!file.exists(files)], collapse = ", "))
  }

  # Update file paths in the object using fromFile slot
  # This avoids the problematic dirname<- approach used by loadXcmsData
  obj@processingData@files <- files

  return(obj)
}

# Create the test data object
# This will be available to all tests
faahko_sub <- create_faahko_sub()
