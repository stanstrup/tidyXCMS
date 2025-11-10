# Test suite for XCMSnExp_CAMERA_peaklist_long and .chromPeaks_classic functions

test_that(".chromPeaks_classic returns a tibble with expected structure", {

  library(xcms)
  library(BiocParallel)

  # Test .chromPeaks_classic
  peaks_tbl <- tidyXCMS:::.chromPeaks_classic(faahko_sub)

  # Check it's a tibble
  expect_s3_class(peaks_tbl, "tbl_df")

  # Check it has rows
  expect_gt(nrow(peaks_tbl), 0)

  # Check for expected columns from chromPeaks
  expect_true("mz" %in% colnames(peaks_tbl))
  expect_true("rt" %in% colnames(peaks_tbl))
  expect_true("into" %in% colnames(peaks_tbl))
  expect_true("sample" %in% colnames(peaks_tbl))

  # Check number of rows matches chromPeaks
  expect_equal(nrow(peaks_tbl), nrow(chromPeaks(faahko_sub)))
})


test_that(".chromPeaks_classic works with XcmsExperiment objects", {


  library(xcms)
  library(MsExperiment)
  library(BiocParallel)

  # Load example files
  fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)

  # Read as XcmsExperiment
  xdata <- readMsExperiment(spectraFiles = fls[1:2], BPPARAM = SerialParam())

  # Peak detection
  cwp <- CentWaveParam(peakwidth = c(5, 20), noise = 100)
  xdata <- findChromPeaks(xdata, param = cwp, BPPARAM = SerialParam())

  # Test .chromPeaks_classic
  peaks_tbl <- tidyXCMS:::.chromPeaks_classic(xdata)

  # Check it's a tibble
  expect_s3_class(peaks_tbl, "tbl_df")

  # Check it has rows
  expect_gt(nrow(peaks_tbl), 0)

  # Check for expected columns
  expect_true("mz" %in% colnames(peaks_tbl))
  expect_true("rt" %in% colnames(peaks_tbl))
})


test_that("XCMSnExp_CAMERA_peaklist_long returns expected structure", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)

  # Load example data
  faahko_sub <- loadXcmsData("faahko_sub")

  # Peak detection
  cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
  xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())

  # Peak grouping
  pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                          minFraction = 0.5)
  xdata <- groupChromPeaks(xdata, param = pdp)

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)
  xs <- findIsotopes(xs)
  xs <- groupCorr(xs)
  xs <- findAdducts(xs, polarity = "positive")

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata, xs)

  # Check it's a tibble
  expect_s3_class(peak_table, "tbl_df")

  # Check it has rows
  expect_gt(nrow(peak_table), 0)

  # Check for expected feature-level columns
  expect_true("feature_id" %in% colnames(peak_table))
  expect_true("f_mzmed" %in% colnames(peak_table))
  expect_true("f_rtmed" %in% colnames(peak_table))
  expect_true("f_mzmin" %in% colnames(peak_table))
  expect_true("f_mzmax" %in% colnames(peak_table))
  expect_true("f_rtmin" %in% colnames(peak_table))
  expect_true("f_rtmax" %in% colnames(peak_table))

  # Check for CAMERA annotation columns
  expect_true("isotopes" %in% colnames(peak_table))
  expect_true("adduct" %in% colnames(peak_table))
  expect_true("pcgroup" %in% colnames(peak_table))

  # Check for sample information columns
  expect_true("filename" %in% colnames(peak_table))
  expect_true("filepath" %in% colnames(peak_table))
  expect_true("fromFile" %in% colnames(peak_table))

  # Check for peak-level columns
  expect_true("mz" %in% colnames(peak_table))
  expect_true("rt" %in% colnames(peak_table))
  expect_true("into" %in% colnames(peak_table))
})


test_that("XCMSnExp_CAMERA_peaklist_long has one row per feature per sample", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)
  library(dplyr)

  # Load example data
  faahko_sub <- loadXcmsData("faahko_sub")

  # Peak detection
  cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
  xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())

  # Peak grouping
  pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                          minFraction = 0.5)
  xdata <- groupChromPeaks(xdata, param = pdp)

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)
  xs <- findIsotopes(xs)

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata, xs)

  # Get number of features and samples
  n_features <- nrow(featureDefinitions(xdata))
  n_samples <- length(fileNames(xdata))

  # Check total rows equals features * samples
  expect_equal(nrow(peak_table), n_features * n_samples)

  # Check no duplicate feature-sample combinations
  duplicates <- peak_table %>%
    group_by(feature_id, filename) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)

  expect_equal(nrow(duplicates), 0)
})


test_that("XCMSnExp_CAMERA_peaklist_long handles missing values correctly", {


  library(xcms)
  library(CAMERA)
  library(BiocParallel)
  library(dplyr)

  # Load example data
  faahko_sub <- loadXcmsData("faahko_sub")

  # Peak detection
  cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
  xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())

  # Peak grouping with minFraction < 1 to ensure some missing values
  pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                          minFraction = 0.3)
  xdata <- groupChromPeaks(xdata, param = pdp)

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata, xs)

  # Check that we have some NA values in intensity columns
  # (features not detected in all samples)
  expect_true(any(is.na(peak_table$into)))

  # Check that feature-level columns are never NA
  expect_false(any(is.na(peak_table$feature_id)))
  expect_false(any(is.na(peak_table$f_mzmed)))
  expect_false(any(is.na(peak_table$f_rtmed)))

  # Check that sample information is never NA
  expect_false(any(is.na(peak_table$filename)))
  expect_false(any(is.na(peak_table$fromFile)))
})


test_that("XCMSnExp_CAMERA_peaklist_long includes pData information", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)
  library(Biobase)

  # Load example data
  faahko_sub <- loadXcmsData("faahko_sub")

  # Add sample metadata
  pd <- data.frame(
    sample_name = basename(fileNames(faahko_sub)),
    sample_group = c("WT", "KO", "WT")[1:length(fileNames(faahko_sub))],
    row.names = basename(fileNames(faahko_sub))
  )
  pData(faahko_sub) <- pd

  # Peak detection
  cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
  xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())

  # Peak grouping
  pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                          minFraction = 0.5)
  xdata <- groupChromPeaks(xdata, param = pdp)

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata, xs)

  # Check that pData columns are present
  expect_true("sample_name" %in% colnames(peak_table))
  expect_true("sample_group" %in% colnames(peak_table))
})


test_that("XCMSnExp_CAMERA_peaklist_long CAMERA annotations are present", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)
  library(dplyr)

  # Load example data
  faahko_sub <- loadXcmsData("faahko_sub")

  # Peak detection
  cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
  xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())

  # Peak grouping
  pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                          minFraction = 0.5)
  xdata <- groupChromPeaks(xdata, param = pdp)

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation with all steps
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)
  xs <- findIsotopes(xs)
  xs <- groupCorr(xs)
  xs <- findAdducts(xs, polarity = "positive")

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata, xs)

  # Check that CAMERA annotations are present
  # Note: Not all features will have annotations, but columns should exist
  expect_true("isotopes" %in% colnames(peak_table))
  expect_true("adduct" %in% colnames(peak_table))
  expect_true("pcgroup" %in% colnames(peak_table))

  # pcgroup should be numeric/integer
  expect_true(is.numeric(peak_table$pcgroup))

  # Check that at least some features have pcgroup assignments
  # (after groupCorr, features should be assigned to pseudospectrum groups)
  non_na_pcgroups <- peak_table %>%
    filter(!is.na(feature_id)) %>%
    distinct(feature_id, pcgroup) %>%
    filter(!is.na(pcgroup))

  expect_gt(nrow(non_na_pcgroups), 0)
})


test_that("XCMSnExp_CAMERA_peaklist_long works without CAMERA annotations", {

  library(xcms)
  library(BiocParallel)

  # Load example data
  faahko_sub <- loadXcmsData("faahko_sub")

  # Peak detection
  cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
  xdata <- findChromPeaks(faahko_sub, param = cwp, BPPARAM = SerialParam())

  # Peak grouping
  pdp <- PeakDensityParam(sampleGroups = rep(1, length(fileNames(xdata))),
                          minFraction = 0.5)
  xdata <- groupChromPeaks(xdata, param = pdp)

  # Create long-format peak table WITHOUT CAMERA
  peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata)

  # Check it's a tibble
  expect_s3_class(peak_table, "tbl_df")

  # Check it has rows
  expect_gt(nrow(peak_table), 0)

  # Check for expected columns
  expect_true("feature_id" %in% colnames(peak_table))
  expect_true("f_mzmed" %in% colnames(peak_table))

  # Check that CAMERA columns are NOT present
  expect_false("isotopes" %in% colnames(peak_table))
  expect_false("adduct" %in% colnames(peak_table))
  expect_false("pcgroup" %in% colnames(peak_table))
})


test_that("XCMSnExp_CAMERA_peaklist_long works with XcmsExperiment", {

  library(xcms)
  library(MsExperiment)
  library(BiocParallel)

  # Load example files
  fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)

  # Read as XcmsExperiment
  xdata <- readMsExperiment(spectraFiles = fls[1:2], BPPARAM = SerialParam())

  # Peak detection
  cwp <- CentWaveParam(peakwidth = c(5, 20), noise = 100)
  xdata <- findChromPeaks(xdata, param = cwp, BPPARAM = SerialParam())

  # Peak grouping
  pdp <- PeakDensityParam(sampleGroups = rep(1, 2), minFraction = 0.5)
  xdata <- groupChromPeaks(xdata, param = pdp)

  # Create long-format peak table
  peak_table <- XCMSnExp_CAMERA_peaklist_long(xdata)

  # Check it's a tibble
  expect_s3_class(peak_table, "tbl_df")

  # Check it has rows
  expect_gt(nrow(peak_table), 0)

  # Check for expected columns
  expect_true("feature_id" %in% colnames(peak_table))
  expect_true("f_mzmed" %in% colnames(peak_table))
  expect_true("filename" %in% colnames(peak_table))
  expect_true("mz" %in% colnames(peak_table))
  expect_true("rt" %in% colnames(peak_table))

  # Check that CAMERA columns are NOT present
  expect_false("isotopes" %in% colnames(peak_table))
  expect_false("adduct" %in% colnames(peak_table))
  expect_false("pcgroup" %in% colnames(peak_table))
})
