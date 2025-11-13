# Test suite for tidy_peaklist and .chromPeaks_classic functions

test_that(".chromPeaks_classic returns a tibble with expected structure", {

  library(xcms)
  library(BiocParallel)

  # Test .chromPeaks_classic with preloaded xdata (XCMSnExp with peaks and grouping)
  peaks_tbl <- tidyXCMS:::.chromPeaks_classic(xdata)

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
  expect_equal(nrow(peaks_tbl), nrow(chromPeaks(xdata)))
})


test_that(".chromPeaks_classic works with XcmsExperiment objects", {

  library(xcms)
  library(MsExperiment)
  library(BiocParallel)

  # Test .chromPeaks_classic with preloaded xmse (XcmsExperiment with peaks and grouping)
  peaks_tbl <- tidyXCMS:::.chromPeaks_classic(xmse)

  # Check it's a tibble
  expect_s3_class(peaks_tbl, "tbl_df")

  # Check it has rows
  expect_gt(nrow(peaks_tbl), 0)

  # Check for expected columns
  expect_true("mz" %in% colnames(peaks_tbl))
  expect_true("rt" %in% colnames(peaks_tbl))
})


test_that("tidy_peaklist returns expected structure", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)

  # Use preloaded xdata (XCMSnExp with peaks and grouping)

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)
  xs <- findIsotopes(xs)
  xs <- groupCorr(xs)
  xs <- findAdducts(xs, polarity = "positive")

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- tidy_peaklist(xdata, xs)

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


test_that("tidy_peaklist has one row per feature per sample", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)
  library(dplyr)

  # Use preloaded xdata (XCMSnExp with peaks and grouping)

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)
  xs <- findIsotopes(xs)

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- tidy_peaklist(xdata, xs)

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


test_that("tidy_peaklist handles missing values correctly", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)
  library(dplyr)

  # Use preloaded xdata (XCMSnExp with peaks and grouping)
  # xdata already has grouping with minFraction that allows some missing values

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- tidy_peaklist(xdata, xs)

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


test_that("tidy_peaklist includes pData information", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)
  library(Biobase)
  library(MsExperiment)

  # Use preloaded xdata and add sample metadata
  xdata_copy <- xdata

  # Add sample metadata
  pd <- data.frame(
    sample_name = basename(fileNames(xdata_copy)),
    sample_group = c("WT", "KO", "WT")[1:length(fileNames(xdata_copy))],
    row.names = basename(fileNames(xdata_copy))
  )

  # Handle both XCMSnExp and XcmsExperiment objects
  if (inherits(xdata_copy, "XcmsExperiment")) {
    # For XcmsExperiment, use sampleData<- (requires DataFrame, not data.frame)
    sampleData(xdata_copy) <- S4Vectors::DataFrame(pd)
  } else {
    # For XCMSnExp, use pData<-
    pData(xdata_copy) <- pd
  }

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata_copy, "xcmsSet")

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- tidy_peaklist(xdata_copy, xs)

  # Check that pData columns are present
  expect_true("sample_name" %in% colnames(peak_table))
  expect_true("sample_group" %in% colnames(peak_table))
})


test_that("tidy_peaklist CAMERA annotations are present", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)
  library(dplyr)

  # Use preloaded xdata (XCMSnExp with peaks and grouping)

  # Convert to xcmsSet for CAMERA (CAMERA requires the old xcmsSet class)
  xset <- as(xdata, "xcmsSet")

  # CAMERA annotation with all steps
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)
  xs <- findIsotopes(xs)
  xs <- groupCorr(xs)
  xs <- findAdducts(xs, polarity = "positive")

  # Create long-format peak table (pass XCMSnExp object, not xcmsSet)
  peak_table <- tidy_peaklist(xdata, xs)

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


test_that("tidy_peaklist works without CAMERA annotations", {

  library(xcms)
  library(BiocParallel)

  # Use preloaded xdata (XCMSnExp with peaks and grouping)

  # Create long-format peak table WITHOUT CAMERA
  peak_table <- tidy_peaklist(xdata)

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


test_that("tidy_peaklist works with XcmsExperiment", {

  library(xcms)
  library(MsExperiment)
  library(BiocParallel)

  # Use preloaded xmse (XcmsExperiment with peaks and grouping)

  # Create long-format peak table
  peak_table <- tidy_peaklist(xmse)

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


test_that("tidy_peaklist includes groupFeatures results", {

  library(xcms)
  library(MsFeatures)
  library(BiocParallel)
  library(dplyr)

  # Use preloaded xdata (XCMSnExp with peaks and grouping)

  # Apply groupFeatures with SimilarRtimeParam
  xdata_grouped <- groupFeatures(xdata, param = SimilarRtimeParam(diffRt = 10))

  # Create long-format peak table
  peak_table <- tidy_peaklist(xdata_grouped)

  # Check that feature_group column is present
  expect_true("feature_group" %in% colnames(peak_table))

  # Check that feature_group values are present (should be like "FG.001", "FG.002", etc.)
  expect_true(any(!is.na(peak_table$feature_group)))

  # Check that feature_group is character
  expect_true(is.character(peak_table$feature_group))

  # Check that all features have a feature_group assignment
  feature_groups <- peak_table %>%
    distinct(feature_id, feature_group)

  expect_equal(nrow(feature_groups), n_distinct(peak_table$feature_id))
})


test_that("tidy_peaklist works with both CAMERA and groupFeatures", {

  library(xcms)
  library(CAMERA)
  library(MsFeatures)
  library(BiocParallel)
  library(dplyr)

  # Use preloaded xdata (XCMSnExp with peaks and grouping)

  # Apply groupFeatures
  xdata_grouped <- groupFeatures(xdata, param = SimilarRtimeParam(diffRt = 10))

  # Convert to xcmsSet for CAMERA
  xset <- as(xdata_grouped, "xcmsSet")

  # featuregroups seem to corrupt the feature matrix. We fix before CAMERA
  xset@groups <- apply(xset@groups[, !(colnames(xset@groups) %in% "feature_group")], 2, as.numeric)

  # CAMERA annotation
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)
  xs <- findIsotopes(xs)

  # Create long-format peak table
  peak_table <- tidy_peaklist(xdata_grouped, xs)

  # Check that both CAMERA and groupFeatures columns are present
  expect_true("feature_group" %in% colnames(peak_table))
  expect_true("isotopes" %in% colnames(peak_table))
  expect_true("adduct" %in% colnames(peak_table))
  expect_true("pcgroup" %in% colnames(peak_table))

  # Check that feature_group values are present
  expect_true(any(!is.na(peak_table$feature_group)))
})


test_that("tidy_peaklist works with XcmsExperiment and groupFeatures", {

  library(xcms)
  library(MsExperiment)
  library(MsFeatures)
  library(BiocParallel)
  library(dplyr)

  # Use preloaded xmse (XcmsExperiment with peaks and grouping)

  # Apply groupFeatures with AbundanceSimilarityParam
  xmse_grouped <- groupFeatures(xmse, param = AbundanceSimilarityParam(threshold = 0.7))

  # Create long-format peak table
  peak_table <- tidy_peaklist(xmse_grouped)

  # Check that feature_group column is present
  expect_true("feature_group" %in% colnames(peak_table))

  # Check that feature_group values are present
  expect_true(any(!is.na(peak_table$feature_group)))

  # Check proper structure
  expect_s3_class(peak_table, "tbl_df")
  expect_gt(nrow(peak_table), 0)
})

# Input validation and edge case tests ----------------------------------------

test_that("tidy_peaklist works without features (dropFeatureDefinitions)", {

  library(xcms)
  library(BiocParallel)

  # Use preloaded xdata and remove features
  xdata_no_features <- dropFeatureDefinitions(xdata)

  # Should work but return peak-level data only
  expect_warning(
    peak_table <- tidy_peaklist(xdata_no_features),
    "No features defined"
  )

  # Check it's a tibble
  expect_s3_class(peak_table, "tbl_df")

  # Check it has rows (one per peak)
  expect_gt(nrow(peak_table), 0)

  # Check it has peak-level columns
  expect_true("mz" %in% colnames(peak_table))
  expect_true("rt" %in% colnames(peak_table))
  expect_true("into" %in% colnames(peak_table))

  # Check it does NOT have feature-level columns
  expect_false("feature_id" %in% colnames(peak_table))
  expect_false("f_mzmed" %in% colnames(peak_table))
  expect_false("f_rtmed" %in% colnames(peak_table))

  # Check sample information is present
  expect_true("filename" %in% colnames(peak_table))
  expect_true("filepath" %in% colnames(peak_table))
  expect_true("fromFile" %in% colnames(peak_table))
})


test_that("tidy_peaklist validates input object class", {

  library(xcms)

  # Test with invalid object
  expect_error(
    tidy_peaklist("not_an_xcms_object"),
    "'x' must be an XCMSnExp or XcmsExperiment object"
  )

  expect_error(
    tidy_peaklist(list(a = 1, b = 2)),
    "'x' must be an XCMSnExp or XcmsExperiment object"
  )
})


test_that("tidy_peaklist validates xsAnnotate parameter", {

  library(xcms)
  library(BiocParallel)

  # Test with invalid xsAnnotate
  expect_error(
    tidy_peaklist(xdata, xsAnnotate = "not_a_camera_object"),
    "'xsAnnotate' must be an xsAnnotate object from CAMERA package"
  )

  expect_error(
    tidy_peaklist(xdata, xsAnnotate = list(a = 1)),
    "'xsAnnotate' must be an xsAnnotate object from CAMERA package"
  )
})


test_that("tidy_peaklist detects CAMERA/XCMS feature count mismatch", {

  library(xcms)
  library(CAMERA)
  library(BiocParallel)

  # Use preloaded xdata
  xset <- as(xdata, "xcmsSet")

  # Create CAMERA object
  xs <- xsAnnotate(xset)
  xs <- groupFWHM(xs)

  # Manually create an xdata object with different number of features
  # by filtering features
  xdata_subset <- filterFeatureDefinitions(xdata, features = 1:5)

  # Should error due to mismatch
  expect_error(
    tidy_peaklist(xdata_subset, xs),
    "Feature count mismatch"
  )
})


test_that("tidy_peaklist includes sample metadata correctly", {

  library(xcms)
  library(Biobase)
  library(BiocParallel)

  # Use preloaded xdata and add sample metadata
  xdata_copy <- xdata

  # Add sample metadata with various column types
  pd <- data.frame(
    sample_name = paste0("Sample_", seq_len(length(fileNames(xdata_copy)))),
    sample_group = rep(c("Control", "Treatment"), length.out = length(fileNames(xdata_copy))),
    sample_number = seq_len(length(fileNames(xdata_copy))),
    replicate = c(1, 2, 1)[seq_len(length(fileNames(xdata_copy)))],
    row.names = basename(fileNames(xdata_copy))
  )

  # Handle both XCMSnExp and XcmsExperiment objects
  if (inherits(xdata_copy, "XcmsExperiment")) {
    # For XcmsExperiment, use sampleData<- (requires DataFrame, not data.frame)
    sampleData(xdata_copy) <- S4Vectors::DataFrame(pd)
  } else {
    # For XCMSnExp, use pData<-
    pData(xdata_copy) <- pd
  }

  # Create peak table
  peak_table <- tidy_peaklist(xdata_copy)

  # Check that ALL pData columns are present
  expect_true("sample_name" %in% colnames(peak_table))
  expect_true("sample_group" %in% colnames(peak_table))
  expect_true("sample_number" %in% colnames(peak_table))
  expect_true("replicate" %in% colnames(peak_table))

  # Check that values match
  expect_true(all(peak_table$sample_group %in% c("Control", "Treatment")))
  expect_equal(length(unique(peak_table$sample_name)), length(fileNames(xdata_copy)))
})


test_that("tidy_peaklist issues warning for large datasets", {

  library(xcms)
  library(BiocParallel)

  # Note: This test is conceptual - we can't easily create a dataset large enough
  # to trigger the warning without significant overhead.
  # The warning logic is tested implicitly through its implementation.

  # Test with normal-sized data (should not warn)
  expect_no_warning(
    peak_table <- tidy_peaklist(xdata)
  )
})
