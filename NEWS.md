## Changes in v0.99.20 (commit: c7c96ef)


### Bug Fixes

* add missing comma in Remotes and add S4Vectors to Suggests ([974a3a8](https://github.com/stanstrup/tidyXCMS/commit/974a3a8344ebd6ce73de67cdf6c7c822fee9859a))
* add missing rlang import for syms function ([b0c0a58](https://github.com/stanstrup/tidyXCMS/commit/b0c0a58bebb497af837c63340163f97fb0c9d65f))
* apply feature_group fix before CAMERA annotation in test ([c6a09e3](https://github.com/stanstrup/tidyXCMS/commit/c6a09e302c95b4f0a904a004e7ca51c88a2ed6c0))
* correct fileNames documentation link ([018a7f6](https://github.com/stanstrup/tidyXCMS/commit/018a7f6465b2f11c8e40f5d50deec2837c16eb02))
* correct sampleData assignment and add alt text to vignette plots ([aa576af](https://github.com/stanstrup/tidyXCMS/commit/aa576af90ea253439eef412393759bbc4660c5a9))
* ensure git tags are pushed to remote during release ([8c9014b](https://github.com/stanstrup/tidyXCMS/commit/8c9014b9d9370ae1a69577c774161fe4b5753bd8))
* exclude .bioc_version file from package build ([38da0f1](https://github.com/stanstrup/tidyXCMS/commit/38da0f1ba19b5249bd9df05e6844d2c90e0c7189))
* handle both XCMSnExp and XcmsExperiment in pData test ([6251bc9](https://github.com/stanstrup/tidyXCMS/commit/6251bc9ce415af0640346c09519ba85020072c55))
* improve path matching in vignette to use basename ([80a3aa1](https://github.com/stanstrup/tidyXCMS/commit/80a3aa1b48d6f46c761b516696523b79c89db235))
* prevent semantic-release from creating 1.0.x tags and releases ([5863644](https://github.com/stanstrup/tidyXCMS/commit/586364490d188501ab997a1c147594ec99f9fff8))
* remove sample_index from documentation ([a63298c](https://github.com/stanstrup/tidyXCMS/commit/a63298c4c0e214099b3459d074ea16bbd84e5cf2))
* resolve NEWS.md duplicate versions and sync with GitHub releases ([90f9180](https://github.com/stanstrup/tidyXCMS/commit/90f91809628bec07e59417ebc81cba4da0a5eac7))
* resolve R CMD check warnings and test failures ([238aaea](https://github.com/stanstrup/tidyXCMS/commit/238aaea68e9de307dfeea8a2ea5bc8254b159dd6))
* resolve semantic-release sync-bioc-tag.sh execution failure ([2609ea2](https://github.com/stanstrup/tidyXCMS/commit/2609ea2231b4a372aae0dc4e55293b1676138e65))
* use list.files with recursive search for CDF path matching ([e85f6f5](https://github.com/stanstrup/tidyXCMS/commit/e85f6f5fd3ec246a17e9de7c28c92accb87f4580))


### Code Refactoring

* comprehensive code improvements and function rename to tidy_peaklist ([079448c](https://github.com/stanstrup/tidyXCMS/commit/079448caa0bebf7d1ac4eafc623d708a8ee2db92))


### Features

* add commit SHAs to NEWS.md for traceability ([07a41ce](https://github.com/stanstrup/tidyXCMS/commit/07a41cec272fb97b3ce3d0256f4173e269fdd590))
* add groupFeatures integration from MsFeatures package ([bea8117](https://github.com/stanstrup/tidyXCMS/commit/bea81177149f638e608c2b21411621fa13e0b776))


### BREAKING CHANGES

* Function renamed from XCMSnExp_CAMERA_peaklist_long to tidy_peaklist. Users must update function calls in existing code.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>

## Changes in v0.99.19 (commit: 624a737)


### Bug Fixes

* add .Rbuildignore to suppress hidden files warning ([11f4586](https://github.com/stanstrup/tidyXCMS/commit/11f458669ad6a6c1d59ab01a22c23882d7dff14c))
* add explicit MSnbase library call in vignette ([6f3c6fb](https://github.com/stanstrup/tidyXCMS/commit/6f3c6fb708ba67c24b893dc71f9b92b9cf0a1f59))
* add explicit namespace references for fileNames in vignette ([7d1e9f9](https://github.com/stanstrup/tidyXCMS/commit/7d1e9f99f6188f790f73561434b9417e7bbd0868))
* add GitHub remote for commonMZ package ([5401ce0](https://github.com/stanstrup/tidyXCMS/commit/5401ce0d4f70441f14051b9bf6a146795dc0a116))
* add missing comma in Remotes and add S4Vectors to Suggests ([974a3a8](https://github.com/stanstrup/tidyXCMS/commit/974a3a8344ebd6ce73de67cdf6c7c822fee9859a))
* add missing rlang import for syms function ([b0c0a58](https://github.com/stanstrup/tidyXCMS/commit/b0c0a58bebb497af837c63340163f97fb0c9d65f))
* apply feature_group fix before CAMERA annotation in test ([c6a09e3](https://github.com/stanstrup/tidyXCMS/commit/c6a09e302c95b4f0a904a004e7ca51c88a2ed6c0))
* convert prepare-news.sh to Unix line endings (LF) ([8abd768](https://github.com/stanstrup/tidyXCMS/commit/8abd768bf177f8d5a5f88f0268ccb7e6b35deb9b))
* correct fileNames documentation link ([018a7f6](https://github.com/stanstrup/tidyXCMS/commit/018a7f6465b2f11c8e40f5d50deec2837c16eb02))
* correct sampleData assignment and add alt text to vignette plots ([aa576af](https://github.com/stanstrup/tidyXCMS/commit/aa576af90ea253439eef412393759bbc4660c5a9))
* correct sampleData import from MsExperiment package ([2f0be73](https://github.com/stanstrup/tidyXCMS/commit/2f0be73cb040242ddfdd308e3adbbe48ab1882a7))
* correct xdata file paths to use system-specific location ([8d60f1a](https://github.com/stanstrup/tidyXCMS/commit/8d60f1ab93fb37b09579cb34157425e51f942dc4))
* ensure CAMERA annotation columns always exist in output ([4baab46](https://github.com/stanstrup/tidyXCMS/commit/4baab46ff4eeeaa0eb785f1a04fd1e100bba64a3))
* ensure git tags are pushed to remote during release ([8c9014b](https://github.com/stanstrup/tidyXCMS/commit/8c9014b9d9370ae1a69577c774161fe4b5753bd8))
* ensure package version stays below 1.0.0 for Bioconductor ([b1bc012](https://github.com/stanstrup/tidyXCMS/commit/b1bc012de886220414ab6302bd3ecafbef71eade))
* exclude .bioc_version file from package build ([38da0f1](https://github.com/stanstrup/tidyXCMS/commit/38da0f1ba19b5249bd9df05e6844d2c90e0c7189))
* handle both XCMSnExp and XcmsExperiment in pData test ([6251bc9](https://github.com/stanstrup/tidyXCMS/commit/6251bc9ce415af0640346c09519ba85020072c55))
* improve path matching in vignette to use basename ([80a3aa1](https://github.com/stanstrup/tidyXCMS/commit/80a3aa1b48d6f46c761b516696523b79c89db235))
* merge conflict ([d5aa8dc](https://github.com/stanstrup/tidyXCMS/commit/d5aa8dc3810f908492ad10d47273039648d9c98a))
* remove sample_index from documentation ([a63298c](https://github.com/stanstrup/tidyXCMS/commit/a63298c4c0e214099b3459d074ea16bbd84e5cf2))
* replace non-existent fileNames<- with proper dataOrigin replacement ([3995b38](https://github.com/stanstrup/tidyXCMS/commit/3995b38641c161957872d9ad48a80e57f2dbeec6))
* resolve devtools::check() test failures with loadXcmsData ([28f60b0](https://github.com/stanstrup/tidyXCMS/commit/28f60b06c9830fc08810dab113552fc23043bf99))
* resolve import issues and vignette build errors ([0a18dec](https://github.com/stanstrup/tidyXCMS/commit/0a18dec8ec30a93bfbc9176e5ca204af328ac209))
* resolve NEWS.md duplicate versions and sync with GitHub releases ([90f9180](https://github.com/stanstrup/tidyXCMS/commit/90f91809628bec07e59417ebc81cba4da0a5eac7))
* resolve package issues and improve CAMERA integration ([cc7a73c](https://github.com/stanstrup/tidyXCMS/commit/cc7a73c25f5b3e00e28023be69425a78a2b299d6))
* resolve R CMD check warnings and test failures ([238aaea](https://github.com/stanstrup/tidyXCMS/commit/238aaea68e9de307dfeea8a2ea5bc8254b159dd6))
* resolve semantic-release sync-bioc-tag.sh execution failure ([2609ea2](https://github.com/stanstrup/tidyXCMS/commit/2609ea2231b4a372aae0dc4e55293b1676138e65))
* resolve XcmsExperiment compatibility and CAMERA column issues ([a2e50ee](https://github.com/stanstrup/tidyXCMS/commit/a2e50ee39f98e99379047a206f20019b23e84552))
* skip groupCorr step in vignette for small example dataset ([561fa6a](https://github.com/stanstrup/tidyXCMS/commit/561fa6ac46b1330f004451d8b0746e7d51cb91e3))
* update author information with correct ORCID ([0af1279](https://github.com/stanstrup/tidyXCMS/commit/0af127966ab813879826dca91bc09f91c1fedfd1))
* use list.files with recursive search for CDF path matching ([e85f6f5](https://github.com/stanstrup/tidyXCMS/commit/e85f6f5fd3ec246a17e9de7c28c92accb87f4580))
* use loadXcmsData() instead of faahKO package ([70c89b0](https://github.com/stanstrup/tidyXCMS/commit/70c89b04150f1f927834253582ddae67413c2df6))


### Code Refactoring

* comprehensive code improvements and function rename to tidy_peaklist ([079448c](https://github.com/stanstrup/tidyXCMS/commit/079448caa0bebf7d1ac4eafc623d708a8ee2db92))


### Features

* add commit SHAs to NEWS.md for traceability ([07a41ce](https://github.com/stanstrup/tidyXCMS/commit/07a41cec272fb97b3ce3d0256f4173e269fdd590))
* add groupFeatures integration from MsFeatures package ([bea8117](https://github.com/stanstrup/tidyXCMS/commit/bea81177149f638e608c2b21411621fa13e0b776))
* add XCMSnExp_CAMERA_peaklist_long function and comprehensive package structure ([c0a7dc8](https://github.com/stanstrup/tidyXCMS/commit/c0a7dc89f2f770e6bfc4427c20dcd58137228775))
* integrate commonMZ for CAMERA adduct/fragment rules ([ababed5](https://github.com/stanstrup/tidyXCMS/commit/ababed546a07201cc0782c72bcecc03d163400b8))
* make CAMERA annotations optional with XcmsExperiment support ([dc802ab](https://github.com/stanstrup/tidyXCMS/commit/dc802abd45298110d9a9e07a2ea7304179c4859d))


### BREAKING CHANGES

* Function renamed from XCMSnExp_CAMERA_peaklist_long to tidy_peaklist. Users must update function calls in existing code.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>

## Changes in v0.99.18 (commit: c3916f6)


### Bug Fixes

* ensure git tags are pushed to remote during release ([8c9014b](https://github.com/stanstrup/tidyXCMS/commit/8c9014b9d9370ae1a69577c774161fe4b5753bd8))


### Features

* add commit SHAs to NEWS.md for traceability ([07a41ce](https://github.com/stanstrup/tidyXCMS/commit/07a41cec272fb97b3ce3d0256f4173e269fdd590))

## Changes in v0.99.17 (commit: d9352af)


### Bug Fixes

* correct fileNames documentation link ([018a7f6](https://github.com/stanstrup/tidyXCMS/commit/018a7f6465b2f11c8e40f5d50deec2837c16eb02))
* remove sample_index from documentation ([a63298c](https://github.com/stanstrup/tidyXCMS/commit/a63298c4c0e214099b3459d074ea16bbd84e5cf2))

## Changes in v0.99.16 (commit: ac33a5d)


### Bug Fixes

* resolve R CMD check warnings and test failures ([238aaea](https://github.com/stanstrup/tidyXCMS/commit/238aaea68e9de307dfeea8a2ea5bc8254b159dd6))

## Changes in v0.99.15 (commit: b839a77)


### Code Refactoring

* comprehensive code improvements and function rename to tidy_peaklist ([079448c](https://github.com/stanstrup/tidyXCMS/commit/079448caa0bebf7d1ac4eafc623d708a8ee2db92))


### BREAKING CHANGES

* Function renamed from XCMSnExp_CAMERA_peaklist_long to tidy_peaklist. Users must update function calls in existing code.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>

## Changes in v0.99.14 (commit: 76e5916)


### Bug Fixes

* use list.files with recursive search for CDF path matching ([e85f6f5](https://github.com/stanstrup/tidyXCMS/commit/e85f6f5fd3ec246a17e9de7c28c92accb87f4580))

## Changes in v0.99.13 (commit: ee81fef)


### Bug Fixes

* improve path matching in vignette to use basename ([80a3aa1](https://github.com/stanstrup/tidyXCMS/commit/80a3aa1b48d6f46c761b516696523b79c89db235))
* resolve semantic-release sync-bioc-tag.sh execution failure ([2609ea2](https://github.com/stanstrup/tidyXCMS/commit/2609ea2231b4a372aae0dc4e55293b1676138e65))

## Changes in v0.99.12 (commit: be6976b)


### Bug Fixes

* add missing comma in Remotes and add S4Vectors to Suggests ([974a3a8](https://github.com/stanstrup/tidyXCMS/commit/974a3a8344ebd6ce73de67cdf6c7c822fee9859a))

## Changes in v0.99.11 (commit: d82953f)


### Bug Fixes

* correct sampleData assignment and add alt text to vignette plots ([aa576af](https://github.com/stanstrup/tidyXCMS/commit/aa576af90ea253439eef412393759bbc4660c5a9))

## Changes in v0.99.10 (commit: 2f68705)


### Bug Fixes

* add missing rlang import for syms function ([b0c0a58](https://github.com/stanstrup/tidyXCMS/commit/b0c0a58bebb497af837c63340163f97fb0c9d65f))
* apply feature_group fix before CAMERA annotation in test ([c6a09e3](https://github.com/stanstrup/tidyXCMS/commit/c6a09e302c95b4f0a904a004e7ca51c88a2ed6c0))
* exclude .bioc_version file from package build ([38da0f1](https://github.com/stanstrup/tidyXCMS/commit/38da0f1ba19b5249bd9df05e6844d2c90e0c7189))
* handle both XCMSnExp and XcmsExperiment in pData test ([6251bc9](https://github.com/stanstrup/tidyXCMS/commit/6251bc9ce415af0640346c09519ba85020072c55))

## Changes in v0.99.9 (commit: d105433)


### Features

* add groupFeatures integration from MsFeatures package ([bea8117](https://github.com/stanstrup/tidyXCMS/commit/bea81177149f638e608c2b21411621fa13e0b776))

## Changes in v0.99.8 (commit: 0e05abc)


### Bug Fixes

* add .Rbuildignore to suppress hidden files warning ([11f4586](https://github.com/stanstrup/tidyXCMS/commit/11f458669ad6a6c1d59ab01a22c23882d7dff14c))
* add explicit MSnbase library call in vignette ([6f3c6fb](https://github.com/stanstrup/tidyXCMS/commit/6f3c6fb708ba67c24b893dc71f9b92b9cf0a1f59))
* add explicit namespace references for fileNames in vignette ([7d1e9f9](https://github.com/stanstrup/tidyXCMS/commit/7d1e9f99f6188f790f73561434b9417e7bbd0868))
* add GitHub remote for commonMZ package ([5401ce0](https://github.com/stanstrup/tidyXCMS/commit/5401ce0d4f70441f14051b9bf6a146795dc0a116))
* convert prepare-news.sh to Unix line endings (LF) ([8abd768](https://github.com/stanstrup/tidyXCMS/commit/8abd768bf177f8d5a5f88f0268ccb7e6b35deb9b))
* correct sampleData import from MsExperiment package ([2f0be73](https://github.com/stanstrup/tidyXCMS/commit/2f0be73cb040242ddfdd308e3adbbe48ab1882a7))
* correct xdata file paths to use system-specific location ([8d60f1a](https://github.com/stanstrup/tidyXCMS/commit/8d60f1ab93fb37b09579cb34157425e51f942dc4))
* ensure CAMERA annotation columns always exist in output ([4baab46](https://github.com/stanstrup/tidyXCMS/commit/4baab46ff4eeeaa0eb785f1a04fd1e100bba64a3))
* ensure package version stays below 1.0.0 for Bioconductor ([b1bc012](https://github.com/stanstrup/tidyXCMS/commit/b1bc012de886220414ab6302bd3ecafbef71eade))
* merge conflict ([d5aa8dc](https://github.com/stanstrup/tidyXCMS/commit/d5aa8dc3810f908492ad10d47273039648d9c98a))
* replace non-existent fileNames<- with proper dataOrigin replacement ([3995b38](https://github.com/stanstrup/tidyXCMS/commit/3995b38641c161957872d9ad48a80e57f2dbeec6))
* resolve devtools::check() test failures with loadXcmsData ([28f60b0](https://github.com/stanstrup/tidyXCMS/commit/28f60b06c9830fc08810dab113552fc23043bf99))
* resolve import issues and vignette build errors ([0a18dec](https://github.com/stanstrup/tidyXCMS/commit/0a18dec8ec30a93bfbc9176e5ca204af328ac209))
* resolve NEWS.md duplicate versions and sync with GitHub releases ([90f9180](https://github.com/stanstrup/tidyXCMS/commit/90f91809628bec07e59417ebc81cba4da0a5eac7))
* resolve package issues and improve CAMERA integration ([cc7a73c](https://github.com/stanstrup/tidyXCMS/commit/cc7a73c25f5b3e00e28023be69425a78a2b299d6))
* resolve XcmsExperiment compatibility and CAMERA column issues ([a2e50ee](https://github.com/stanstrup/tidyXCMS/commit/a2e50ee39f98e99379047a206f20019b23e84552))
* skip groupCorr step in vignette for small example dataset ([561fa6a](https://github.com/stanstrup/tidyXCMS/commit/561fa6ac46b1330f004451d8b0746e7d51cb91e3))
* update author information with correct ORCID ([0af1279](https://github.com/stanstrup/tidyXCMS/commit/0af127966ab813879826dca91bc09f91c1fedfd1))
* use loadXcmsData() instead of faahKO package ([70c89b0](https://github.com/stanstrup/tidyXCMS/commit/70c89b04150f1f927834253582ddae67413c2df6))


### Features

* add XCMSnExp_CAMERA_peaklist_long function and comprehensive package structure ([c0a7dc8](https://github.com/stanstrup/tidyXCMS/commit/c0a7dc89f2f770e6bfc4427c20dcd58137228775))
* integrate commonMZ for CAMERA adduct/fragment rules ([ababed5](https://github.com/stanstrup/tidyXCMS/commit/ababed546a07201cc0782c72bcecc03d163400b8))
* make CAMERA annotations optional with XcmsExperiment support ([dc802ab](https://github.com/stanstrup/tidyXCMS/commit/dc802abd45298110d9a9e07a2ea7304179c4859d))

## Changes in v0.99.7 (commit: 08d950f)


### Bug Fixes

* replace non-existent fileNames<- with proper dataOrigin replacement ([3995b38](https://github.com/stanstrup/tidyXCMS/commit/3995b38641c161957872d9ad48a80e57f2dbeec6))

## Changes in v0.99.6 (commit: 42429bf)


### Bug Fixes

* add explicit MSnbase library call in vignette ([6f3c6fb](https://github.com/stanstrup/tidyXCMS/commit/6f3c6fb708ba67c24b893dc71f9b92b9cf0a1f59))
* add explicit namespace references for fileNames in vignette ([7d1e9f9](https://github.com/stanstrup/tidyXCMS/commit/7d1e9f99f6188f790f73561434b9417e7bbd0868))

## Changes in v0.99.5 (commit: f952c46)


### Bug Fixes

* add GitHub remote for commonMZ package ([5401ce0](https://github.com/stanstrup/tidyXCMS/commit/5401ce0d4f70441f14051b9bf6a146795dc0a116))

## Changes in v0.99.4 (commit: 41f32fd)


### Features

* integrate commonMZ for CAMERA adduct/fragment rules ([ababed5](https://github.com/stanstrup/tidyXCMS/commit/ababed546a07201cc0782c72bcecc03d163400b8))

## Changes in v0.99.3 (commit: 729c750)


### Bug Fixes

* correct xdata file paths to use system-specific location ([8d60f1a](https://github.com/stanstrup/tidyXCMS/commit/8d60f1ab93fb37b09579cb34157425e51f942dc4))

## Changes in v0.99.2 (commit: cfed50a)


### Bug Fixes

* ensure package version stays below 1.0.0 for Bioconductor ([b1bc012](https://github.com/stanstrup/tidyXCMS/commit/b1bc012de886220414ab6302bd3ecafbef71eade))

## Changes in v0.99.1 (commit: b1ef759)


### Bug Fixes

* add .Rbuildignore to suppress hidden files warning ([11f4586](https://github.com/stanstrup/tidyXCMS/commit/11f458669ad6a6c1d59ab01a22c23882d7dff14c))
* convert prepare-news.sh to Unix line endings (LF) ([8abd768](https://github.com/stanstrup/tidyXCMS/commit/8abd768bf177f8d5a5f88f0268ccb7e6b35deb9b))
* correct sampleData import from MsExperiment package ([2f0be73](https://github.com/stanstrup/tidyXCMS/commit/2f0be73cb040242ddfdd308e3adbbe48ab1882a7))
* ensure CAMERA annotation columns always exist in output ([4baab46](https://github.com/stanstrup/tidyXCMS/commit/4baab46ff4eeeaa0eb785f1a04fd1e100bba64a3))
* merge conflict ([d5aa8dc](https://github.com/stanstrup/tidyXCMS/commit/d5aa8dc3810f908492ad10d47273039648d9c98a))
* resolve devtools::check() test failures with loadXcmsData ([28f60b0](https://github.com/stanstrup/tidyXCMS/commit/28f60b06c9830fc08810dab113552fc23043bf99))
* resolve import issues and vignette build errors ([0a18dec](https://github.com/stanstrup/tidyXCMS/commit/0a18dec8ec30a93bfbc9176e5ca204af328ac209))
* resolve package issues and improve CAMERA integration ([cc7a73c](https://github.com/stanstrup/tidyXCMS/commit/cc7a73c25f5b3e00e28023be69425a78a2b299d6))
* resolve XcmsExperiment compatibility and CAMERA column issues ([a2e50ee](https://github.com/stanstrup/tidyXCMS/commit/a2e50ee39f98e99379047a206f20019b23e84552))
* skip groupCorr step in vignette for small example dataset ([561fa6a](https://github.com/stanstrup/tidyXCMS/commit/561fa6ac46b1330f004451d8b0746e7d51cb91e3))
* update author information with correct ORCID ([0af1279](https://github.com/stanstrup/tidyXCMS/commit/0af127966ab813879826dca91bc09f91c1fedfd1))
* use loadXcmsData() instead of faahKO package ([70c89b0](https://github.com/stanstrup/tidyXCMS/commit/70c89b04150f1f927834253582ddae67413c2df6))


### Features

* add XCMSnExp_CAMERA_peaklist_long function and comprehensive package structure ([c0a7dc8](https://github.com/stanstrup/tidyXCMS/commit/c0a7dc89f2f770e6bfc4427c20dcd58137228775))
* make CAMERA annotations optional with XcmsExperiment support ([dc802ab](https://github.com/stanstrup/tidyXCMS/commit/dc802abd45298110d9a9e07a2ea7304179c4859d))
