## exprso 0.4.7
---------------------
* Update `exprso` import
    * Correctly apply `make.names` to feature names
* Improvements to `calcStats`
    * Add correlation and R-squared measures

## exprso 0.4.6
---------------------
* Update `mod` modules
    * New `modPermute` will permute predictor variables
* Update `getWeights` method
    * Extend support to `randomForest` objects

## exprso 0.4.5
---------------------
* Update `fs` modules
    * Extend `fsBalance` to support `sbp.fromPropd`

## exprso 0.4.4
---------------------
* Update `mod` modules
    * Extend `modHistory` to support balance deployment
* Update `fs` modules
    * New `fsBalance` to convert features into balances
    * Uses `sbp.subset` if `ternary` or `ratios` argument is provided
* Package maintenance
    * Turn README into FAQ vignette

## exprso 0.4.3
---------------------
* New `getWeights` method
    * New method gets LASSO weights for `ExprsModel` and related objects

## exprso 0.4.2
---------------------
* Update `build` modules
    * Extend `buildGLM` to work for any outcome (continuous and factor)
    * New `buildLASSO` module for LASSO and Ridge regression

## exprso 0.4.1
---------------------
* Update `mod` modules
    * Add `alpha` argument to `modRatios` function

## exprso 0.4.0
---------------------
* Package maintenance
    * Clean up message dump within vignettes
    * Update `fsPropd` for compatibility with `propr` version 4.0
    * Update tests for compatibility with `e1071` version 1.7
    * Remove the `gtools` import

## exprso 0.3.9
---------------------
* Update `mod` modules
    * Simplify `modSkew` function to add feature based only on `alpha`
* Update `fs` modules
    * Add `fsPathClassRFE` to the `?fs` help file

## exprso 0.3.8
---------------------
* Update `mod` modules
    * New `modSkew` function skews sample vectors by a factor range

## exprso 0.3.7
---------------------
* Update `exprso` import
    * Fix bug caused by importing an annotation table

## exprso 0.3.6
---------------------
* General changes
    * `exprso` can now swap "Control" and "Case" labels with `switch` argument
    * `exprso` no longer gives error when importing tabular annotations
    * `exprso` can now predict any column from tabular annotations

## exprso 0.3.5
---------------------
* Update `split` modules
    * Add `splitBalanced` to balance both the training and validation sets

## exprso 0.3.4
---------------------
* Update `mod` modules
    * Add `modInclude` function to select specific features from a data set
    * Rename `modShuffle` to `modSample`

## exprso 0.3.3
---------------------
* General changes
    * `exprso` no longer alters labels already called "Control" and "Case"
* Update `split` modules
    * `splitStratify` now runs quietly
    * `splitSample` NULL validation set `cat` changed to a warning

## exprso 0.3.2
---------------------
* Update `fs` modules
    * Add `keep` argument to subset features after performing selection
    * Use of `keep` can speed up model deployment for lage data sets

## exprso 0.3.1
---------------------
* Update `mod` modules
    * Add `base = exp(1)` argument as default for `modTransform` method
* Update `fs` modules
    * Add `modRatios` argument to `fsPropd` for data provided as ratios
* Update `split` modules
    * Add `splitBy` function to split data by user-defined group
* General changes
    * Improve performance by making fewer subsets with `fs.` and `build.`
    * Improve performance by making fewer subsets with `modHistory`

## exprso 0.3.0
---------------------
* Update `mod` modules
    * Replace `modScale` with a new procedure that uses `runif`
    * Extend `modScale` with another method that uses `rnorm`
    * Add `size` argument to `modShuffle`
* Update `fs` modules
    * New `fsRankProd` method ranks features using `RankProd` package
    * Extend `fsStats` to use `wilcox.test` and `f.test`
    * Add `suppressMessages` to `fsPropd` module
* Update `build` modules
    * Sink printed output from `nnet` function
    * Tweak default `buildANN` arguments
* Major `verbose = FALSE` update
    * Add `verbose` to toggle `defaultArg` and `forceArg` messages
    * Force `verbose = FALSE` for `calcStats` calls from `plCV`
    * Make `verbose = FALSE` the default for `plGrid`
* General changes
    * Add `progress` function to list of global functions
    * Add `suppressMessages` to `plMonteCarlo` and `plNested` calls
    * Embed `progress` bar in `plMonteCarlo` and `plNested` calls
    * Add default arguments to argument handler functions

## exprso 0.2.9
---------------------
* Update `mod` modules
    * New `modShuffle` function shuffles feature order by random sampling
    * New `modScale` function scales sample vectors by a factor range
* Update `build` modules
    * New `buildDT` function builds decision tree model
        * Add `predict.ExprsMachine` and `predict.RegrsModel` support
    * New `buildFRB` function builds fuzzy rule based model
        * Add `predict.ExprsMachine` and `predict.RegrsModel` support
    * New `buildLR` function builds logistic regression model
        * Add `predict.ExprsMachine` support
    * New `buildLM` and `buildGLM` functions
        * Add `predict.RegrsModel` support
* General changes
    * `plMonteCarlo` now accepts a list of `ctrlModSet` handlers
    * `plMonteCarlo` now skips `fs` if `ctrlFS = NULL`
    * `plNested` now skips `fs` if `ctrlFS = NULL`

## exprso 0.2.8
---------------------
* New `modRatios` function recasts data as feature ratios

## exprso 0.2.7
---------------------
* General changes
    * Fix bug in `show.ExprsPredict` method

## exprso 0.2.6
---------------------
* General changes
    * New `show` method for `RegrsArray` objects
    * New `show` method for `RegrsModel` objects
    * New `show` method for `RegrsPredict` objects
    * Review all documentation

## exprso 0.2.5
---------------------
* General changes
    * Validate performance of regression `pl` pipelines
    * Users can now deploy multi-class classifier ensembles
    * Users can now deploy regression ensembles
    * Trim documentation

## exprso 0.2.4
---------------------
* General changes
    * Move `modHistory` into script with other mod modules
    * Move `predict` and `calcStats` into separate files
* Revise split modules
    * Make `splitSample` work with `RegrsArray` objects
    * Note `splitStratify` is classification only
* Revise build modules
    * Fix bug when making `RegrsModel` objects
    * Note `buildLDA` is classification only
    * Note `buildNB` is classification only
* Revise predict
    * New predict method for `RegrsModel` objects
    * New calcStats method for `RegrsPredict` objects
    * Trim documentation

## exprso 0.2.3
---------------------
* Revise build modules
    * All build wrapper functions now support continuous outcomes
* Revise pl modules
    * New `ctrlModSet` function to handle pre-process arguments
    * Let `plMonteCarlo` apply pre-process at each bootstrap

## exprso 0.2.2
---------------------
* General changes
    * Have `exprso` throw error if passed boolean outcomes
    * Add new `RegrsPredict` class with two slots
* Revise pre-filter modules
    * New `modAcomp` makes all sample vectors have same total sum
    * New `modCLR` applies centered log-ratio transformation
* Revise fs modules
    * Pass additional arguments to `fsPropd` function
    * Make `fsStats` run a bit faster
* Revise build modules
    * Add `RegrsModel` support to `build.` constructor
* Revise pl modules
    * Add default argument to `calcMonteCarlo`
    * Add default argument to `calcNested`

## exprso 0.2.1
---------------------
* General changes
    * Fix NOTE associated with `magrittr` import
    * Revise `?fs` and `?build` help files
    * Rename `ExprsCont` class to `RegrsArray`
    * New `RegrsModel` class
* Revise fs modules
    * New `fsCor` function ranks features for `RegrsArray` objects

## exprso 0.2.0
---------------------
* General changes
    * Remove all "Golub_Merge" use in vignettes
    * Revise the exprso overview figure
    * Completely revise all vignettes
    * Completely revise README file
    * New `?pipe` help file
* Revise pre-filter modules
    * `modTMM` now returns effective library size instead of cpm
* Revise fs modules
    * New `fsEdger` function ranks features using exact test
    * `fsEbayes` now uses the correct design matrix
* Revise build modules
    * `build` modules now have own documentation
    * New `?build` help file

## exprso 0.1.9
---------------------
* General changes
    * Reformat NOTES and vignettes for CRAN
    * All Bioconductor packages now listed in Suggests:
    * Remove `golubEsets` dependency
* Improve data import
    * New `exprso` function provides `x` and `y` style import
    * New `ExprsCont` class handles continuous outcomes
* Revise pre-filter modules
    * New `modTMM` function performs TMM normalization
    * Replace formal methods with class checks
    * New `?mod` help file
* Revise split modules
    * `split` now defaults to `percent.include = 67`
    * `split` modules now have own documentation
    * New `?split` help file
* Revise fs modules
    * New `fsPropd` function selects features using `propr` package
    * `fs` modules now have own documentation
    * New `?fs` help file

## exprso 0.1.8
---------------------
* Improvements to `arrayExprs`
    * Fixed an import error when reading directly from file
* Improvements to `calcStats`
    * All `NaN` and `NA` performance metrics now replaced with 0
* Improvements to `fs` functions
    * New `fsInclude` ranks stated features above all others
    * Add `fs.` check for non-numeric features
* Improvements to `pl` functions
    * `pl` functions now store run parameters in results
    * Embedding `plMonteCarlo` now calculates outer-fold accuracy
    * Embedding `plNested` now calculates outer-fold accuracy
* Improvements to `pipeFilter`
    * `top` now chooses simplest classifier in case of tie
* Added F1000Research citation to `inst/CITATION`

## exprso 0.1.7
---------------------
* 2.1-import.R
    * Make `Biobase` and `GEOquery` optional packages
* 5.1-fs.R
    * Make `limma` an optional package
    * Renamed from 5.1-fs-binary.R
    * Add `fsNULL` method to return feature input unaltered
    * Add `fsANOVA` feature selection method
* 5.2-build.R
    * Renamed from 5.2-build-binary.R
* 7.2-plGrid.R
    * Add back-end argument handler `makeGridFromArgs`
    * Add default arguments for `how = buildSVM`
* 7.3-plGridMulti.R
    * New function for 1-vs.all classification with `ctrlFS`
* 7.4-plMonteCarlo.R
    * Renamed from 7.3-plMonteCarlo.R
    * `ctrlGS` will now accept any `pl` function
    * Use `rbind.fill` instead of `rbind` to prevent error
* 7.5-plNested.R
    * Renamed from 7.4-plNested.R
    * `ctrlGS` will now accept any `pl` function
    * Use `rbind.fill` instead of `rbind` to prevent error
    * Add `plGridMulti` check to `check.ctrlGS`

## exprso 0.1.6
---------------------
* 1.2-methods.R
    * Rename `getProbeSet` to `getFeatures`
* 2.1-import.R
    * Rename `probes.begin` argument to `begin`
* 4.1-modSwap.R
    * Have `fsPrcomp` use `top` argument
* 4.2-modCluster.R
    * Change `probes` argument to `top`
* 5.1-fs-binary.R
    * Change `probes` argument to `top`
    * Make `fsSample` an `ExprsArray` method
* 5.2-build-binary.R
    * Change `probes` argument to `top`
    * Add `ExprsArray` signature to all `build` methods
    * Add `doMulti` call to `build.` function
* 5.3-doMulti.R
    * Change `probes` argument to `top`
    * Change `what` agument to `method` to allow for `do.call`
    * Use `NA` instead of `NULL` for empty `ExprsMachine`
    * Remove `ExprsMulti` method for `fsSample`
    * Remove `ExprsMulti` method for `fsStats`
* 6-predict.R
    * Use `NA` instead of `NULL` for empty `ExprsMachine`
    * Fix `ExprsMulti` tie breaker warning
* 7.1-plCV.R
    * Change `probes` argument to `top`
* 7.2-plGrid.R
    * Change `probes` argument to `top`
* 7.3-plMonteCarlo.R
    * Change `probes` argument to `top`
* 8.1-pipe.R
    * Change `top.N` argument to `top`
* 8.2-ens.R
    * Change `top.N` argument to `top`
* 9-global.R
    * Add documentation for `data`

## exprso 0.1.5
---------------------
* 1.1-classes.R
    * Add `actual` slot to `ExprsPredict` object
* 1.2-methods.R
    * Add 2D plotting to `plot` method
* 2.1-import.R
    * Force `stringsAsFactors = FALSE`
* 4.1-modSwap.R
    * Clean up plot calls using new `plot` method
* 6-predict.R
    * Add class check for `modHistory` `@reductionModel`
    * Add class check for `predict` `@mach`
    * Pass along known class values to `ExprsPredict` result
    * Remove `calcStats` `array` argument
* 7.1-plCV.R
    * Tidy `calcStats` calls
* 7.2-plGrid.R
    * Tidy `calcStats` calls
* 8.2-ens.R
    * Pass along known class values to `ExprsPredict` result
* 9-global.R
    * Find optimal import combination
* 9-tidy.R
    * Add `pipeSubset` function

## exprso 0.1.4
---------------------
* 2.1-import.R
    * Deprecated `arrayRead` and `arrayEset` functions
    * New `arrayExprs` function adds `data.frame` support
* 3.1-split.R
    * Remove `splitSample` warning
    * Add details to the "Please Read" vignette
* 4.2-cluster.R
    * Add support for numeric vector probes
* 4.3-compare.R
    * Convert `compare` warning into error
* 5.1-fs-binary.R
    * Add imports to NAMESPACE
    * Add `fs.` method to wrap repetitive code
    * Add support for numeric vector probes
    * Remove the `fsStats` "ks-boot" method
    * Remove `fsPenalizedSVM`
* 5.2-build-binary.R
    * Add imports to NAMESPACE
    * Add `build.` method to wrap repetitive code
    * Add support for numeric vector probes
* 7.1-plCV.R
    * Remove `plCV` warning
    * Add details to the "Please Read" vignette
* 7.2-plGrid.R
    * Convert `plGrid` warning into a message
    * Consolidate numeric probes handling
    * Add handling for a list of numeric or character probes
    * Add details to the "Please Read" vignette
* 8.1-pipe.R
    * Convert `pipeFilter` warning into a message
    * Add details to the "Please Read" vignette
* 9-deprecated.R
    * Contains deprecated functions
* 9-tidy.R
    * Add `getArgs`, `defaultArg`, and `forceArg` functions
    * Add `trainingSet`, `validationSet`, and `testSet` functions
    * Add `modSubset` wrapper for `subset` method

## exprso 0.1.3
---------------------
* 1.1-classes.R
    * Fixed warnings and notes
* 1.2-methods.R
    * Fixed warnings and notes
* 2.1-import.R
    * Fixed warnings and notes
* 2.2-misc.R
    * Code renamed to file 2.2-process.R
    * Fixed warnings and notes
* 3.1-split.R
    * Fixed warnings and notes
* 4.1-modSwap.R
    * Store mutated annotation as boolean (not factor)
    * Fixed warnings and notes
* 4.2-modCluster.R
    * Fixed warnings and notes
* 4.3-compare.R
    * Fixed warnings and notes
* 5.1-fs-binary.R
    * Fixed warnings and notes
* 5.2-build-binary.R
    * Fixed warnings and notes
* 5.3-doMulti.R
    * Fixed warnings and notes
* 6-predict.R
    * Fixed warnings and notes
* 8.1-pipe.R
    * Fixed warnings and notes
* 8.2-ens.R
    * Fixed warnings and notes
* 9-global.R
    * Contains global imports

## exprso 0.1.2
---------------------
* 1.2-methods.R
    * Added `subset` method for `ExprsArray` objects
    * Added `subset` method for `ExprsPipeline` objects
* 3-split.R
    * Code renamed to file 3.1-split.R
* 4-conjoin.R
    * Code renamed to file 3.2-conjoin.R
    * `modMutate` merged with `modSwap` as 4.1-modSwap.R
    * `modCluster` method added as 4.2-modCluster.R
    * `compare` method added as 4.3-compare.R
    * `compare` now handles `ExprsMulti` objects
    * `compare` test added to validate method
* 5.1-fs-binary.R
    * `fsStats` has tryCatch to address rare error

## exprso 0.1.1
---------------------
* 1.2-methods.R
    * `summary` method now accommodates lists of vector arguments
* 5.2-build-binary.R
    * Added `buildDNN.ExprsBinary` method
    * Removed `e1071` cross-validation
* 5.3-doMulti.R
    * Added `buildDNN.ExprsMulti` method
* 6-build.R
    * Added `buildDNN` predict clause
* 7.2-plGrid.R
    * `plGrid` method now accommodates lists of vector arguments
    * Removed `e1071` cross-validation

## exprso 0.1.0
---------------------
* Project now organized in a package distribution format
* 0-misc.R
    * Temporarily removed `compare` function
    * Code renamed to file 2.2-misc.R
* 1-classes.R
    * Code divided into files 1.1-classes.R and 1.2-methods.R
    * Remaining 1-classes.R code renamed to 4-conjoin.R
    * Removed `getCases` and `getConts`. Use `[` and `$` instead
    * `getProbeSet` extended to replace `getProbeSummary`
* 2-import.R
    * Code renamed to file 2.1-import.R
* 3-split.R
    * `arraySubset` replaced with `[` and `$` in 1.1-classes.R
    * `splitSample` code heavily edited, including an `all.in` bug fix
    * `splitStratify` now handles `ExprsMulti` objects
* 4-speakEasy.R
    * Temporarily removed `speakEasy` and `abridge` functions
* 5-fs.R
    * Code renamed to file 5.1-fs-binary.R
* 6-build.R
    * `reRank` function serializes `doMulti` fs added to 5.3-doMulti.R
    * `fsSample` and `fsStats` now have `ExprsMulti` methods
    * Some code move to 5.2-build-binary.R and 5.3-doMulti.R
    * Remaining 6-build.R code renamed to 6-predict.R
* 7-pl.R
    * Code divided into a separate file for each `pl` method
    * Replaced ctrlGS (ctrlGridSearch) with ctrlPL (ctrlPipeLine)
    * `plCV`, fixed "1-subject artifact"" with `drop = FALSE`
    * `plNested`, fixed "1-subject artifact" with `drop = FALSE`
    * `plNested` argument checks moved to separate function
* 8-ens.R
    * Code divided into files 8.1-pipe.R and 8.2-ens.R
    * Removed `pipeSubset`. Use `[` and `$` instead
