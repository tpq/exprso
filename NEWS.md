## exprso 0.1.2.9000
---------------------
* 0-global.R
  * Includes all global function definitions.
* 1.1-classes.R
  * Revised documentation.
* 1.2-methods.R
  * Revised documentation.
* 2.1-import.R
  * Revised documentation.
* 2.2-misc.R
  * Code renamed to file 2.2-process.R.
  * Revised documentation.
* 3.1-split.R
  * Revised documentation.
* 4.1-modSwap.R
  * Store mutated annotation as boolean (not factor).
  * Revised documentation.
* 4.2-modCluster.R
  * Revised documentation.
* 4.3-compare.R
  * Revised documentation.
* 5.1-fs-binary.R
  * Revised documentation.
* 5.2-build-binary.R
  * Revised documentation.
* 5.3-doMulti.R
  * Revised documentation.

## exprso 0.1.2
---------------------
* 1.2-methods.R
  * Added `subset` method for `ExprsArray` objects.
  * Added `subset` method for `ExprsPipeline` objects.
* 3-split.R
  * Code renamed to file 3.1-split.R.
* 4-conjoin.R
  * Code renamed to file 3.2-conjoin.R.
  * `modMutate` merged with `modSwap` as 4.1-modSwap.R.
  * `modCluster` method added as 4.2-modCluster.R.
  * `compare` method added as 4.3-compare.R.
  * `compare` now handles `ExprsMulti` objects.
  * `compare` test added to validate method.
* 5.1-fs-binary.R
  * `fsStats` has tryCatch to address rare error.

## exprso 1.0.1
---------------------
* 1.2-methods.R
  * `summary` method now accommodates lists of vector arguments.
* 5.2-build-binary.R
  * Added `buildDNN.ExprsBinary` method.
  * Removed `e1071` cross-validation.
* 5.3-doMulti.R
  * Added `buildDNN.ExprsMulti` method.
* 6-build.R
  * Added `buildDNN` predict clause.
* 7.2-plGrid.R
  * `plGrid` method now accommodates lists of vector arguments.
  * Removed `e1071` cross-validation.

## exprso 0.1.0
---------------------
* Project now organized in a package distribution format.
* 0-misc.R
  * Temporarily removed `compare` function.
  * Code renamed to file 2.2-misc.R.
* 1-classes.R
  * Code divided into files 1.1-classes.R and 1.2-methods.R.
  * Remaining 1-classes.R code renamed to 4-conjoin.R.
  * Removed `getCases` and `getConts`. Use `[` and `$` instead.
  * `getProbeSet` extended to replace `getProbeSummary`.
* 2-import.R
  * Code renamed to file 2.1-import.R.
* 3-split.R
  * `arraySubset` replaced with `[` and `$` in 1.1-classes.R.
  * `splitSample` code heavily edited, including an `all.in` bug fix.
  * `splitStratify` now handles `ExprsMulti` objects.
* 4-speakEasy.R
  * Temporarily removed `speakEasy` and `abridge` functions.
* 5-fs.R
  * Code renamed to file 5.1-fs-binary.R.
* 6-build.R
  * `reRank` function serializes `doMulti` fs added to 5.3-doMulti.R.
  * `fsSample` and `fsStats` now have `ExprsMulti` methods.
  * Some code move to 5.2-build-binary.R and 5.3-doMulti.R.
  * Remaining 6-build.R code renamed to 6-predict.R.
* 7-pl.R
  * Code divided into a separate file for each `pl` method.
  * Replaced ctrlGS (ctrlGridSearch) with ctrlPL (ctrlPipeLine).
  * `plCV`, fixed "1-subject artifact"" with `drop = FALSE`.
  * `plNested`, fixed "1-subject artifact" with `drop = FALSE`.
  * `plNested` argument checks moved to separate function.
* 8-ens.R
  * Code divided into files 8.1-pipe.R and 8.2-ens.R.
  * Removed `pipeSubset`. Use `[` and `$` instead.
