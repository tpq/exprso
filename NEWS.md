## exprso 1.0.0.9000
---------------------
* 5.2-build-binary.R
  * Added `buildDNN.ExprsBinary` method
* 5.3-doMulti.R
  * Added `buildDNN.ExprsMulti` method
* 6-build.R
  * Added `buildDNN` predict clause

## exprso 1.0.0.0
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
  * `reRank` function to serialize `doMulti` fs added to 5.3-doMulti.R.
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
