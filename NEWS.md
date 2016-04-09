## exprso 0.0.0.9000
---------------------
* Project now organized in a package distribution format.
* 0-misc.R
..* Temporarily commented out until appropriate documentation gets added.
* 1-classes.R
..* Code divided into files 1.1-classes.R and 1.2-conjoin.R.
..* Removed `getCases` and `getConts`. Use `[` and `$` instead.
..* `getProbeSet` extended to replace `getProbeSummary`.
* 3-split.R
..* `arraySubset` replaced with `[` and `$` in 1.1-classes.R.
..* `splitSample` code heavily edited, including an `all.in` bug fix.
..* `splitStratify` now properly designated an ExprsBinary method.
* 4-speakEasy.R
..* Temporarily commented out until appropriate documentation gets added.
* 5-fs.R
..* Code renamed into file 5.1-fs-binary.R.
* 6-build.R
..* Code divided into files 5.2-build-binary.R and 6- scripts.
# ..* `reRank` added to 6.1-doMulti.R.
# * 7-pl.R
# ..* Code divided into four separate files.
# * 8-ens.R
# ..* `pipeSubset` replaced with `[` and `$` in 1.1-classes.R.
