## exprso 0.0.0.9000
---------------------
* Project now organized in a package distribution format.
* 1-classes.R
..* Code divided into files 1.1-classes.R and 1.2-conjoin.R.
..* Removed `getCases` and `getConts`. Use `[` and `$` instead.
..* `getProbeSet` extended to replace `getProbeSummary`.
* 3-split.R
..* `arraySubset` replaced with `[` and `$` in 1.1-classes.R.
..* `splitStratify` now properly designated an ExprsBinary method.
