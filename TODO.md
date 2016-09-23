## Works in progress
---------------------
* `regrso` expansion:
  * [] Add `RegrsArray` class
  * [] Add `RegrsModel` class
  * [] Watch out for `x[,,drop=TRUE]` errors
  * [] Make ExprsBinary $defineCase a factor?
    * This change alone crashes downstream code
  * [x] Make `calcStats` tidy
  * [x] Add 2D plot
  * [x] Rename `getProbeSet` function.
  * [x] Rename `probes` as `top`
  * [x] Rename `top.N` as `top`
  
* `ExprsMulti` expansion:
  * [x] Fix `splitStratify` for multi-class
  * [x] Fix `compare` for multi-class
  * [x] Implement all build_ ExprsMulti methods
    * Uses `doMulti` for 1-vs-all classification
    * Add `ExprsModule` predict method
    * Add multi-class `calcStats`
    
* `pl` expansion:
  * [] Consider "random Plains" wrapper
  * [] Consider `plRFE` with embedded RFE
  * [] `plGridMulti`
    * Performs fs before each 1-vs-all build
  * [x] Disable ROC plotting during high-throughput `pl`
  * [] `plGridPlus` with "better" plCV?
    * Would feature select anew at each fold
  
* `fs` expansion:
  * [x] Remove `doMulti` fs methods?
  * [] F-test
  * [] ANOVA
    * Add as a true multi-class method
    
* `build` expansion:
  * [] Logistic regression
  * [] Democratic SVM
  * [] KNN
