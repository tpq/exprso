<!-- README.md is generated from README.Rmd. Please edit that file -->
Quick start
-----------

Welcome to the `exprso` GitHub page! Let's get started.

``` r
library(devtools)
devtools::install_github("tpq/exprso")
library(exprso)
```

Importing data
--------------

To import data, we use the `exprso` function. This function has two arguments.

``` r
data(iris)
array <- exprso(iris[1:80, 1:4], iris[1:80, 5])
```

    ## [1] "Preparing data for binary classification."

Pre-processing data
-------------------

Functions with a `mod` prefix pre-process the data.

``` r
array <- modTransform(array)
array <- modNormalize(array, c(1, 2))
```

Split data
----------

Functions with a `split` prefix split the data into training and test sets.

``` r
arrays <- splitSample(array, percent.include = 67)
array.train <- arrays$array.train
array.test <- arrays$array.valid
```

Select features
---------------

Functions with a `fs` prefix select features.

``` r
array.train <- fsStats(array.train, top = 0, how = "t.test")
```

Build models
------------

Functions with a `build` prefix build models.

``` r
mach <- buildSVM(array.train,
                 top = 50,
                 kernel = "linear",
                 cost = 1)
```

    ## Setting probability to TRUE (forced behavior, cannot override)...
    ## Setting cross to 0 (forced behavior, cannot override)...

``` r
pred <- predict(mach, array.train)
```

    ## Individual classifier performance:
    ## Arguments not provided in an ROCR AUC format. Calculating accuracy outside of ROCR...
    ## Classification confusion table:
    ##          actual
    ## predicted Control Case
    ##   Control      29    0
    ##   Case          0   25
    ##   acc sens spec
    ## 1   1    1    1

``` r
pred <- predict(mach, array.test)
```

    ## Individual classifier performance:
    ## Arguments not provided in an ROCR AUC format. Calculating accuracy outside of ROCR...
    ## Classification confusion table:
    ##          actual
    ## predicted Control Case
    ##   Control      21    0
    ##   Case          0    5
    ##   acc sens spec
    ## 1   1    1    1

``` r
calcStats(pred)
```

Deploy pipelines
----------------

Functions with a `pl` prefix deploy high-throughput learning pipelines.

``` r
pl <- plGrid(array.train,
             array.test,
             how = "buildSVM",
             top = c(2, 4),
             kernel = "linear",
             cost = 10^(-3:3),
             fold = NULL)
```

``` r
pl
```

    ## Accuracy summary (complete summary stored in @summary slot):
    ## 
    ##      build top kernel  cost train.acc train.sens train.spec train.auc
    ## 1 buildSVM   2 linear 0.001  0.537037          0          1         0
    ## 2 buildSVM   4 linear 0.001  0.537037          0          1         0
    ## 3 buildSVM   2 linear 0.010  1.000000          1          1         1
    ## 4 buildSVM   4 linear 0.010  1.000000          1          1         1
    ##   valid.acc valid.sens valid.spec valid.auc
    ## 1 0.8076923          0          1         0
    ## 2 0.8076923          0          1         0
    ## 3 1.0000000          1          1         1
    ## 4 1.0000000          1          1         1
    ## ...
    ##       build top kernel cost train.acc train.sens train.spec train.auc
    ## 11 buildSVM   2 linear  100         1          1          1         1
    ## 12 buildSVM   4 linear  100         1          1          1         1
    ## 13 buildSVM   2 linear 1000         1          1          1         1
    ## 14 buildSVM   4 linear 1000         1          1          1         1
    ##    valid.acc valid.sens valid.spec valid.auc
    ## 11         1          1          1         1
    ## 12         1          1          1         1
    ## 13         1          1          1         1
    ## 14         1          1          1         1
    ## 
    ## Machine summary (all machines stored in @machs slot):
    ## 
    ## ##Number of classes: 2 
    ## @preFilter summary: 4 2 
    ## @reductionModel summary: logical logical 
    ## @mach class: svm.formula svm 
    ## ...
    ## ##Number of classes: 2 
    ## @preFilter summary: 4 4 
    ## @reductionModel summary: logical logical 
    ## @mach class: svm.formula svm

Read the exprso vignettes for more details.
