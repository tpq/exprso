###########################################################
### Cross-validation argument handlers

ctrlSplitSet <- function(func, percent.include, ...){
  
  list("func" = func,
       "percent.include" = percent.include,
       ...)
}

ctrlFeatureSelect <- function(func, probes, ...){
  
  list("func" = func,
       "probes" = probes,
       ...)
}

ctrlGridSearch <- function(func, probes, ...){
  
  list("func" = func,
       "probes" = probes,
       ...)
}

###########################################################
### Inner-loop cross-validation

# Calculates v-fold or leave-one-out cross-validation
# NOTE: set 'fold' = 0 to perform leave-one-out cross-validation
# NOTE: set 'fold' = v to perform v-fold cross-validation
plCV <- function(array, probes, how, fold, ...){ # args to get(how)
  
  if(!is.null(array@preFilter) | !is.null(array@reductionModel)){
    
    warning("plCV can help inform parameter selection but will provide overly optimistic cross-validation accuracies!")
  }
  
  # Extract args from ...
  args <- as.list(substitute(list(...)))[-1]
  
  # Perform LOOCV if 0 fold
  if(fold == 0) fold <- nrow(array@annot)
  
  if(fold > nrow(array@annot)){
    
    warning("Insufficient subjects for v-fold cross-validation. Performing LOOCV instead.")
    fold <- nrow(array@annot)
  }
  
  # Prepare list to receive per-fold subject IDs
  subjects <- vector("list", fold)
  
  # Randomly sample subject IDs
  ids <- sample(rownames(array@annot))
  
  # Initialize while loop
  i <- 1
  
  # Add the ith subject ID to the vth fold
  while(i <= nrow(array@annot)){
    
    subjects[[i %% fold + 1]] <- c(subjects[[i %% fold + 1]], ids[i])
    i <- i + 1
  }
  
  # Prepare vector to receive cv accs
  accs <- vector("numeric", fold)
  
  # Build a machine against the vth fold
  for(v in 1:length(subjects)){
    
    # The leave one out
    array.train <- new(class(array),
                       exprs = as.matrix(array@exprs[, !colnames(array@exprs) %in% subjects[[v]]]),
                       annot = array@annot[!rownames(array@annot) %in% subjects[[v]], ],
                       preFilter = array@preFilter,
                       reductionModel = array@reductionModel)
    
    # Clean up 1-subject artifact
    colnames(array.train@exprs) <- rownames(array.train@annot)
    
    # The left out one
    array.valid <- new(class(array),
                       exprs = as.matrix(array@exprs[, colnames(array@exprs) %in% subjects[[v]]]),
                       annot = array@annot[rownames(array@annot) %in% subjects[[v]], ],
                       preFilter = array@preFilter,
                       reductionModel = array@reductionModel)
    
    # Clean up 1-subject artifact
    colnames(array.valid@exprs) <- rownames(array.valid@annot)
    
    # Prepare args for do.call
    args.v <- append(list("object" = array.train, "probes" = probes), args)
    
    # Build machine
    mach <- do.call(what = how, args = args.v)
    
    # Deploy
    pred <- predict(mach, array.valid, verbose = FALSE)
    
    # Save accuracy
    accs[v] <- calcStats(pred, array.valid, aucSkip = TRUE, plotSkip = TRUE)$acc
    
    cat("Inner-fold", v, "accuracy:", accs[v], "\n")
  }
  
  acc <- mean(accs)
  
  return(acc)
}

###########################################################
### High-throughput classification

# NOTE: buildSVM will automatically calculate AUC unless 'probability' argument is explicitly FALSE
# NOTE: SVM predicted class membership does seem to differ slightly when 'probability' = TRUE
#   "On the one hand, when probability=FALSE, predict() uses the signs of the decision values.
#   On the other hand, when probability=TRUE, predict() uses a fitted logistic model."
# NOTE: set 'fold' = 0 to perform leave-one-out cross-validation
# NOTE: set 'fold' = v to perform v-fold cross-validation
# NOTE: set 'fold' = NULL to skip cross-validation
plGrid <- function(array.train, array.valid = NULL, probes, how, fold = 10, aucSkip = FALSE, verbose = TRUE, ...){ # args to how function
  
  args <- as.list(substitute(list(...)))[-1]
  
  # Prepare the use of a numeric probe vector
  if(is.numeric(probes)){
    
    # Exclude 'probes' values less than or equal total available probes
    probes <- probes[!probes > nrow(array.train@exprs)]
    
    # If no valid 'probes' values remain, use all probes
    if(length(probes) == 0){
      
      warning("Each supplied 'probes' values too large. Using all probes instead.")
      probes <- 0
    }
    
    # Replace any 'probes = 0' arguments with the maximum number of probes
    probes[probes == 0] <- nrow(array.train@exprs)
    
  }else{
    
    # Turn character vector into single list entry
    probes <- as.list(probes)
  }
  
  # Build grid
  grid <- expand.grid(append(list("probes" = probes), lapply(args, eval)), stringsAsFactors = FALSE)
  
  # Refine grid for buildSVM
  if(how == "buildSVM"){
    
    if(!"kernel" %in% names(args)) stop("Uh oh! 'kernel' argument missing!")
    if(!"cost" %in% names(args)) stop("Uh oh! 'cost' argument missing!")
    if("radial" %in% eval(args$kernel)){
      
      if(!"gamma" %in% names(args)) stop("Uh oh! 'gamma' argument missing!")
      grid[grid$kernel %in% "linear", "gamma"] <- NA
    }
    
    if("polynomial" %in% eval(args$kernel)){
      
      if(!"degree" %in% names(args)) stop("Uh oh! 'degree' argument missing!")
      if(!"coef0" %in% names(args)) stop("Uh oh! 'coef0' argument missing!")
      grid[grid$kernel %in% c("linear", "kernel"), "degree"] <- NA
      grid[grid$kernel %in% c("linear", "kernel"), "coef0"] <- NA
    }
    
    grid <- unique(grid)
  }
  
  # Initialize ExprsPipeline summary container
  statistics <- vector("list", nrow(grid))
  
  # Initialize ExprsPipeline machs container
  models <- vector("list", nrow(grid))
  
  # For each gridpoint in grid
  for(i in 1:nrow(grid)){
    
    cat("Now building machine at gridpoint:\n")
    print(grid[i, , drop = FALSE])
    
    # Format gridpoint args to pass along to build do.call
    args <- append(list("object" = array.train), as.list(grid[i, , drop = FALSE]))
    
    # Build and save model
    model <- do.call(what = how, args = args[!is.na(args)])
    models[[i]] <- model
    
    # Predict class labels using the provided training set and calculate accuracy
    pred.train <- predict(model, array.train, verbose = verbose)
    stats <- calcStats(pred.train, array.train, aucSkip = aucSkip, plotSkip = TRUE)
    colnames(stats) <- paste0("train.", colnames(stats))
    acc <- stats
    
    # Extract cross-validation accuracy if available (should work even if how = "buildANN"?)
    if(!is.null(model@mach$tot.accuracy)) acc <- data.frame("cv.train.acc" = model@mach$tot.accuracy / 100, acc)
    
    # If a validation set is provided
    if(!is.null(array.valid)){
      
      # Predict class labels using the provided validation set and calculate accuracy
      pred.valid <- predict(model, array.valid, verbose = verbose)
      stats <- calcStats(pred.valid, array.valid, aucSkip = aucSkip, plotSkip = TRUE)
      colnames(stats) <- paste0("valid.", colnames(stats))
      acc <- data.frame(acc, stats)
    }
    
    # If 'fold' argument is provided
    if(!is.null(fold)){
      
      # Perform leave-one-out or v-fold cross-validation
      args <- append(list("how" = how, "fold" = fold), args)
      names(args)[names(args) == "object"] <- "array"
      cv <- do.call(what = plCV, args = args[!is.na(args)])
      acc <- data.frame("fold" = fold, "train.plCV" = cv, acc)
    }
    
    # Save summary statistics
    statistics[[i]] <- data.frame(grid[i, , drop = FALSE], acc)
  }
  
  pl <- new("ExprsPipeline",
            summary = do.call(rbind, statistics),
            machs = models
  )
  
  return(pl)
}

###########################################################
### Monte Carlo cross-validation

# Performs a sophisticated Monte Carlo style cross-validation
# Function returns a single ExprsPipeline object result
# For single accuracy statistic, call calcMonteCarlo
plMonteCarlo <- function(array, B = 10, ctrlSS, ctrlFS, ctrlGS, save = FALSE){
  
  if(!is.null(array@preFilter) | !is.null(array@reductionModel)){
    
    warning("Prior use of feature selection may result in overly optimistic cross-validation accuracies!")
  }
  
  # For each bootstrap
  pls <- lapply(1:B,
                function(boot){
                  
                  # Perform some split function (e.g. splitStrat)
                  func <- ctrlSS$func
                  args <- append(list("object" = array), ctrlSS[!ctrlSS %in% func])
                  arrays <- do.call(what = func, args = args)
                  array.boot <- arrays[[1]]
                  array.demi <- arrays[[2]]
                  
                  # Save files
                  if(save){
                    
                    save(array.boot, file = paste0("plMonteCarlo ", boot, " (", gsub(":", ".", Sys.time()), ") bootstrap.RData"))
                    save(array.demi, file = paste0("plMonteCarlo ", boot, " (", gsub(":", ".", Sys.time()), ") demi-holdout.RData"))
                  }
                  
                  # Perform fs_ function for each argument set in ctrlFS
                  if(!"list" %in% lapply(ctrlFS, class)) ctrlFS <- list(ctrlFS)
                  for(i in 1:length(ctrlFS)){
                    
                    func <- ctrlFS[[i]]$func
                    args <- append(list("object" = array.boot), ctrlFS[[i]][!ctrlFS[[i]] %in% func])
                    array.boot <- do.call(what = func, args = args)
                  }
                  
                  # Perform some gridsearch function (e.g. plGrid)
                  func <- ctrlGS$func
                  args <- append(list("array.train" = array.boot, "array.valid" = array.demi), ctrlGS[!ctrlGS %in% func])
                  pl <- do.call(what = func, args = args)
                  
                  # Append pl@summary
                  pl@summary <- cbind(boot, pl@summary)
                  
                  return(pl)
                }
  )
  
  pl <- new("ExprsPipeline",
            summary = do.call(rbind, lapply(pls, function(obj) obj@summary)),
            machs = unlist(lapply(pls, function(obj) obj@machs))
  )
  
  return(pl)
}

# Returns a single accuracy statistic based on plMonteCarlo output
calcMonteCarlo <- function(pl, colBy){
  
  if(missing(colBy)) stop("Uh oh! Missing 'colBy' argument.")
  
  if("boot" %in% colnames(pl@summary)){
    
    if("train.plCV" %in% colnames(pl@summary)){
      
      # Prepare container to store validation accuracy
      acc <- vector("numeric", length(unique(pl@summary$boot)))
      
      for(b in 1:length(unique(pl@summary$boot))){
        
        cat("Retrieving best accuracy for boot", b, "...\n")
        
        # Subset only boot 'b'
        boot <- pl@summary[pl@summary$boot == b, ]
        
        # Select best model based on cross-validation accuracy
        best <- boot[which.max(boot$train.plCV), ]
        
        # Save validation accuracy as colBy product
        acc[b] <- apply(best[colBy], MARGIN = 1, prod)
      }
      
    }else{
      
      stop("Uh oh! Supplied data not in expected format. Cannot calculate this cross-validation accuracy.")
    }
    
  }else{
    
    stop("Uh oh! Supplied data not in expected format. Cannot calculate this cross-validation accuracy.")
  }
  
  # Return average validation accuracy
  cat("Averaging best accuracies across all boots...\n")
  return(mean(acc))
}

###########################################################
### Nested cross-validation

# Performs a v-fold or leave-one-out nested cross-validation
# NOTE: set 'fold' = 0 to perform leave-one-out cross-validation
# NOTE: set 'fold' = v to perform v-fold cross-validation
# Function returns a single ExprsPipeline object result
# For single accuracy statistic, call calcNested
# NOTE: will not allow a NULL 'fold' in ctrlGS
# NOTE: will not allow FALSE aucSkip in ctrlGS
plNested <- function(array, fold = 10, ctrlFS, ctrlGS, save = FALSE){
  
  if(!is.null(array@preFilter) | !is.null(array@reductionModel)){
    
    warning("Prior use of feature selection may result in overly optimistic cross-validation accuracies!")
  }
  
  # Perform LOOCV if 0 fold
  if(fold == 0) fold <- nrow(array@annot)
  
  if(fold > nrow(array@annot)){
    
    warning("Insufficient subjects for v-fold cross-validation. Performing LOOCV instead.")
    fold <- nrow(array@annot)
  }
  
  # Prepare list to receive per-fold subject IDs
  subjects <- vector("list", fold)
  
  # Randomly sample subject IDs
  ids <- sample(rownames(array@annot))
  
  # Initialize while loop
  i <- 1
  
  # Add the ith subject ID to the vth fold
  while(i <= nrow(array@annot)){
    
    subjects[[i %% fold + 1]] <- c(subjects[[i %% fold + 1]], ids[i])
    i <- i + 1
  }
  
  # Prepare list to receive ExprsPipeline objects
  pls <- vector("list", fold)
  
  # Perform nested cross-validation
  for(v in 1:length(subjects)){
    
    # The v-th fold
    array.boot <- new(class(array),
                      exprs = as.matrix(array@exprs[, !colnames(array@exprs) %in% subjects[[v]]]),
                      annot = array@annot[!rownames(array@annot) %in% subjects[[v]], ],
                      preFilter = array@preFilter,
                      reductionModel = array@reductionModel)
    
    # Clean up 1-subject artifact
    colnames(array.boot@exprs) <- rownames(array.boot@annot)
    
    # The leave out
    array.demi <- new(class(array),
                      exprs = as.matrix(array@exprs[, colnames(array@exprs) %in% subjects[[v]]]),
                      annot = array@annot[rownames(array@annot) %in% subjects[[v]], ],
                      preFilter = array@preFilter,
                      reductionModel = array@reductionModel)
    
    # Clean up 1-subject artifact
    colnames(array.demi@exprs) <- rownames(array.demi@annot)
    
    # Save files
    if(save){
      
      save(array.boot, file = paste0("plNested ", v, " (", gsub(":", ".", Sys.time()), ") bootstrap.RData"))
      save(array.demi, file = paste0("plNested ", v, " (", gsub(":", ".", Sys.time()), ") demi-holdout.RData"))
    }
    
    # Perform fs_ function for each argument set in ctrlFS
    if(!"list" %in% lapply(ctrlFS, class)) ctrlFS <- list(ctrlFS)
    for(i in 1:length(ctrlFS)){
      
      func <- ctrlFS[[i]]$func
      args <- append(list("object" = array.boot), ctrlFS[[i]][!ctrlFS[[i]] %in% func])
      array.boot <- do.call(what = func, args = args)
    }
    
    # Perform some gridsearch function (e.g. plGrid)
    func <- ctrlGS$func
    args <- append(list("array.train" = array.boot, "array.valid" = array.demi), ctrlGS[!ctrlGS %in% func])
    
    # Mandate v-fold or leave-one-out cross-validation
    if(!"fold" %in% names(args)){
      
      cat("Setting 'fold' to 10 (default behavior, override explicitly)...\n")
      args <- append(args, list("fold" = 10))
    }
    
    # Mandate v-fold or leave-one-out cross-validation
    if(is.null(args$fold)){
      
      cat("Uh oh! This function requires non-NULL 'fold'. Setting 'fold' to 10...\n")
      args$fold <- 10
    }
    
    # Mandate skipping AUC performance estimates
    if(!"aucSkip" %in% names(args)){
      
      cat("Setting 'aucSkip' to TRUE (default behavior, override explicitly)...\n")
      args <- append(args, list("aucSkip" = TRUE))
    }
    
    # Mandate skipping AUC performance estimates
    if(!args$aucSkip){
      
      cat("Uh oh! This function requires TRUE 'aucSkip'. Setting 'aucSkip' to TRUE...\n")
      args$aucSkip <- TRUE
    }
    
    # Save pl object
    pl <- do.call(what = func, args = args)
    pl@summary <- cbind(v, pl@summary)
    pls[[v]] <- pl
  }
  
  pl <- new("ExprsPipeline",
            summary = do.call(rbind, lapply(pls, function(obj) obj@summary)),
            machs = unlist(lapply(pls, function(obj) obj@machs))
  )
}

# Returns a single accuracy statistic based on plNested output
calcNested <- function(pl, colBy){
  
  if(missing(colBy)) stop("Uh oh! Missing 'colBy' argument.")
  
  if("v" %in% colnames(pl@summary)){
    
    if("train.plCV" %in% colnames(pl@summary)){
      
      # Prepare container to store validation accuracy
      acc <- vector("numeric", length(unique(pl@summary$v)))
      
      for(b in 1:length(unique(pl@summary$v))){
        
        cat("Retrieving best accuracy for fold", b, "...\n")
        
        # Subset only fold 'b'
        fold <- pl@summary[pl@summary$v == b, ]
        
        # Select best model based on cross-validation accuracy
        best <- fold[which.max(fold$train.plCV), ]
        
        # Save validation accuracy as colBy product
        acc[b] <- apply(best[colBy], MARGIN = 1, prod)
      }
      
    }else{
      
      stop("Uh oh! Supplied data not in expected format. Cannot calculate this cross-validation accuracy.")
    }
    
  }else{
    
    stop("Uh oh! Supplied data not in expected format. Cannot calculate this cross-validation accuracy.")
  }
  
  # Return average validation accuracy
  cat("Averaging best accuracies across all folds...\n")
  return(mean(acc))
}
