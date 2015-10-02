###########################################################
### Define Exprso package specific object classes

setClass("ExprsArray",
         slots = c(exprs = "matrix", # Stores most current expression matrix
                   annot = "data.frame", # Stores most current annotation data.frame
                   preFilter = "ANY", # Catalogs a history of changes
                   reductionModel = "ANY", # Catalogs a history of changes
                   n.comp = "ANY" # Catalogs a history of changes
         )
)

setClass("ExprsMachine",
         slots = c(preFilter = "ANY",
                   reductionModel = "ANY",
                   n.comp = "ANY",
                   mach = "ANY"
         )
)

setClass("ExprsPipeline",
         slots = c(summary = "ANY",
                   machs = "ANY"
         )
)

setClass("ExprsEnsemble",
         slots = c(machs = "ANY")
)

setClass("ExprsPredict",
         slots = c(pred = "factor",
                   decision.values = "ANY",
                   probability = "ANY"
         )
)

setMethod("show", "ExprsArray",
          function(object){
            
            cat("@exprs summary:",
                nrow(object@exprs), "probes by", ncol(object@exprs), "subjects\n")
            
            cat("@annot summary:",
                nrow(object@annot), "subjects by", ncol(object@annot), "annotations\n")
            
            cat("@preFilter summary:",
                unlist(lapply(object@preFilter, length)), "\n")
            
            cat("@reductionModel summary:",
                unlist(lapply(object@reductionModel, class)), "\n")
            
            cat("@n.comp summary:", unlist(object@n.comp), "\n")
          }
)

setMethod("show", "ExprsMachine",
          function(object){
            
            cat("@preFilter summary:",
                unlist(lapply(object@preFilter, length)), "\n")
            
            cat("@reductionModel summary:",
                unlist(lapply(object@reductionModel, class)), "\n")
            
            cat("@n.comp summary:", unlist(object@n.comp), "\n")
            
            cat("@mach class:", class(object@mach), "\n")
          }
)

setMethod("show", "ExprsPipeline",
          function(object){
            
            cat("Accuracy summary (complete summary stored in @summary slot):\n\n")
            if(nrow(object@summary) > 8){
              
              print(head(object@summary, 4))
              cat("...\n")
              print(tail(object@summary, 4))
              cat("\n")
              
            }else{
              
              print(object@summary)
              cat("\n")
            }
            
            cat("Machine summary (all machines stored in @machs slot):\n\n")
            if(length(object@machs) > 1){
              
              show(object@machs[[1]])
              cat("...\n")
              show(object@machs[[length(object@machs)]])
              cat("\n")
              
            }else{
              
              show(object@machs[[1]])
              cat("\n")
            }
          }
)

setMethod("show", "ExprsEnsemble",
          function(object){
            
            cat("Machine summary (all machines stored in @machs slot):\n\n")
            if(length(object@machs) > 1){
              
              show(object@machs[[1]])
              cat("...\n")
              show(object@machs[[length(object@machs)]])
              cat("\n")
              
            }else{
              
              show(object@machs[[1]])
              cat("\n")
            }
          }
)

setMethod("show", "ExprsPredict",
          function(object){
            
            cat("@pred summary:", as.numeric(object@pred == "Case"), "\n")
            cat("@decision.values summary:", colnames(object@decision.values), "\n")
            cat("@probability summary:", colnames(object@probability), "\n")
          }
)

###########################################################
### Define functions for plotting and summarizing data

setMethod("plot", signature(x = "ExprsArray", y = "missing"),
          function(x, i = 1, j = 2, k = 3, colors, shapes){
            
            require(lattice)
            
            # Extract components i, j, and k
            df <- data.frame("i" = t(x@exprs)[, i], "j" = t(x@exprs)[, j], "k" = t(x@exprs)[, k])
            colnames(df) <- colnames(t(x@exprs)[, c(i, j, k)])
            colnames(df) <- paste0(c("i_", "j_", "k_"), colnames(df))
            
            if(missing(colors)){
              
              cases <- x@annot$defineCase %in% "Case"
              colors <- ifelse(cases, "red", "black")
            }
            
            if(missing(shapes)){
              
              shapes <- 19
            }
            
            # Plot k ~ i + j in 3D
            func <- as.formula(paste(colnames(df)[3], "~", colnames(df)[1], "+", colnames(df)[2], collapse = ""))
            print(cloud(func, data = df, col = colors, pch = shapes))
          }
)

setMethod("summary", "ExprsArray",
          function(object){
            
            # For large datasets, plot only a sub-sample
            if(length(object@exprs) > 1000){
              
              v <- sample(as.vector(object@exprs), 1000)
              
            }else{
              
              v <- as.vector(object@exprs)
            }
            
            # Prepare layout for multiple plots
            layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
            
            # Plot charts
            qqnorm(v)
            plot(density(v), main = "Density Plot")
            boxplot(v, horizontal = TRUE, main = "Box Plot", xlab = "Expression Values")
            
            # Print per-subject summary
            summary(object@exprs)
          }
)

setMethod("summary", "ExprsPipeline",
          function(object){
            
            # Index which columns contain performance metrics
            index <- grepl("train.", colnames(object@summary)) | grepl("valid.", colnames(object@summary))
            
            # Calculate means and sds for performance metrics
            performance <- list("means" = apply(object@summary[, index], MARGIN = 2, mean),
                                "sds" = apply(object@summary[, index], MARGIN = 2, sd))
            
            # Tabulate parameter frequency
            parameters <- lapply(object@summary[, !index], table)
            
            # Append performance metric summary with parameter frequencies
            summary <- append(performance, parameters)
            
            return(summary)
          }
)

###########################################################
### Define functions for retrieving and combining data

getCases <- function(array){
  
  array <- new("ExprsArray",
               exprs = as.matrix(array@exprs[, array@annot$defineCase %in% "Case"]),
               annot = array@annot[array@annot$defineCase %in% "Case", ],
               preFilter = array@preFilter,
               reductionModel = array@reductionModel,
               n.comp = array@n.comp)
  
  return(array)
}

getConts <- function(array){
  
  array <- new("ExprsArray",
               exprs = as.matrix(array@exprs[, array@annot$defineCase %in% "Control"]),
               annot = array@annot[array@annot$defineCase %in% "Control", ],
               preFilter = array@preFilter,
               reductionModel = array@reductionModel,
               n.comp = array@n.comp)
  
  return(array)
}

setGeneric("getProbeSet",
           function(object, ...){
             standardGeneric("getProbeSet")
           }
)

setMethod("getProbeSet", "ExprsArray",
          function(object){
            
            return(rownames(object@exprs))
          }
)

setMethod("getProbeSet", "ExprsMachine",
          function(object){
            
            return(object@preFilter[[length(object@preFilter)]])
          }
)

setMethod("getProbeSet", "ExprsPipeline",
          function(object, index){
            
            return(object@machs[[index]]@preFilter[[length(object@machs[[index]]@preFilter)]])
          }
)

setMethod("getProbeSet", "ExprsEnsemble",
          function(object, index){
            
            return(object@machs[[index]]@preFilter[[length(object@machs[[index]]@preFilter)]])
          }
)

setGeneric("getProbeSummary",
           function(object, ...){
             standardGeneric("getProbeSummary")
           }
)

setMethod("getProbeSummary", "ExprsPipeline",
          function(object){
            
            probes <- unlist(lapply(object@machs, getProbeSet))
            probes <- table(probes)[order(table(probes), decreasing = TRUE)]
            
            return(probes)
          }
)

setMethod("getProbeSummary", "ExprsEnsemble",
          function(object){
            
            probes <- unlist(lapply(object@machs, getProbeSet))
            probes <- table(probes)[order(table(probes), decreasing = TRUE)]
            
            return(probes)
          }
)

setGeneric("conjoin",
           function(object, ...){
             standardGeneric("conjoin")
           }
)

setMethod("conjoin", "ExprsPipeline",
          function(object, ...){
            
            args <- list(...)
            index <- unlist(lapply(args, function(arg) class(arg) == "ExprsPipeline"))
            
            args.summary <- append(list(object@summary), lapply(args[index], function(pl) pl@summary))
            args.machs <- append(list(object@machs), lapply(args[index], function(pl) pl@machs))
            
            pl <- new("ExprsPipeline",
                      summary = do.call(rbind, args.summary),
                      machs = unlist(args.machs)
            )
            
            return(pl)
          }
)

setMethod("conjoin", "ExprsPipeline",
          function(object, ...){
            
            args <- list(...)
            index <- unlist(lapply(args, function(arg) class(arg) == "ExprsEnsemble"))
            
            args.machs <- append(list(object@machs), lapply(args[index], function(pl) pl@machs))
            
            new("ExprsEnsemble",
                machs = unlist(args.machs)
            )
          }
)

###########################################################
### Read, subset, and split ExprsArray objects

arrayEset <- function(eSet, col.defineBy, case.include, cont.include){
  
  require(Biobase)
  
  # Build an ExprsArray object from the eSet object
  array <- new("ExprsArray",
               exprs = exprs(eSet),
               annot = eSet@phenoData@data,
               preFilter = NULL,
               reductionModel = NULL,
               n.comp = NULL)
  
  # Filter out samples in @annot not matching case.include or cont.include
  array@annot <- array@annot[array@annot[, col.defineBy] %in% c(case.include, cont.include), ]
  
  # Collapse multiple case and cont labels into a binary
  array@annot[, "defineCase"] <- ifelse(array@annot[, col.defineBy] %in% case.include, "Case", "Control")
  
  # Subset @exprs to filter samples no longer found in @annot
  array@exprs <- array@exprs[, rownames(array@annot)]
  
  # Change hyphens for underscores as necessary for some downstream functions
  rownames(array@exprs) <- gsub("-", "_", rownames(array@exprs))
  
  return(array)
}

arrayRead <- function(file, probes.begin, col.subjectID, col.defineBy, case.include, cont.include, ...){ # args to read.delim
  
  # Pass any additional arguments to read.delim
  table <- read.delim(file, stringsAsFactors = FALSE, ...)
  
  # Label table by subject ID
  rownames(table) <- table[, col.subjectID]
  
  # Separate the expression data from annotation data
  array <- new("ExprsArray",
               exprs = t(table[, probes.begin:ncol(table)]),
               annot = table[, 1:(probes.begin-1)],
               preFilter = NULL,
               reductionModel = NULL,
               n.comp = NULL)
  
  # Filter out samples in @annot not matching case.include or cont.include
  array@annot <- array@annot[array@annot[, col.defineBy] %in% c(case.include, cont.include), ]
  
  # Collapse multiple case and cont labels into a binary
  array@annot[, "defineCase"] <- ifelse(array@annot[, col.defineBy] %in% case.include, "Case", "Control")
  
  # Subset @exprs to filter samples no longer found in @annot
  array@exprs <- array@exprs[, rownames(array@annot)]
  
  # Change hyphens for underscores as necessary for some downstream functions
  rownames(array@exprs) <- gsub("-", "_", rownames(array@exprs))
  
  return(array)
}

arraySubset <- function(array, col.subsetBy, set.include){
  
  # Filter out samples not matching set.include
  array@annot <- array@annot[array@annot[, col.subsetBy] %in% set.include, ]
  
  # Subset @exprs to filter samples no longer provided by @annot
  array@exprs <- array@exprs[, rownames(array@annot)]
  
  return(array)
}

# NOTE: to stratify by CASE/CONTROL status only, set col.stratBy = NULL
# NOTE: validated to work with c("Sex", "AGE") on 2015/02/08
splitStrat <- function(array, percent.include = 67, col.stratBy = NULL, bin, breaks, ...){ # args to cut
  
  require(sampling)
  
  if(percent.include < 1 | percent.include > 100) stop("You must choose an inclusion percentage between 1-100!")
  
  if(!is.null(col.stratBy)){
    
    # Perform pre-stratification binning
    vals <- lapply(seq_along(bin),
                   function(i){
                     
                     # For each col.stratBy argument, bin when appropriate
                     if(bin[i]) return(cut(array@annot[, col.stratBy[i]], breaks = breaks[[i]], ...))
                     if(!bin[i]) return(array@annot[, col.stratBy[i]])
                   })
    
    # Add "defineCase" values to the vals bin list
    vals <- c(list(array@annot[, "defineCase"]), vals)
    
    # Name vals according to source
    names(vals) <- c("defineCase", col.stratBy)
    
    # Build a sampling data.frame
    df <- data.frame(vals, row.names = rownames(array@annot), stringsAsFactors = FALSE)
    
    # The order function does not take a list as an argument, so we will use do.call to build an index from vals
    index <- do.call(order, vals)
    df <- df[index, ]
    
    # Remove any row with NAs
    index.NAs <- apply(df, 1, function(row) !NA %in% row)
    df <- df[index.NAs, ]
    
    cat("\nPre-stratification table:\n\n")
    print(table(df))
    
    # Manipulate order of apply(df) so that size vector matches strata expectations
    sizes <- apply(table(df[, c("defineCase", rev(col.stratBy))]), MARGIN = -1, FUN = min)
    sizes <- round(sizes * percent.include/100)
    sizes <- rep(sizes, 2)
    
    # Provide error for anticipated stratum of size 0
    if(0 %in% sizes) stop("This function cannot create a stratum of size 0. You must subset before stratification!")
    
    # Perform stratification with remaining non-NA col.stratBy values (e.g. those introduced by binning)
    s <- sampling::strata(df, stratanames = colnames(df), size = sizes, method = "srswor")
    if(!identical(rownames(s), rownames(df)[s$ID_unit])) stop("Uh-oh, something went wrong! Error Code: 500")
    
    cat("\nWeighted stratification results:\n\n")
    print(table(s[, colnames(df)]))
  }
  
  if(is.null(col.stratBy)){
    
    # Build a sampling data.frame
    df <- data.frame("defineCase" = array@annot$defineCase, row.names = rownames(array@annot), stringsAsFactors = FALSE)
    
    cat("\nPre-stratification table:\n\n")
    print(table(df))
    
    # Compute strata sizes
    sizes <- min(table(df))
    sizes <- round(sizes * percent.include/100)
    sizes <- rep(sizes, 2)
    
    # Provide error for anticipated stratum of size 0
    if(0 %in% sizes) stop("This function cannot create a stratum of size 0. You must subset before stratification!")
    
    # Stratify
    s <- sampling::strata(df, stratanames = colnames(df), size = sizes, method = "srswor")
    rownames(s) <- rownames(df)[s$ID_unit]
    
    cat("\nWeighted stratification results:\n\n")
    print(table(s[, colnames(df)]))
  }
  
  # Build array.train from the stratified sample
  cat("\nBuilding the stratified training set...\n\n")
  array.train <- new("ExprsArray",
                     exprs = array@exprs[, rownames(s)],
                     annot = array@annot[rownames(s), ],
                     preFilter = array@preFilter,
                     reductionModel = array@reductionModel,
                     n.comp = array@n.comp)
  
  # Dump remaining samples into array.valid
  cat("\nBuilding a validation set with remaining samples...\n\n")
  array.valid <- new("ExprsArray",
                     exprs = array@exprs[, !colnames(array@exprs) %in% colnames(array.train@exprs)],
                     annot = array@annot[!rownames(array@annot) %in% rownames(array.train@annot), ],
                     preFilter = array@preFilter,
                     reductionModel = array@reductionModel,
                     n.comp = array@n.comp)
  
  return(list("array.train" = array.train, "array.valid" = array.valid))
}

# Set bag = FALSE to sample training set WITHOUT replacement
# Set bag = TRUE to sample training set WITH replacement (will likely include less than percent.include)
#   NOTE: The demi-holdout will contain all subjects not randomly sampled
# Set percent.include = 100 to return a NULL demi-holdout array
splitBag <- function(array, percent.include = 67, bag = TRUE){
  
  warning("This method is not truly random; at least one case and one control will always appear in demi-holdout!")
  
  if(percent.include < 1 | percent.include > 100){
    
    stop("You must choose an inclusion percentage between 1-100!")
  }
  
  if(!"Case" %in% array@annot$defineCase & "Control" %in% array@annot$defineCase){
    
    stop("Provided ExprsArray object must contain both 'Case' and 'Control' subjects!")
  }
  
  # Calculate size of bootstrap set
  size <- round((ncol(array@exprs) * percent.include)/100, digits = 0)
  
  at.least.1 <- FALSE # Set at.least.1 to FALSE by default
  counter <- 1
  while(!at.least.1){ # If there is not at least one CASE and at least one CONTROL...
    
    # Terminate after 10 iterations
    counter <- counter + 1
    if(counter > 10) stop("splitRandom could not find a solution. Check the supplied parameters.")
    
    # Attempt a simple random sample with or without replacement depending on 'bag' argument
    boot <- sample(colnames(array@exprs), size = size, replace = bag)
    
    # Check whether to terminate the while loop
    if("Case" %in% array@annot[boot, "defineCase"] & "Control" %in% array@annot[boot, "defineCase"]){
      
      at.least.1 <- TRUE
      
    }else{
      
      at.least.1 <- FALSE
    }
  }
  
  # Build bootstrap set
  cat("\nBuilding the random training set...\n\n")
  array.boot <- new("ExprsArray",
                    exprs = array@exprs[, boot], # as.matrix not needed because always at least 2 columns
                    annot = array@annot[boot, ],
                    preFilter = array@preFilter,
                    reductionModel = array@reductionModel,
                    n.comp = array@n.comp)
  
  # Build demi-holdout set based on the size of the bootstrap set
  if(ncol(array.boot@exprs) < ncol(array@exprs)){
    
    cat("\nBuilding a validation set with remaining samples...\n\n")
    array.demi <- new("ExprsArray",
                      exprs = as.matrix(array@exprs[, !colnames(array@exprs) %in% boot]),
                      annot = array@annot[!rownames(array@annot) %in% boot, ],
                      preFilter = array@preFilter,
                      reductionModel = array@reductionModel,
                      n.comp = array@n.comp)
    
  }else{
    
    cat("\nBuilding a NULL validation set...\n\n")
    array.demi <- NULL
  }
  
  return(list("array.train" = array.boot, "array.valid" = array.demi))
}

###########################################################
### Functions for pre-processing ExprsArray objects

# Sets minimum and maximum thresholds, then removes genes with low MAX:MIN relationship
modFilter <- function(array, threshold, maximum, beta1, beta2, displayAll = FALSE){
  
  # Prepare layout for multiple plots
  layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))
  
  # Plot before-filter charts
  if(!displayAll){
    
    if(length(array@exprs) > 1000){
      
      v <- sample(as.vector(array@exprs), 1000)
      
    }else{
      
      v <- as.vector(array@exprs)
    }
    
  }else{
    
    v <- as.vector(array@exprs)
  }
  
  qqnorm(v, main = "Normal Q-Q Plot (Before)")
  plot(density(v), main = "Density Plot (Before)")
  boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
  
  # Set minimum and maximum threshold
  array@exprs[array@exprs < threshold] <- threshold
  array@exprs[array@exprs > maximum] <- maximum
  
  # Calculate RANGE and MAX:MIN
  index.beta1 <- apply(array@exprs, MARGIN = 1, function(gene) max(gene) - min(gene))
  index.beta2 <- apply(array@exprs, MARGIN = 1, function(gene) max(gene) / min(gene))
  
  # INCLUDE those with variance GREATER THAN beta1 AND beta2
  array@exprs <- array@exprs[index.beta1 > beta1 & index.beta2 > beta2, ]
  
  # Plot after-filter charts
  if(!displayAll){
    
    if(length(array@exprs) > 1000){
      
      v <- sample(as.vector(array@exprs), 1000)
      
    }else{
      
      v <- as.vector(array@exprs)
    }
    
  }else{
    
    v <- as.vector(array@exprs)
  }
  
  qqnorm(v, main = "Normal Q-Q Plot (After)")
  plot(density(v), main = "Density Plot (After)")
  boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")
  
  return(array)
}

modTransform <- function(array, displayAll = FALSE){
  
  # Prepare layout for multiple plots
  layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))
  
  # Plot before-log charts
  if(!displayAll){
    
    if(length(array@exprs) > 1000){
      
      v <- sample(as.vector(array@exprs), 1000)
      
    }else{
      
      v <- as.vector(array@exprs)
    }
    
  }else{
    
    v <- as.vector(array@exprs)
  }
  
  qqnorm(v, main = "Normal Q-Q Plot (Before)")
  plot(density(v), main = "Density Plot (Before)")
  boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
  
  # Perform log transformation
  array@exprs <- log(array@exprs, base = 2)
  
  # Plot after-log charts
  if(!displayAll){
    
    if(length(array@exprs) > 1000){
      
      v <- sample(as.vector(array@exprs), 1000)
      
    }else{
      
      v <- as.vector(array@exprs)
    }
    
  }else{
    
    v <- as.vector(array@exprs)
  }
  
  qqnorm(v, main = "Normal Q-Q Plot (After)")
  plot(density(v), main = "Density Plot (After)")
  boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")
  
  return(array)
}

# MARGIN = 1 implies normalize by gene vector, MARGIN = 2 implies normalize by subject vector
# Provide MARGIN = c(1, 2) for both
modNormalize <- function(array, MARGIN = c(1, 2), displayAll = FALSE){
  
  # Prepare layout for multiple plots
  layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))
  
  # Plot before-normalize charts
  if(!displayAll){
    
    if(length(array@exprs) > 1000){
      
      v <- sample(as.vector(array@exprs), 1000)
      
    }else{
      
      v <- as.vector(array@exprs)
    }
    
  }else{
    
    v <- as.vector(array@exprs)
  }
  
  qqnorm(v, main = "Normal Q-Q Plot (Before)")
  plot(density(v), main = "Density Plot (Before)")
  boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
  
  # Perform normalization
  if(1 %in% MARGIN) array@exprs <- t(apply(array@exprs, MARGIN = 1, function(gene) (gene - mean(gene)) / sd(gene)))
  if(2 %in% MARGIN) array@exprs <- apply(array@exprs, MARGIN = 2, function(samp) (samp - mean(samp)) / sd(samp))
  
  # Plot after-normalize charts
  if(!displayAll){
    
    if(length(array@exprs) > 1000){
      
      v <- sample(as.vector(array@exprs), 1000)
      
    }else{
      
      v <- as.vector(array@exprs)
    }
    
  }else{
    
    v <- as.vector(array@exprs)
  }
  
  qqnorm(v, main = "Normal Q-Q Plot (After)")
  plot(density(v), main = "Density Plot (After)")
  boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")
  
  return(array)
}

# Experimental function to mutate case subjects into noisy positives, false positives, or fixed out-groups
modMutate <- function(array, percent = 10, how = "fp", theta = 1){
  
  if(percent < 1 | percent > 100) stop("You must choose an inclusion percentage between 1-100!")
  
  # Mutate a percent of case subjects
  cases <- array@annot$defineCase %in% "Case"
  mut.size <- round(ncol(array@exprs[, cases]) * percent/100)
  mut.name <- sample(colnames(array@exprs[, cases]), size = mut.size, replace = FALSE)
  temp1 <- fsPrcomp(array, probes = 0, n.comp = 3)
  
  if(how == "rp.1"){
    
    for(mut.col in mut.name){
      
      # Randomize the per-subject expression values to create "subject noise"
      array@exprs[, mut.col] <- sample(array@exprs[, mut.col])
    }
    
  }else if(how == "rp.2"){
    
    # Weigh probability of noise contribution by CASE/CONTROL ratio
    probs <- ifelse(cases, 1/(2*ncol(array@exprs[, cases])), 1/(2*ncol(array@exprs[, !cases])))
    
    for(mut.col in mut.name){
      
      # For each probe, draw a random expression value as expressed by any subject
      array@exprs[, mut.col] <- apply(array@exprs, MARGIN = 1, sample, size = 1, prob = probs)
    }
    
  }else if(how == "fp"){
    
    for(mut.col in mut.name){
      
      # For each probe, draw a random expression value as expressed by a control subject
      array@exprs[, mut.col] <- apply(array@exprs[, !cases], MARGIN = 1, sample, size = 1)
    }
    
  }else if(how == "ng"){
    
    # Calculate per-probe case means and sds
    means <- apply(array@exprs[, cases], MARGIN = 1, mean)
    sds <- apply(array@exprs[, cases], MARGIN = 1, sd)
    
    # Calculate fixed "third group" per-probe means and sds as falling on a distrubtion of case means and sds
    # By this method, the mean of means remains the same while still imposing "third group" differences
    ng.means <- unlist(lapply(1:length(means), function(i) rnorm(1, mean = means[i], sd = sds[i]/theta)))
    
    # For every mutant
    for(mut.col in mut.name){
      
      # For each probe, draw a random value based on the per-probe "third group" means and case sds
      array@exprs[, mut.col] <- unlist(lapply(1:length(ng.means), function(i) rnorm(1, mean = ng.means[i], sd = sds[i])))
    }
    
  }else if(how == "tg"){
    
    # Weigh third group mean contribution by CASE/CONTROL ratio
    probs <- ifelse(cases, 1/(2*ncol(array@exprs[, cases])), 1/(2*ncol(array@exprs[, !cases])))
    
    # Calculate per-probe means and sds using ALL subjects
    means <- apply(array@exprs, MARGIN = 1, weighted.mean, w = probs)
    sds <- apply(array@exprs, MARGIN = 1, sd)
    
    # Calculate fixed "third group" per-probe means and sds as falling on a distrubtion of case means and sds
    # By this method, the mean of means remains the same while still imposing "third group" differences
    tg.means <- unlist(lapply(1:length(means), function(i) rnorm(1, mean = means[i], sd = sds[i]/theta)))
    
    # For every mutant
    for(mut.col in mut.name){
      
      # For each probe, draw a random value based on the per-probe "third group" means and case sds
      array@exprs[, mut.col] <- unlist(lapply(1:length(tg.means), function(i) rnorm(1, mean = tg.means[i], sd = sds[i])))
    }
    
  }else{
    
    stop("Provided how not recognized. Select from 'rp.1', 'rp.2', 'fp', 'ng', or 'tg'.\n")
  }
  
  # Store Boolean index of mutated subjects in @annot
  array@annot$mutated <- rownames(array@annot) %in% mut.name
  
  # Visualize results
  temp2 <- fsPrcomp(array, probes = 0, n.comp = 3)
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  plot(temp1@exprs[1, ], temp1@exprs[2, ], col = ifelse(cases, "red", "black"), pch = ifelse(array@annot$mutated, 22, 16), main = "Before", xlab = "PCA1", ylab = "PCA2")
  plot(temp2@exprs[1, ], temp2@exprs[2, ], col = ifelse(cases, "red", "black"), pch = ifelse(array@annot$mutated, 22, 16), main = "After", xlab = "PCA1", ylab = "PCA2")
  plot(temp1@exprs[1, ], temp1@exprs[3, ], col = ifelse(cases, "red", "black"), pch = ifelse(array@annot$mutated, 22, 16), main = "Before", xlab = "PCA1", ylab = "PCA3")
  plot(temp2@exprs[1, ], temp2@exprs[3, ], col = ifelse(cases, "red", "black"), pch = ifelse(array@annot$mutated, 22, 16), main = "After", xlab = "PCA1", ylab = "PCA3")
  
  return(array)
}

###########################################################
### Feature selection methods for ExprsArray objects

# NOTE: @preFilter stores probe selection history and exists to recreate fs during svmPredict
# NOTE: IF probes = 0, include ALL probes when building data; otherwise, select top N occurring first
# NOTE: The function fsStats() will auto-sort array@exprs by p-value before supplying output
# NOTE: IF probes is a character vector, include only these probes when building data
# NOTE: 'probes' argument refers to what you feed INTO the function, not what you expect OUT

fsPrcomp <- function(array, probes, n.comp, ...){ # args to prcomp
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
    data <- t(array@exprs[probes, ])
  }
  
  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){
    
    data <- t(array@exprs[probes, ])
  }
  
  # ATTENTION: We want dependent variables as columns
  reductionModel <- prcomp(data, ...)
  
  cat("\nDimension reduction model summary:\n\n")
  print(summary(reductionModel))
  
  # ATTENTION: The value of predict(reductionModel, data) equals $x
  # @preFilter stores probes used to build reductionModel (i.e. as passed on by 'probes' argument)
  # @n.comp stores the number of top components to include after building reductionModel
  # This information will automatically distill the data when calling svmPredict
  array <- new("ExprsArray",
               exprs = t(reductionModel$x[, 1:n.comp]),
               annot = array@annot,
               preFilter = append(array@preFilter, list(probes)),
               reductionModel = append(array@reductionModel, list(reductionModel)),
               n.comp = append(array@n.comp, list(n.comp))
  )
  
  return(array)
}

fsStats <- function(array, probes, how, ...){ # args to ks.test, ks.boot, or t-test
  
  require(plyr)
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
  }
  
  # Retrieve case subjectIDs
  cases <- array@annot$defineCase %in% "Case"
  
  # Retrieve control subjectIDs
  conts <- array@annot$defineCase %in% "Control"
  
  if(how == "ks"){
    
    # Perform KS pre-filter and return p-values as a data.frame
    vals <- ldply(probes,
                  function(probe){
                    
                    result <- ks.test(array@exprs[probe, cases], array@exprs[probe, conts], ...)
                    data.frame(probe, "p.value" = result$p.value)
                  }
    )
  }else if(how == "ks.boot"){
    
    require(Matching)
    
    # Perform KS pre-filter and return p-values as a data.frame
    vals <- ldply(probes,
                  function(probe){
                    
                    result <- ks.boot(array@exprs[probe, cases], array@exprs[probe, conts], ...)
                    data.frame(probe, "p.value" = result$ks.boot.pvalue)
                  }
    )
  }else if(how == "t-test"){
    
    # Perform t-test pre-filter and return p-values as a data.frame
    vals <- ldply(probes,
                  function(probe){
                    
                    result <- t.test(array@exprs[probe, cases], array@exprs[probe, conts], ...)
                    data.frame(probe, "p.value" = result$p.value)
                  }
    )
  }else{
    
    stop("Selected 'how' argument not found. See ?fsStats for available methods.")
  }
  
  final <- as.character(vals[order(vals$p.value), "probe"])
  array <- new("ExprsArray",
               exprs = array@exprs[final,],
               annot = array@annot,
               preFilter = append(array@preFilter, list(final)),
               reductionModel = append(array@reductionModel, list(NA)),
               n.comp = append(array@n.comp, list(0))
  )
  
  return(array)
}

fsPenalizedSVM <- function(array, probes, ...){ # args to svm.fs
  
  require(penalizedSVM)
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
    data <- t(array@exprs[probes, ])
  }
  
  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){
    
    data <- t(array@exprs[probes, ])
  }
  
  # Set labels as factor
  labels <- factor(array@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
  levels(labels) <- c(-1, 1)
  
  # Run svm.fs()
  fs <- svm.fs(x = data, y = labels, ...)
  
  final <- names(fs$model$w)
  array <- new("ExprsArray",
               exprs = array@exprs[final, ],
               annot = array@annot,
               preFilter = append(array@preFilter, list(final)),
               reductionModel = append(array@reductionModel, list(NA)),
               n.comp = append(array@n.comp, list(0))
  )
  
  return(array)
}

fsPathClassRFE <- function(array, probes, ...){ # args to fit.rfe
  
  require(pathClass)
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
    data <- t(array@exprs[probes, ])
  }
  
  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){
    
    data <- t(array@exprs[probes, ])
  }
  
  # Set labels as factor
  labels <- factor(array@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
  
  # NOTE: RFE as assembled by pathClass is via a linear kernel only
  # NOTE: By default, fit.rfe iterates through C = 10^c(-3:3)
  # Run fit.rfe()
  rfe <- fit.rfe(x = data, y = labels, ...)
  
  final <- rfe$features
  array <- new("ExprsArray",
               exprs = array@exprs[final, ],
               annot = array@annot,
               preFilter = append(array@preFilter, list(final)),
               reductionModel = append(array@reductionModel, list(NA)),
               n.comp = append(array@n.comp, list(0))
  )
  
  return(array)
}

fsEbayes <- function(array, probes, ...){ # args to ebayes
  
  require(limma)
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
  }
  
  # Set up and perform eBayes
  design <- as.matrix(ifelse(array@annot$defineCase == "Case", 1, 0))
  colnames(design) <- "CaseVCont"
  fit <- lmFit(array@exprs[probes, ], design)
  ebaye <- ebayes(fit, ...)
  
  # Sort probes
  vals <- data.frame("probe" = rownames(ebaye$p.value), "p.value" = ebaye$p.value[, 1])
  final <- as.character(vals[order(vals$p.value), "probe"])
  
  array <- new("ExprsArray",
               exprs = array@exprs[final,],
               annot = array@annot,
               preFilter = append(array@preFilter, list(final)),
               reductionModel = append(array@reductionModel, list(NA)),
               n.comp = append(array@n.comp, list(0))
  )
  
  return(array)
}

# NOTE: mRMR.classic is prone to crashing when supplied a very large 'feature_count' argument
fsMrmre <- function(array, probes, feature_count, ...){ # args to mRMR.classic
  
  require(mRMRe)
  
  if(missing(feature_count)){
    
    feature_count <- 64
  }
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
    data <- t(array@exprs[probes, ])
  }
  
  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){
    
    data <- t(array@exprs[probes, ])
  }
  
  # Set up and perform mRMR
  labels <- as.numeric(array@annot$defineCase == "Case")
  mRMRdata <- mRMR.data(data = data.frame(labels, data))
  mRMRout <- mRMR.classic(data = mRMRdata, target_indices = 1, feature_count = feature_count, ...)
  
  # Sort probes
  final <- as.vector(apply(solutions(mRMRout)[[1]], 2, function(x, y) { return(y[x]) }, y = mRMRe::featureNames(mRMRdata)))
  
  # Remove the X appended by the mRMRe::featureNames function to all features named with a number
  final[!final %in% rownames(array@exprs)] <- unlist(lapply(final[!final %in% rownames(array@exprs)], function(id) substr(id, start = 2, stop = nchar(id))))
  
  array <- new("ExprsArray",
               exprs = array@exprs[final,],
               annot = array@annot,
               preFilter = append(array@preFilter, list(final)),
               reductionModel = append(array@reductionModel, list(NA)),
               n.comp = append(array@n.comp, list(0))
  )
  
  return(array)
}

# NOTE: col.modelBy is a character string of columns to use in the model (i.e. between '+' signs)
# NOTE: the column defineCase is AUTOMATICALLY included in model
# NOTE: to denote random-effects, include "|" or "||" in string element
# e.g. if you want lmer(probe ~ Days + (Days | Subject), sleepstudy)
#   let col.modelBy = c("Days", "(Days | Subject)", "sleepstudy")
fsLmerTest <- function(array, probes, col.modelBy = colnames(array@annot), ...){ # args to lmer
  
  require(lmerTest)
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
    data <- t(array@exprs[probes, ])
  }
  
  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){
    
    data <- t(array@exprs[probes, ])
  }
  
  # Join annotations with prepared data
  if(identical(rownames(array@annot), rownames(data))){ df <- cbind(array@annot, data)
  }else{ stop("Uh oh! We have a mismatch between annotations and prepared data!")}
  
  # Prepare col.modelBy for formula
  x <- paste(unique(c("defineCase", col.modelBy)), collapse = " + ")
  
  # Build a model for each element of probes
  models <- lapply(probes, function(y) lmerTest::lmer(formula(paste(y, "~", x)), df, ...))
  
  # Retrieve coefficients
  results <- lapply(models, summary)
  finals <- lapply(results, coef)
  
  # Retrieve p-values
  px <- unlist(lapply(finals, function(final) final[grepl("defineCase", rownames(final)), "Pr(>|t|)"]))
  
  # Sort probes
  vals <- data.frame("probe" = probes, "p.value" = px)
  final <- as.character(vals[order(vals$p.value), "probe"])
  
  array <- new("ExprsArray",
               exprs = array@exprs[final,],
               annot = array@annot,
               preFilter = append(array@preFilter, list(final)),
               reductionModel = append(array@reductionModel, list(NA)),
               n.comp = append(array@n.comp, list(0))
  )
  
  return(array)
}

###########################################################
### Build and deploy individual classifiers

# NOTE: User has one more opportunity to subset the probe set before building machine
# NOTE: @preFilter, @reductionModel, and @n.comp data will get passed to ExprsMachine object
# NOTE: Once predict is called, validation set will get processed accordingly
buildSVM <- function(array, probes, ...){ # args to svm
  
  require(e1071)
  
  args <- as.list(substitute(list(...)))[-1]
  
  if(!"cross" %in% names(args)){
    
    cat("Setting 'cross' to 10 (default behavior, override explicitly)...\n")
    args <- append(args, list("cross" = 10))
  }  
  
  if(!"probability" %in% names(args)){
    
    cat("Setting 'probability' to TRUE (default behavior, override explicitly)...\n")
    args <- append(args, list("probability" = TRUE))
  }
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
    data <- t(array@exprs[probes, ])
  }
  
  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){
    
    data <- t(array@exprs[probes, ])
  }
  
  # Set labels as factor
  labels <- factor(array@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
  
  # Perform SVM via ~ method (permits plotting)
  df <- data.frame(data, "defineCase" = labels)
  args <- append(list("formula" = defineCase ~ ., "data" = df), args)
  model <- do.call(svm, args)
    
  # Carry through and append fs history as stored in the ExprsArray object
  # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
  machine <- new("ExprsMachine",
                 preFilter = append(array@preFilter, list(probes)),
                 reductionModel = append(array@reductionModel, list(NA)),
                 n.comp = append(array@n.comp, list(0)),
                 mach = model
  )
  
  return(machine)
}

# NOTE: Consider varying decay = c(.5, 1) and size = c(1:7)
buildANN <- function(array, probes, ...){
  
  require(nnet)
  
  args <- as.list(substitute(list(...)))[-1]
  
  if(!"size" %in% names(args)){
    
    cat("Setting 'size' to 1 (default behavior, override explicitly)...\n")
    args <- append(args, list("size" = 1))
  }
  
  if(!"range" %in% names(args)){
    
    cat("Setting 'range' to 1/max(|x|) (default behavior, override explicitly)...\n")
    args <- append(args, list("range" = 1/max(abs(as.vector(array@exprs)))))
  }
  
  if(!"decay" %in% names(args)){
    
    cat("Setting 'decay' to 0.5 (default behavior, override explicitly)...\n")
    args <- append(args, list("decay" = 0.5))
  }
  
  if(!"maxit" %in% names(args)){
    
    cat("Setting 'maxit' to 1000 (default behavior, override explicitly)...\n")
    args <- append(args, list("maxit" = 1000))
  }
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
    data <- t(array@exprs[probes, ])
  }
  
  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){
    
    data <- t(array@exprs[probes, ])
  }
  
  # Set labels as factor
  labels <- factor(array@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
  
  # Perform ANN via ~ method
  df <- data.frame(data, "defineCase" = labels)
  args <- append(list("formula" = defineCase ~ ., "data" = df), args)
  model <- do.call(nnet, args)
    
  # Carry through and append fs history as stored in the ExprsArray object
  # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
  machine <- new("ExprsMachine",
                 preFilter = append(array@preFilter, list(probes)),
                 reductionModel = append(array@reductionModel, list(NA)),
                 n.comp = append(array@n.comp, list(0)),
                 mach = model
  )
  
  return(machine)
}

buildRF <- function(array, probes, ...){
  
  require(randomForest)
  
  args <- as.list(substitute(list(...)))[-1]
  
  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){
    
    if(probes == 0) probes <- nrow(array@exprs)
    probes <- rownames(array@exprs[1:probes, ])
    data <- t(array@exprs[probes, ])
  }
  
  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){
    
    data <- t(array@exprs[probes, ])
  }
  
  # Set labels as factor
  labels <- factor(array@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
  
  # Perform RF via ~ method
  df <- data.frame(data, "defineCase" = labels)
  args <- append(list("formula" = defineCase ~ ., "data" = df), args)
  model <- do.call(randomForest, args)
  
  # Carry through and append fs history as stored in the ExprsArray object
  # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
  machine <- new("ExprsMachine",
                 preFilter = append(array@preFilter, list(probes)),
                 reductionModel = append(array@reductionModel, list(NA)),
                 n.comp = append(array@n.comp, list(0)),
                 mach = model
  )
  
  return(machine)
}

modHistory <- function(object, reference){
  
  if(class(object) != "ExprsArray"){
    
    stop("Uh oh! You can only modify the history of an ExprsArray object.")
  }
  
  if(!is.null(object@preFilter)){
    
    # If reference@preFilter has less (or equal) history than object@preFilter
    if(length(reference@preFilter) <= length(object@preFilter)){
      
      stop("The provided object does not have less history than the reference object.")
    }
    
    # If history of reference@preFilter is not also in the object@preFilter
    if(!identical(object@preFilter, reference@preFilter[1:length(object@preFilter)])){
      
      stop("The provided object history does not match the reference history.")
    }
  }
  
  # Index first non-overlapping history
  index <- length(object@reductionModel) + 1
  
  # Manipulate object according to history stored in reference
  for(i in index:length(reference@reductionModel)){
    
    # If the i-th fs did NOT involve dimension reduction
    if(any(is.na(reference@reductionModel[[i]]))){
      
      # Build new object
      object <- new("ExprsArray",
                    exprs = object@exprs[reference@preFilter[[i]], ], # Update @exprs
                    annot = object@annot, # Preserve @annot
                    preFilter = append(object@preFilter, list(reference@preFilter[[i]])), # Append history
                    reductionModel = append(object@reductionModel, list(reference@reductionModel[[i]])),
                    n.comp = append(object@n.comp, list(reference@n.comp[[i]])))
      
    }else{
      
      # Build data according to the i-th probe set
      data <- t(object@exprs[reference@preFilter[[i]], ])
      
      # Then, apply the i-th reduction model
      comps <- predict(reference@reductionModel[[i]], data)
      
      # Build new object
      object <- new("ExprsArray",
                    exprs = t(comps[, 1:reference@n.comp[[i]]]), # Update @exprs
                    annot = object@annot, # Preserve @annot
                    preFilter = append(object@preFilter, list(reference@preFilter[[i]])), # Append history
                    reductionModel = append(object@reductionModel, list(reference@reductionModel[[i]])),
                    n.comp = append(object@n.comp, list(reference@n.comp[[i]])))
    }
  }
  
  return(object)
}

# Build predict function for all ExprsMachine objects regardless of @mach class
# NOTE: The validation set should not get modified once separated from the training set
# NOTE: This function performs only rudimentary checks for this type of error
# NOTE: predict.nnet class predictions match @decision.values, but differ when using ROCR
# NOTE: predict.svm class predictions do not match @decision.values, but match ROCR
setMethod("predict", "ExprsMachine",
          function(object, array, ...){ # args to predict(ExprsMachine, array, ...)
            
            if(class(array) != "ExprsArray"){
              
              stop("Uh oh! You can only use an ExprsMachine to predict on an ExprsArray object.")
            }
            
            # Reproduce the history stored in the ExprsMachine object
            array <- modHistory(array, reference = object)
            
            # Build data from modHistory output
            data <- t(array@exprs)
            
            # Set labels as factor
            labels <- factor(array@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
            
            if("svm" %in% class(object@mach)){
              
              require(e1071)
              
              args <- as.list(substitute(list(...)))[-1]
              
              if(!"decision.values" %in% names(args)){
                
                cat("Setting 'decision.values' to TRUE (default behavior, override explicitly)...\n")
                args <- append(args, list("decision.values" = TRUE))
              }
              
              if(!"probability" %in% names(args)){
                
                cat("Setting 'probability' to TRUE (default behavior, override explicitly)...\n")
                args <- append(args, list("probability" = TRUE))
              }
              
              # Predict SVM via ~ method (permits plotting)
              df <- data.frame(data, "defineCase" = labels)
              args <- append(list("object" = object@mach, "newdata" = df), args)
              pred <- do.call(predict, args)
              
              # Extract 'probabilities' attribute
              px <- attr(pred, "probabilities")
              px <- px[, c("Control", "Case")]
              
              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              if(!is.null(px)){
                
                dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
                colnames(dv) <- "Case/Control"
                
              }else{
                
                dv <- NULL
              }
              
              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))
              names(pred) <- rownames(df)
            }
            
            if("nnet" %in% class(object@mach)){
              
              require(nnet)
              
              # Predict ANN via ~ method
              df <- data.frame(data, "defineCase" = labels)
              pred <- predict(object@mach, df, type = "class")
              
              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              px <- predict(object@mach, df, type = "raw")
              if(px[which.min(px)] > .5){ labMin <- ifelse(pred[which.min(px)] == "Case", "Control", "Case")
              }else{ labMin <- pred[which.min(px)] }
              if(px[which.max(px)] < .5){ labMax <- ifelse(pred[which.max(px)] == "Case", "Control", "Case")
              }else{ labMax <- pred[which.max(px)] }
              px <- cbind(1 - px, px)
              colnames(px) <- c(labMin, labMax)
              px <- px[, c("Control", "Case")] # order columns
              dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
              colnames(dv) <- "Case/Control"
              
              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))
              names(pred) <- rownames(df)
            }
            
            if("randomForest" %in% class(object@mach)){
              
              require(randomForest)
              
              # Predict RF via ~ method
              df <- data.frame(data, "defineCase" = labels)
              pred <- predict(object@mach, df, type = "response")
              
              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              px <- unclass(predict(object@mach, df, type = "prob"))
              px <- px[, c("Control", "Case")] # order columns
              dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
              colnames(dv) <- "Case/Control"
              
              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))
              names(pred) <- rownames(df)
            }
            
            final <- new("ExprsPredict", pred = pred, decision.values = dv, probability = px)
            
            # NOTE: validated to work on 2015/06/30
            # NOTE: rebuilt on 2015/08/12
            cat("Individual classifier performance:\n")
            print(calcStats(final, array))
            
            return(final)
          }
)

calcStats <- function(pred, array, aucSkip = FALSE){
  
  require(ROCR)
  
  layout(matrix(c(1), 1, 1, byrow = TRUE))
  
  # Build 'actual' object as binary
  actual <- ifelse(array@annot$defineCase == "Case", 1, 0)
  
  # If pred contains decision.values and aucSkip = FALSE
  if(!is.null(pred@probability) & !aucSkip){
    
    p <- prediction(pred@probability[, "Case"], as.numeric(actual))
    
    # Plot AUC curve
    perf <- performance(p, measure = "tpr", x.measure = "fpr")
    plot(perf, col = rainbow(10))
    
    # Index optimal cutoff based on Euclidean distance
    index <- which.min(sqrt((1 - perf@y.values[[1]])^2 + (0 - perf@x.values[[1]])^2))
    
    # Calculate performance metrics
    acc <- performance(p, "acc")@y.values[[1]][index]
    sens <- performance(p, "sens")@y.values[[1]][index]
    spec <- performance(p, "spec")@y.values[[1]][index]
    auc <- performance(p, "auc")@y.values[[1]]
    
    return(data.frame(acc, sens, spec, auc))
    
  }else{
    
    # Build 'p' object as binary
    p <- prediction(as.numeric(pred@pred == "Case"), as.numeric(actual))
    
    # Plot AUC curve
    perf <- performance(p, measure = "tpr", x.measure = "fpr")
    plot(perf, col = rainbow(10))
    
    # Index optimal cutoff based on Euclidean distance
    index <- which.min(sqrt((1 - perf@y.values[[1]])^2 + (0 - perf@x.values[[1]])^2))
    
    # Calculate performance metrics
    acc <- performance(p, "acc")@y.values[[1]][index]
    sens <- performance(p, "sens")@y.values[[1]][index]
    spec <- performance(p, "spec")@y.values[[1]][index]
    
    return(data.frame(acc, sens, spec))
  }
}

###########################################################
### Build plGrid

# NOTE: buildSVM will automatically calculate AUC unless 'probability' argument is explicitly FALSE
# NOTE: SVM predicted class membership does seem to differ slightly when 'probability' = TRUE
#   "On the one hand, when probability=FALSE, predict() uses the signs of the decision values.
#   On the other hand, when probability=TRUE, predict() uses a fitted logistic model."
plGrid <- function(array.train, array.valid = NULL, probes, how, ...){ # args to how function
  
  require(plyr)
  
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
  
  # Initialize ExprsMachine container
  models <- NULL
  
  # Build grid
  grid <- expand.grid(append(list("probes" = probes), lapply(args, eval)), stringsAsFactors = FALSE)
  
  # Refine grid for svmBuild
  if(how == "svmBuild"){
    
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
  
  summary <- ddply(grid,
                   .variables = colnames(grid),
                   .fun = function(gridpoint){
                     
                     cat("Now building machine at gridpoint:\n")
                     print(gridpoint)
                     
                     # Format gridpoint args to pass along to build do.call
                     args <- append(list("array" = array.train), as.list(gridpoint))
                     
                     # Build model then add to machs container
                     model <- do.call(what = how, args = args[!is.na(args)])
                     models <<- append(models, list(model))
                     
                     # Predict class labels using the provided training set and calculate accuracy
                     pred.train <- predict(model, array.train)
                     acc <- data.frame("train" = calcStats(pred.train, array.train))
                     
                     # Extract cross-validation accuracy if available (should work even if how = "buildANN"?)
                     if(!is.null(model@mach$tot.accuracy)) acc <- data.frame("cv.train.acc" = model@mach$tot.accuracy / 100, acc)
                     
                     # If a validation set is provided
                     if(!is.null(array.valid)){
                       
                       # Predict class labels using the provided validation set and calculate accuracy
                       pred.valid <- predict(model, array.valid)
                       acc <- data.frame(acc, "valid" = calcStats(pred.valid, array.valid))
                     }
                     
                     df <- data.frame(gridpoint, acc)
                   }
  )
  
  pl <- new("ExprsPipeline",
            summary = summary,
            machs = models
  )
  
  return(pl)
}

###########################################################
### Build plBoot

ctrlSplitSet <- function(func, percent.include, ...){
  
  list("func" = func,
       "percent.include" = percent.include,
       ...)
}

# NOTE: ctrlFS is a way of organizing arguments supplied to fs_ functions
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

plBoot <- function(array, B, ctrlSS, ctrlFS, ctrlGS, save = FALSE){
  
  require(plyr)
  
  # Perform check to make sure the supplied array has enough cases and controls to work meaningfully well
  if(sum(array@annot$defineCase == "Case") < 5 | sum(array@annot$defineCase == "Control") < 5){
    
    stop("Use at least 5 cases and at least 5 controls!\n")
  }
  
  # For each bootstrap
  pls <- lapply(1:B,
                function(boot){
                  
                  ###
                  # (I) For each B boot, index the demi-holdout
                  
                  # Perform some split function (e.g. splitStrat)
                  func <- ctrlSS$func
                  args <- append(list("array" = array), ctrlSS[!ctrlSS %in% func])
                  arrays <- do.call(what = func, args = args)
                  
                  # Build demi-holdout
                  array.boot <- arrays[[1]]
                  
                  # Build bootstrap
                  array.demi <- arrays[[2]]
                  
                  # ERROR CHECK: Does array.demi overlap with array.boot?
                  if(sum(colnames(array.demi@exprs) %in% colnames(array.boot@exprs)) > 0){
                    
                    stop("Uh oh, it seems that the demi-holdout contains at least one non-unique subject! Contact developer to debug.\n")
                  }
                  
                  # Save files
                  if(save){
                    
                    save(array.boot, file = paste0("plBoot ", boot, " (", gsub(":", ".", Sys.time()), ") bootstrap.RData"))
                    save(array.demi, file = paste0("plBoot ", boot, " (", gsub(":", ".", Sys.time()), ") demi-holdout.RData"))
                  }
                  
                  ###
                  # (II) Process array.boot via each ctrlFS in ctrlFS
                  
                  # If ctrlFS is not a list, make it a list
                  if(!"list" %in% lapply(ctrlFS, class)) ctrlFS <- list(ctrlFS)
                  
                  # Perform fs_ function for each argument set in ctrlFS
                  for(i in 1:length(ctrlFS)){
                    
                    func <- ctrlFS[[i]]$func
                    args <- append(list("array" = array.boot), ctrlFS[[i]][!ctrlFS[[i]] %in% func])
                    array.boot <- do.call(what = func, args = args)
                  }
                  
                  ###
                  # (III) Run GridSearch via provided ctrlGS
                  
                  # Perform some gridsearch function (e.g. plSVM)
                  func <- ctrlGS$func
                  args <- append(list("array.train" = array.boot, "array.valid" = array.demi), ctrlGS[!ctrlGS %in% func])
                  pl <- do.call(what = func, args = args)
                  
                  ###
                  # (IV) Append pl@summary with additional information
                  
                  # Calculate the percent of cases in the bootstrap and demi-holdout cohorts
                  tab.boot <- table(array.boot@annot$defineCase)
                  tab.demi <- table(array.demi@annot$defineCase)
                  cohort.stats <- data.frame("boot.cases" = as.vector(tab.boot["Case"] / (tab.boot["Case"] + tab.boot["Control"])),
                                             "demi.cases" = as.vector(tab.demi["Case"] / (tab.demi["Case"] + tab.demi["Control"])))
                  
                  # Append pl@summary
                  pl@summary <- cbind(cohort.stats, pl@summary)
                  pl@summary <- cbind(boot, pl@summary)
                  
                  return(pl)
                }
  )
  
  ###
  # (V) Build ExprsPipeline
  
  pl <- new("ExprsPipeline",
            summary = ldply(pls, function(obj) obj@summary),
            machs = unlist(lapply(pls, function(obj) obj@machs))
  )
  
  return(pl)
}


###########################################################
### Build functions for modifying ExprsPipeline objects

setGeneric("pipeSubset",
           function(object, ...){
             standardGeneric("pipeSubset")
           }
)

setMethod("pipeSubset", "ExprsPipeline",
          function(object, col.subsetBy, set.include){
            
            # Build an index to filter out results not matching set.include
            index <- object@summary[, col.subsetBy] %in% set.include
            
            # Subset @summary using index
            object@summary <- object@summary[index, ]
            
            # Subset @machs using index
            object@machs <- object@machs[index]
            
            return(object)
          }
)

setGeneric("pipeFilter",
           function(object, ...){
             standardGeneric("pipeFilter")
           }
)

# NOTE: IF is.vector(col.accBy), use product of accuracy measures
# NOTE: Increase weight of col.accBy elements with multiple mentions
# NOTE: IF how = 0: do not impose any threshold filter
# NOTE: IF .01 < how < 1: define threshold equal to 'how'
# NOTE: IF 1 < how < 100: define threshold equal to 'how' as a quantile
# NOTE: IF top.N = 0, include ALL @machs when building ensemble
setMethod("pipeFilter", "ExprsPipeline",
          function(object, col.accBy, how = 0, top.N = 0){
            
            # Check if ExprsPipeline contains all col.accBy
            if(!all(col.accBy %in% colnames(object@summary))){
              
              stop("Uh oh! This ExprsPipeline does not contain all provided 'col.accBy' accuracy measures.")
            }
            
            # Calculate emergent top accuracy measure as product of col.accBy columns
            accMeasures <- apply(object@summary[col.accBy], MARGIN = 1, prod)
            
            # Set threshold equal to 'how'
            if(0 <= how & how <= 1){ threshold <- how
            
            # Set threshold equal to 'how' as a quantile
            }else if(1 < how & how <= 100){ threshold <- quantile(accMeasures, how / 100)
            
            }else if(how == "midrange"){ threshold <- (max(accMeasures) + min(accMeasures)) / 2
            
            }else if(how == "median"){ threshold <- median(accMeasures)
            
            }else if(how == "mean"){ threshold <- mean(accMeasures)
            
            }else{ stop("Uh oh! Selected 'how' not recognized!")}
            
            # Filter ExprsPipeline object
            object@summary <- object@summary[accMeasures >= threshold, ]
            object@machs <- object@machs[accMeasures >= threshold]
            
            # Filter top accuracy measure as product of col.accBy columns
            accMeasures <- accMeasures[accMeasures >= threshold]
            
            # If ExprsPipeline comes from plBoot
            if("boot" %in% colnames(object@summary)){
              
              # For each B boot, select 'top.N' @machs
              index <- unlist(
                
                lapply(unique(object@summary$boot),
                       function(boot){
                         
                         # Calculate total number of gridpoints for this boot
                         if(top.N > sum(object@summary$boot == boot)){
                           
                           warning("Provided 'top.N' too large for boot ", boot, ". Using all gridpoints instead.")
                           top.N <- 0
                         }
                         
                         if(top.N == 0) top.N <- sum(object@summary$boot == boot)
                         
                         # Order 'top.N' accMeasures for this boot
                         topMachs <- order(accMeasures[object@summary$boot == boot], decreasing = TRUE)[1:top.N]
                         
                         # Index by rowname for this boot
                         rownames(object@summary[object@summary$boot == boot,])[topMachs]
                       }
                )
              )
              
            }else{
              
              # Calculate total number of gridpoints for ExprsPipeline object
              if(top.N > nrow(object@summary)){
                
                warning("Provided 'top.N' too large for this ExprsPipeline object. Using all gridpoints instead.")
                top.N <- 0
              }
              
              if(top.N == 0) top.N <- nrow(object@summary)
              
              # Order 'top.N' accMeasures for entire object
              topMachs <- order(accMeasures, decreasing = TRUE)[1:top.N]
              
              # Index by rowname
              index <- rownames(object@summary)[topMachs]
            }
            
            # Filter ExprsPipeline object
            final <- rownames(object@summary) %in% index # 'final' depends on initial object@summary definition
            object@summary <- object@summary[final, ]
            object@machs <- object@machs[final]
            
            return(object)
          }
)

###########################################################
### Build and deploy ensemble classifiers

setGeneric("buildEnsemble",
           function(object, ...){
             standardGeneric("buildEnsemble")
           }
)

# Build ExprsEnsemble using an explicit set of ExprsMachine objects
setMethod("buildEnsemble", "ExprsMachine",
          function(object, ...){ # args to include additional ExprsMachine objects
            
            args <- list(...)
            index <- unlist(lapply(args, function(arg) class(arg) == "ExprsMachine"))
            
            new("ExprsEnsemble",
                machs = append(object, args[index])
            )
          }
)

# Build ExprsEnsemble from an ExprsPipeline object
setMethod("buildEnsemble", "ExprsPipeline",
          function(object, col.accBy, how = 0, top.N = 0){
            
            pf <- pipeFilter(object, col.accBy = col.accBy, how = how, top.N = top.N)
            
            new("ExprsEnsemble",
                machs = unlist(pf@machs)
            )
          }
)

# NOTE: As of 08/13/2015, 'simple' voting is the only prediction method implemented
setMethod("predict", "ExprsEnsemble",
          function(object, array, how = "simple", ...){ # args to predict(ExprsMachine, array, ...)
            
            if(class(array) != "ExprsArray"){
              
              stop("Uh oh! You can only use an ExprsEnsemble to predict on an ExprsArray object.")
            }
            
            # Deploy each machine in @machs on the provided ExprsArray
            results <- lapply(object@machs, function(mach) predict(mach, array, ...))
            
            if(how == "simple"){
              
              # Cast a vote (1 for Case, -1 for Control)
              votes <- lapply(results, function(result) ifelse(result@pred == "Case", 1, -1))
              final <- rowSums(data.frame(votes, row.names = names(results[[1]]@pred)))
              
              # Randomly assign ties as "Case" or "Control" based on the proportion of cases
              cases <- sum(array@annot$defineCase == "Case")
              conts <- sum(array@annot$defineCase == "Control")
              tieBreaker <- cases/(cases + conts)
              final[final == 0] <- sample(c(1, -1), length(final[final == 0]), replace = TRUE, prob = c(tieBreaker, 1 - tieBreaker))
              
              # Prepare ExprsPredict object
              pred <- factor(ifelse(final > 0, "Case", "Control"), levels = c("Control", "Case"))
              
              # If there exists non-NULL @decision.values slots
              if(!any(is.null(lapply(results, function(result) result@probability)))){
                
                # Calculate average probabilities
                pxs <- lapply(results, function(result) result@probability)
                Case <- lapply(pxs, function(px) px[, "Case"])
                Case <- rowMeans(as.data.frame(Case))
                Control <- lapply(pxs, function(px) px[, "Control"])
                Control <- rowMeans(as.data.frame(Control))
                px <- cbind(Control, Case)
                
                # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
                dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
                colnames(dv) <- "Case/Control"
                
              }else{
                
                dv <- NULL
                px <- NULL
              }
            }
            
            final <- new("ExprsPredict", pred = pred, decision.values = dv, probability = px)
            
            cat("Ensemble classifier performance:\n")
            print(calcStats(final, array))
            
            return(final)
          }
)
