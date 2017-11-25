#' Filter \code{ExprsPipeline} Object
#'
#' \code{pipeFilter} subsets an \code{ExprsPipeline} object.
#'
#' The filter process occurs in three steps. However, the user may skip
#'  any one of these steps by setting the respective argument to \code{0}.
#'  First, a threshold filter gets imposed. Any model with a performance
#'  less than the threshold filter, \code{how}, gets excluded. Second,
#'  a ceiling filter gets imposed. Any model with a performance greater
#'  than the ceiling filter, \code{gate}, gets excluded. Third, an
#'  arbitrary subset occurs. The top N models in the \code{ExprsPipeline}
#'  object get selected based on the argument \code{top}. However,
#'  in the case that the \code{@@summary} slot contains the column "boot",
#'  \code{pipeFilter} selects the top N models per bootstrap.
#'
#' \code{pipeFilter} will apply this filter based on the performance
#'  metrics listed in the \code{colBy} argument. Listing multiple columns
#'  will result in a filter based on the product of all listed columns.
#'  To weigh one metric over another, list that column more times.
#'
#' @param object An \code{\link{ExprsPipeline-class}} object.
#' @param colBy A character vector or string. Specifies column(s) to use when
#'  filtering by model performance. Listing multiple columns will result
#'  in a filter based on the product all listed columns.
#' @param how,gate A numeric scalar. Arguments between 0 and 1 will impose
#'  a threshold or ceiling filter, respectively, based on the raw value of
#'  \code{colBy}. Arguments between 1 and 100 will impose a filter based on
#'  the percentile of \code{colBy}. The user may also provide "midrange",
#'  "median", or "mean" as an argument for these filters.
#' @param top A numeric scalar. Determines the top N models based on
#'  \code{colBy} to include after the threshold and ceiling filters.
#'  In the case that the \code{@@summary} slot contains the column "boot",
#'  this selects the top N models for each unique bootstrap.
#' @return An \code{\link{ExprsPipeline-class}} object.
#' @export
pipeFilter <- function(object, colBy = "valid.acc", how = 0, gate = 0, top = 0){

  classCheck(object, "ExprsPipeline",
             "This function is applied to the results of ?pl.")

  if(identical(colBy, 0)) return(object) # for ExprsEnsemble
  if(!all(colBy %in% colnames(object@summary))){
    stop("Some 'colBy' measures not in data.")
  }

  # Calculate emergent top accuracy measure as product of colBy columns
  accMeasures <- apply(object@summary[colBy], MARGIN = 1, prod)

  # Impose threshold filter if how != 0
  if(how != 0){

    if(0 < how & how <= 1){ threshold <- how
    }else if(1 < how & how <= 100){ threshold <- quantile(accMeasures, how / 100)
    }else if(how == "midrange"){ threshold <- (max(accMeasures) + min(accMeasures)) / 2
    }else if(how == "median"){ threshold <- median(accMeasures)
    }else if(how == "mean"){ threshold <- mean(accMeasures)
    }else{ stop("Uh oh! Selected 'how' not recognized!")}

    # Filter top accuracy measure as product of colBy columns
    object@summary <- object@summary[accMeasures >= threshold, ]
    object@machs <- object@machs[accMeasures >= threshold]
    accMeasures <- accMeasures[accMeasures >= threshold]
  }

  # Impose ceiling filter if gate != 0
  if(gate != 0){

    if(0 < gate & gate <= 1){ ceiling <- gate
    }else if(1 < gate & gate <= 100){ ceiling <- quantile(accMeasures, gate / 100)
    }else if(gate == "midrange"){ ceiling <- (max(accMeasures) + min(accMeasures)) / 2
    }else if(gate == "median"){ ceiling <- median(accMeasures)
    }else if(gate == "mean"){ ceiling <- mean(accMeasures)
    }else{ stop("Uh oh! Selected 'gate' not recognized!")}

    # Filter top accuracy measure as product of colBy columns
    object@summary <- object@summary[accMeasures <= ceiling, ]
    object@machs <- object@machs[accMeasures <= ceiling]
    accMeasures <- accMeasures[accMeasures <= ceiling]
  }

  # Select top based on presence of a boot column
  if("boot" %in% colnames(object@summary)){

    # For each B boot, select 'top' @machs
    index <- unlist(

      lapply(unique(object@summary$boot),
             function(boot){

               # Calculate total number of gridpoints for this boot
               if(top > sum(object@summary$boot == boot)){

                 message("Provided 'top' too large for boot ", boot,
                         ". Using all gridpoints instead.")
                 top <- 0
               }

               # Order 'top' accMeasures for this boot
               if(top == 0) top <- sum(object@summary$boot == boot)
               topMachs <- order(accMeasures[object@summary$boot == boot],
                                 decreasing = TRUE)[1:top]
               rownames(object@summary[object@summary$boot == boot,])[topMachs]
             }
      )
    )

  }else{

    # Calculate total number of gridpoints for ExprsPipeline object
    if(top > nrow(object@summary)){

      message("Provided 'top' too large for this ExprsPipeline object.",
              "Using all gridpoints instead.")
      top <- 0
    }

    # Order 'top' accMeasures for entire object
    if(top == 0) top <- nrow(object@summary)
    topMachs <- order(accMeasures,
                      decreasing = TRUE)[1:top]
    index <- rownames(object@summary)[topMachs]
  }

  # Filter ExprsPipeline object
  final <- rownames(object@summary) %in% index
  object@summary <- object@summary[final, ]
  object@machs <- object@machs[final]

  return(object)
}

#' Rename "boot" Column
#'
#' \code{pipeUnboot} renames the "boot" column summary to "unboot".
#'
#' This method provides a convenient adjunct to \code{\link{pipeFilter}} owing to
#'  how \code{exprso} handles \code{ExprsPipeline} objects with a "boot" column.
#'
#' @param object An \code{\link{ExprsPipeline-class}} object.
#' @return An \code{\link{ExprsPipeline-class}} object.
#' @export
pipeUnboot <- function(object){

  classCheck(object, "ExprsPipeline",
             "This function is applied to the results of ?pl.")

  if("boot" %in% colnames(object@summary)){

    # Rename 'boot' column to 'unboot'
    colnames(object@summary)[colnames(object@summary) == "boot"] <- "unboot"
  }

  return(object)
}

#' Build Ensemble
#'
#' \code{buildEnsemble} builds an ensemble from \code{ExprsModel} or
#'  \code{ExprsPipeline} objects. See Details.
#'
#' This function can combine any number of model objects into an ensemble.
#'  These models do not necessarily have to derive from the same \code{build}
#'  method. In this way, it works like \code{\link{conjoin}}.
#'
#' This function can also build an ensemble from pipeline objects. It does
#'  this by calling \code{\link{pipeFilter}}, then joining the remaining models
#'  into an ensemble. As an adjunct to this method, consider first combining
#'  multiple pipeline objects with \code{\link{conjoin}}.
#'
#' @inheritParams pipeFilter
#' @param ... Additional \code{ExprsModel} objects to use in the ensemble.
#'  Argument applies to the \code{\link{ExprsModel-class}} method only.
#' @return An \code{\link{ExprsEnsemble-class}} object.
#' @export
setGeneric("buildEnsemble",
           function(object, ...) standardGeneric("buildEnsemble")
)

#' @describeIn buildEnsemble Method to build ensemble from \code{ExprsModel} objects.
#' @export
setMethod("buildEnsemble", "ExprsModel",
          function(object, ...){ # args to include additional ExprsMachine objects

            conjoin(object, ...)
          }
)

#' @describeIn buildEnsemble Method to build ensemble from \code{ExprsPipeline} objects.
#' @export
setMethod("buildEnsemble", "ExprsPipeline",
          function(object, colBy = 0, how = 0, gate = 0, top = 0){

            object <- pipeFilter(object, colBy = colBy, how = how, gate = gate, top = top)
            new("ExprsEnsemble", machs = unlist(object@machs))
          }
)
