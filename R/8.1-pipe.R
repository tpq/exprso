###########################################################
### Filter ExprsPipeline object

#' Filter \code{ExprsPipeline} Object
#'
#' \code{pipeFilter} subsets an \code{ExprsPipeline} object.
#'
#' The filter process occurs in three steps. However, the user may skip
#'  any one of these steps by setting the respective argument to \code{0}.
#'
#' First, a threshold filter gets imposed. Any model with a performance
#'  less than the threshold filter, \code{how}, gets excluded. Second,
#'  a ceiling filter gets imposed. Any model with a performance less
#'  than the ceiling filter, \code{gate}, gets excluded. Third, an
#'  arbitrary subset occurs. The top N models in the \code{ExprsPipeline}
#'  object get selected based on the argument \code{top}. However,
#'  in the case that the \code{@@summary} slot contains the column
#'  "boot", \code{pipeFilter} selects the top N models for each unique
#'  bootstrap.
#'
#' \code{pipeFilter} will apply this filter for one or more performance
#'  metrics listed in the \code{colBy} argument. Listing multiple columns
#'  will result in a filter based on a performance metric equal to the
#'  product of all listed performance metrics. To more heavily weigh
#'  one performance metric over another, consider listing that column
#'  more than once.
#'
#' @param object An \code{\link{ExprsPipeline-class}} object.
#' @param colBy A character vector or string. Specifies column(s) to use when
#'  filtering by classifier performance. Listing multiple columns will result
#'  in a filter based on a performance metric equal to the product of those
#'  listed columns.
#' @param how,gate A numeric scalar. Arguments between 0 and 1 will impose
#'  a threshold or ceiling filter, respectively, based on the raw value of
#'  \code{colBy}. Arguments between 1 and 100 will impose a filter based on
#'  the percentile of \code{colBy}. The user may also provide "midrange",
#'  "median", or "mean" as an argument for these filters. Set \code{how = 0}
#'  or \code{gate = 0}, to skip the threshold or ceiling filter,
#'  respectively.
#' @param top A numeric scalar. Determines the top N models based on
#'  \code{colBy} to include after the threshold and ceiling filters.
#'  In the case that the \code{@@summary} slot contains the column
#'  "boot", this determines the top N models for each unique bootstrap.
#'  Set \code{top = 0} to skip this subset.
#'
#' @return An \code{\link{ExprsPipeline-class}} object.
#'
#' @seealso
#' \code{\link{pipeFilter}}\cr
#' \code{\link{pipeUnboot}}\cr
#' \code{\link{plCV}}\cr
#' \code{\link{plGrid}}\cr
#' \code{\link{plGridMulti}}\cr
#' \code{\link{plMonteCarlo}}\cr
#' \code{\link{plNested}}
#'
#' @export
setGeneric("pipeFilter",
           function(object, colBy, how = 0, gate = 0, top = 0) standardGeneric("pipeFilter")
)

#' @describeIn pipeFilter Method to filter \code{ExprsPipeline} objects.
#' @importFrom stats median quantile
#' @export
setMethod("pipeFilter", "ExprsPipeline",
          function(object, colBy, how, gate, top){

            # Check if ExprsPipeline contains all colBy
            if(!all(colBy %in% colnames(object@summary))){

              stop("Uh oh! This ExprsPipeline does not contain all provided 'colBy' accuracy measures.")
            }

            # Calculate emergent top accuracy measure as product of colBy columns
            accMeasures <- apply(object@summary[colBy], MARGIN = 1, prod)

            # Impose threshold filter if how != 0
            if(how != 0){

              # Set threshold equal to 'how'
              if(0 < how & how <= 1){ threshold <- how

              # Set threshold equal to 'how' as a quantile
              }else if(1 < how & how <= 100){ threshold <- quantile(accMeasures, how / 100)

              }else if(how == "midrange"){ threshold <- (max(accMeasures) + min(accMeasures)) / 2

              }else if(how == "median"){ threshold <- median(accMeasures)

              }else if(how == "mean"){ threshold <- mean(accMeasures)

              }else{ stop("Uh oh! Selected 'how' not recognized!")}

              # Filter ExprsPipeline object
              object@summary <- object@summary[accMeasures >= threshold, ]
              object@machs <- object@machs[accMeasures >= threshold]

              # Filter top accuracy measure as product of colBy columns
              accMeasures <- accMeasures[accMeasures >= threshold]
            }

            # Impose ceiling filter if gate != 0
            if(gate != 0){

              # Set ceiling equal to 'gate'
              if(0 < gate & gate <= 1){ ceiling <- gate

              # Set ceiling equal to 'gate' as a quantile
              }else if(1 < gate & gate <= 100){ ceiling <- quantile(accMeasures, gate / 100)

              }else if(gate == "midrange"){ ceiling <- (max(accMeasures) + min(accMeasures)) / 2

              }else if(gate == "median"){ ceiling <- median(accMeasures)

              }else if(gate == "mean"){ ceiling <- mean(accMeasures)

              }else{ stop("Uh oh! Selected 'gate' not recognized!")}

              # Filter ExprsPipeline object
              object@summary <- object@summary[accMeasures <= ceiling, ]
              object@machs <- object@machs[accMeasures <= ceiling]

              # Filter top accuracy measure as product of colBy columns
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

                         if(top == 0) top <- sum(object@summary$boot == boot)

                         # Order 'top' accMeasures for this boot
                         topMachs <- rev(order(accMeasures[object@summary$boot == boot]))[1:top]

                         # Index by rowname for this boot
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

              if(top == 0) top <- nrow(object@summary)

              # Order 'top' accMeasures for entire object
              topMachs <- order(accMeasures, decreasing = TRUE)[1:top]

              # Index by rowname
              index <- rownames(object@summary)[topMachs]
            }

            # Filter ExprsPipeline object
            final <- rownames(object@summary) %in% index
            object@summary <- object@summary[final, ]
            object@machs <- object@machs[final]

            return(object)
          }
)

###########################################################
### Rename ExprsPipeline "boot" column

#' Rename "boot" Column
#'
#' \code{pipeUnboot} renames the "boot" column of an \code{ExprsPipeline} object
#'  to "unboot".
#'
#' This method provides a convenient adjunct to \code{\link{pipeFilter}} owing to
#'  how \code{exprso} handles \code{ExprsPipeline} objects with a "boot" column.
#'
#' @param object An \code{\link{ExprsPipeline-class}} object.
#' @return An \code{\link{ExprsPipeline-class}} object.
#'
#' @seealso
#' \code{\link{pipeFilter}}\cr
#' \code{\link{pipeUnboot}}\cr
#' \code{\link{plCV}}\cr
#' \code{\link{plGrid}}\cr
#' \code{\link{plGridMulti}}\cr
#' \code{\link{plMonteCarlo}}\cr
#' \code{\link{plNested}}
#'
#' @export
setGeneric("pipeUnboot",
           function(object) standardGeneric("pipeUnboot")
)

#' @describeIn pipeUnboot Method to rename \code{ExprsPipeline} object columns.
#' @export
setMethod("pipeUnboot", "ExprsPipeline",
          function(object){

            if("boot" %in% colnames(object@summary)){

              # Rename 'boot' column to 'unboot'
              colnames(object@summary)[colnames(object@summary) == "boot"] <- "unboot"
            }

            return(object)
          }
)
