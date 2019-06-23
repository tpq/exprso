#' Combine \code{exprso} Objects
#'
#' \code{conjoin} combines two or more \code{exprso} objects based on their class.
#'
#' When joining \code{ExprsArray} objects, this function returns one
#'  \code{ExprsArray} object as output. This only works on \code{ExprsArray} objects
#'  that have not undergone feature selection. Any missing annotations in \code{@annot}
#'  will get replaced with \code{NA} values.
#'
#' When joining \code{ExprsModel} or \code{ExprsEnsemble} objects,
#'  this function returns an ensemble.
#'
#' When joining \code{ExprsPipeline} objects, this function returns one
#'  \code{ExprsPipeline} object as output. To track which \code{ExprsPipeline}
#'  objects contributed to the resultant object, the source gets flagged
#'  with a \code{boot} column. If a pipeline already has a \code{boot} column,
#'  the original boot tracker will receive an offset (and the old \code{boot}
#'  column will get renamed to \code{unboot}). This system ensures that all
#'  models deriving from the same training set will get handled as a
#'  "pseudo-bootstrap" by downstream \code{\link{pipe}} functions.
#'
#' @param object Any \code{exprso} object.
#' @param ... More objects of the same class.
#' @return See Details.
#' @export
setGeneric("conjoin",
           function(object, ...) standardGeneric("conjoin")
)

#' @describeIn conjoin Method to join \code{ExprsArray} objects.
#' @export
setMethod("conjoin", "ExprsArray",
          function(object, ...){

            # Prepare list of objects
            args <- list(...)
            args <- append(list(object), args)

            if(!lequal(lapply(args, class))){
              stop("All provided objects must have the same class.")
            }

            if(!lequal(lapply(args, function(o) rownames(o@exprs)))){
              stop("All provided objects must have the same class.")
            }

            if(any(!sapply(args, function(e) is.null(e@preFilter) | is.null(e@reductionModel)))){
              stop("All provided objects must not have undergone feature selection.")
            }

            # Prepare single matrix for @exprs and @annot each
            exprs <- as.matrix(do.call(cbind, lapply(args, function(a) a@exprs)))
            annot <- do.call(plyr::rbind.fill, lapply(args, function(a) a@annot))
            rownames(annot) <- unlist(lapply(args, function(a) rownames(a@annot)))

            # Return single ExprsArray object
            new(class(object), exprs = exprs, annot = annot,
                preFilter = NULL, reductionModel = NULL)
          }
)

#' @describeIn conjoin Method to join \code{ExprsModel} objects.
#' @export
setMethod("conjoin", "ExprsModel",
          function(object, ...){

            # Prepare list of objects
            args <- list(...)
            args <- append(list(object), args)

            if(!lequal(lapply(args, class))){
              stop("All provided objects must have the same class.")
            }

            # Return single ExprsEnsemble object
            new("ExprsEnsemble", machs = args)
          }
)

#' @describeIn conjoin Method to join \code{ExprsPipeline} objects.
#' @export
setMethod("conjoin", "ExprsPipeline",
          function(object, ...){

            # Prepare list of objects
            args <- list(...)
            args <- append(list(object), args)

            if(!lequal(lapply(args, class))){
              stop("All provided objects must have the same class.")
            }

            # Get @summary and @machs from list
            args.summary <- lapply(args, function(pl) pl@summary)
            args.machs <- lapply(args, function(pl) pl@machs)

            # When joining pl objects, treat each input as a new "boot"
            b <- 1 # <-- the conjoin boot counter
            pls <- lapply(args.summary,
                          function(pl){

                            pl <- cbind("join" = 0, pl) # <-- set "join" as first column
                            if(!"boot" %in% colnames(pl)){

                              # Add conjoin boot counter
                              pl$join <- b
                              b <<- b + 1

                            }else{

                              # For each boot in $boot
                              for(i in 1:length(unique(pl$boot))){

                                # Change each unique boot to conjoin boot counter
                                pl$join[pl$boot == i] <- b
                                b <<- b + 1
                              }

                              # Rename $boot to $unboot
                              colnames(pl)[colnames(pl) == "boot"] <- "unboot"
                            }

                            # Rename $join to $boot
                            colnames(pl)[colnames(pl) == "join"] <- "boot"

                            return(pl)
                          }
            )

            # Return single ExprsPipeline object
            new("ExprsPipeline", summary = do.call(plyr::rbind.fill, pls), machs = unlist(args.machs))
          }
)

#' @describeIn conjoin Method to join \code{ExprsEnsemble} objects.
#' @export
setMethod("conjoin", "ExprsEnsemble",
          function(object, ...){

            # Prepare list of objects
            args <- list(...)
            args <- append(list(object), args)

            if(!lequal(lapply(args, class))){
              stop("All provided objects must have the same class.")
            }

            # Return single ExprsEnsemble object
            machs <- unlist(lapply(args, function(pl) pl@machs), recursive = TRUE)
            new("ExprsEnsemble", machs = as.list(machs))
          }
)
