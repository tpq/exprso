# # NOTE: supplied 'probes' argument does not change @exprs composition
# # NOTE: use "binary" distance when clustering SNP data
# clustHclust <- function(object, probes = 0, col.clustBy = "defineCase", set.include = "Case", k = 2, dist = "euclidean", ...){ # args to hclust
#
#   # Convert 'numeric' probe argument to 'character' probe vector
#   if(class(probes) == "numeric"){
#
#     if(probes == 0) probes <- nrow(object@exprs)
#     probes <- rownames(object@exprs[1:probes, ])
#     data <- t(object@exprs[probes, ])
#   }
#
#   # Build data using supplied 'character' probe vector
#   if(class(probes) == "character"){
#
#     data <- t(object@exprs[probes, ])
#   }
#
#   # Subset 'data' by 'col.clustBy'
#   data <- data[object@annot[rownames(data), col.clustBy] %in% set.include, ]
#
#   # Cluster subsetted data
#   d <- dist(data, method = dist)
#   hc <- hclust(d, ...)
#
#   # Plot dendrogram
#   layout(matrix(c(1), 1, 1, byrow = TRUE))
#   plot(as.dendrogram(hc))
#
#   # Cut into k clusters
#   clusts <- cutree(hc, k = k)
#
#   # Add cluster membership to $cluster
#   object@annot[names(clusts), "cluster"] <- clusts
#   object@annot[!rownames(object@annot) %in% names(clusts), "cluster"] <- 0
#
#   return(object)
# }

#' Cluster Subjects
#'
#' This method clusters subjects based on expression data using any of the
#'  several available clustering algorithms.
#'
#' @param object An \code{ExprsArray object}. The object containing
#'  the subject data to cluster.
#' @param probes A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should guide clustering. A character vector
#'  indicates specifically which features by name should guide clustering.
#'  Set \code{probes = 0} to include all features.
#' @param how A character string. The name of the function used to cluster.
#'  Select from "hcluster", "kmeans", "agnes", "clara", "diana", "fanny",
#'  "pam", "som", "Mclust", or "sota".
#' @param onlyCluster A logical sclar. Toggles whether to return a processed
#'  cluster object or an updated \code{ExprsArray} object.
#' @param ... Additional arguments to the cluster function and/or
#'  other functions used for clustering (e.g., \code{dist} or
#'  \code{cutree}).
#'
#' @importFrom stats hclust kmeans
#' @importFrom cluster agnes clara diana fanny pam
#' @importFrom kohonen som
#' @importFrom mclust Mclust
#' @importFrom clValid sota
#'
#' @export
modCluster <- function(object, probes, how, onlyCluster = FALSE, ...){

  args <- as.list(substitute(list(...)))[-1]

  # Convert 'numeric' probe argument to 'character' probe vector
  if(class(probes) == "numeric"){

    if(probes == 0) probes <- nrow(object@exprs)
    probes <- rownames(object@exprs[1:probes, ])
    data <- t(object@exprs[probes, ])
  }

  # Build data using supplied 'character' probe vector
  if(class(probes) == "character"){

    data <- t(object@exprs[probes, ])
  }

  # Set the argument 'k' aside to use in 'cutree'
  # NOTE: hclust, agnes, and daisy use 'k' via 'cutree'
  if(how %in% c("hclust", "agnes", "daisy")){

    if(!"k" %in% names(args)){

      cat("Setting 'k' to '2' (default behavior, override explicitly)...\n")
      args <- append(args, list("k" = 2))
    }

    # Store 'k' outside of 'args' list, then remove 'k' argument
    arg.k <- args$k
    args <- args[!"k" %in% names(args)]

    if(how == "hclust"){

      # Add 'metric' to args if not already included (i.e., for dist())
      if(!"metric" %in% names(args)){

        cat("Setting 'metric' to 'euclidean' (default behavior, override explicitly)...\n")
        args <- append(args, list("metric" = "euclidean"))
      }

      # Add 'p' to args if not already included (i.e., for dist())
      if(!"p" %in% names(args)){

        cat("Setting 'p' to '2' (default behavior, override explicitly)...\n")
        args <- append(args, list("p" = 2))
      }

      # Build distance matrix, then remove dist() args
      cat("Calculating dissimilarity matrix using 'metric' and 'p' arguments...")
      d <- dist(data, method = args$metric, p = args$p)
      args <- args[!"metric" == names(args)]
      args <- args[!"p" == names(args)]

      args <- append(args, list("d" = d))
      result <- do.call(how, args)

    }else{

      args <- append(args, list("x" = data))
      result <- do.call(how, args)
    }
  }

  if(how == "kmeans"){

    # Add 'centers' to args (instead of 'k')
    if(!"centers" %in% names(args)){

      cat("Setting 'centers' to '2' (default behavior, override explicitly)...\n")
      args <- append(args, list("centers" = 2))
    }

    args <- append(args, list("x" = data))
    result <- do.call(how, args)
  }

  if(how %in% c("clara", "fanny")){

    # Add 'k' to args if not already included
    if(!"k" %in% names(args)){

      cat("Setting 'k' to '2' (default behavior, override explicitly)...\n")
      args <- append(args, list("k" = 2))
    }

    # clara: cluster large applications
    # fanny: fuzzy analysis
    args <- append(args, list("x" = data))
    result <- do.call(how, args)
  }

  if(how == "som"){

    # som: self-organizing map
    # From Wikipedia: The usual arrangement of nodes is a two-dimensional
    #  regular spacing in a hexagonal or rectangular grid. The self-organizing
    #  map describes a mapping from a higher-dimensional input space to a
    #  lower-dimensional map space. The procedure for placing a vector from
    #  data space onto the map is to find the node with the closest (smallest
    #  distance metric) weight vector to the data space vector.
    args <- append(args, list("data" = data))
    result <- do.call(how, args)
  }

  if(how == "Mclust"){

    # mclust: normal mixture modeling
    args <- append(args, list("data" = data))
    result <- do.call(how, args)
  }

  if(how == "sota"){

    # Add 'maxCycles' to args if not already included
    if(!"maxCycles" %in% names(args)){

      cat("Setting 'maxCycles' to '10' (default behavior, override explicitly)...\n")
      args <- append(args, list("maxCycles" = 10))
    }

    # sota: self-organizing tree
    args <- append(args, list("data" = data))
    result <- do.call(how, args)
  }

  if(onlyCluster){

    cat("Returning the processed cluster object...\n")
    return(result)
  }

  cat("Updating ExprsArray object...\n")

  # [[PLACEHOLDER]]
}
