###########################################################
### Import data from file or eSet

#' Import Data from File
#'
#' A convenience function that builds an \code{ExprsArray} object from a tab-delimited file.
#'
#' \code{arrayRead} helps build an \code{ExprsArray} object from a tab-delimited
#'  data file, passing along the \code{file} and \code{...} argument(s) to
#'  \code{\link{read.delim}}. This function expects that the delimited data file has
#'  the following format: rows indicate subject entries while columns indicate measured
#'  variables.
#'
#' The first several columns should contain annotation information (e.g., age,
#'  sex, diagnosis). The remaining columns should contain feature data (e.g. expression
#'  values). The argument \code{probes.begin} defines the j-th column at which the feature
#'  data starts. By default, \code{arrayRead} forces \code{stringsAsFactors = FASE}.
#'
#' This function automatically removes any features with \code{NA} values.
#'
#' @param file A character string. Argument passed along to \code{read.delim}.
#' @param probes.begin A numeric scalar. The j-th column at which feature data starts.
#' @param colID A numeric or character index. The column used to name subjects.
#' @param colBy A numeric or character index. The column that contains group annotations.
#' @param include A list of character vectors. Specifies which annotations in \code{colBy}
#'  to include into which groups. Each element of a list specifies a unique group while
#'  each element of the character vector specifies an annotation to fit to that group. For
#'  binary classification, the first list element defines the negative or control group.
#' @param ... Additional arguments passed along to \code{read.delim}.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{arrayEset}}, \code{\link{GSE2eSet}}
#' @importFrom utils read.delim
#' @export
arrayRead <- function(file, probes.begin, colID, colBy, include, ...){ # args to read.delim

  if(!class(include) == "list") stop("Uh oh! User must provide 'include' argument as list!")
  if(length(include) < 2) stop("Uh oh! User must provide at least two classes!")

  # Pass any additional arguments to read.delim
  table <- read.delim(file, stringsAsFactors = FALSE, ...)

  # Label table by proper subject ID names
  rownames(table) <- make.names(table[, colID], unique = TRUE)

  # Separate the expression data from annotation data
  array <- new(ifelse(length(include) == 2, "ExprsBinary", "ExprsMulti"),
               exprs = t(table[, probes.begin:ncol(table)]),
               annot = table[, 1:(probes.begin-1)],
               preFilter = NULL,
               reductionModel = NULL)

  # Filter out samples in @annot not in 'include'
  array@annot <- array@annot[array@annot[, colBy] %in% unlist(include), ]

  # Number classes in order of 'include'
  for(i in 1:length(include)){

    array@annot[array@annot[, colBy] %in% include[[i]], "defineCase"] <- i
  }

  # Name classes if binary
  if(class(array) == "ExprsBinary"){

    array@annot$defineCase <- ifelse(array@annot$defineCase == 1, "Control", "Case")
  }

  # Set factor if multi
  if(class(array) == "ExprsMulti"){

    array@annot$defineCase <- factor(array@annot$defineCase)
  }

  # Filter @exprs
  array@exprs <- array@exprs[, rownames(array@annot)]

  # Remove probes with missing values
  if(any(is.na(array@exprs))){

    cat("Removing probes with missing values...\n")
    array@exprs <- array@exprs[apply(array@exprs, 1, function(x) !any(is.na(x))), ]
  }

  return(array)
}

#' Import Data from eSet
#'
#' A convenience function that builds an \code{ExprsArray} object from an \code{eSet} object.
#'
#' The package Biobase maintains a popular class object called \code{ExpressionSet} that
#'  often gets used to store expression data. This function converts this \code{eSet}
#'  object into an \code{ExprsArray} object. This function automatically removes any
#'  features with \code{NA} values.
#'
#' @param eSet An \code{eSet} object.
#' @inheritParams arrayRead
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{arrayRead}}, \code{\link{GSE2eSet}}
#' @import Biobase
#' @export
arrayEset <- function(eSet, colBy, include){

  if(!class(include) == "list") stop("Uh oh! User must provide 'include' argument as list!")
  if(length(include) < 2) stop("Uh oh! User must provide at least two classes!")

  # Build an ExprsArray object from the eSet object
  array <- new(ifelse(length(include) == 2, "ExprsBinary", "ExprsMulti"),
               exprs = exprs(eSet),
               annot = eSet@phenoData@data,
               preFilter = NULL,
               reductionModel = NULL)

  # Force @annot rownames to mirror proper @exprs colnames
  colnames(array@exprs) <- make.names(colnames(array@exprs), unique = TRUE)
  rownames(array@annot) <- colnames(array@exprs)

  # Filter out samples in @annot not in 'include'
  array@annot <- array@annot[array@annot[, colBy] %in% unlist(include), ]

  # Number classes in order of 'include'
  for(i in 1:length(include)){

    array@annot[array@annot[, colBy] %in% include[[i]], "defineCase"] <- i
  }

  # Name classes if binary
  if(class(array) == "ExprsBinary"){

    array@annot$defineCase <- ifelse(array@annot$defineCase == 1, "Control", "Case")
  }

  # Set factor if multi
  if(class(array) == "ExprsMulti"){

    array@annot$defineCase <- factor(array@annot$defineCase)
  }

  # Filter @exprs
  array@exprs <- array@exprs[, rownames(array@annot)]

  # Remove probes with missing values
  if(any(is.na(array@exprs))){

    cat("Removing probes with missing values...\n")
    array@exprs <- array@exprs[apply(array@exprs, 1, function(x) !any(is.na(x))), ]
  }

  return(array)
}

#' Convert GSE to eSet
#'
#' A convenience function that builds an \code{eSet} object from a GSE data source.
#'
#' The NCBI GEO hosts files in GSE or GDS format, the latter of which exists as a curated version
#'  the former. These GDS data files easily convert to an  \code{ExpressionSet} (abbreviated
#'  \code{eSet}) object using the \code{GDS2eSet} function available from the GEOquery package.
#'  However, not all GSE data files have a corresponding GDS data file available. To convert GSE
#'  data files into \code{eSet} objects, \code{exprso} provides this convenience function.
#'
#' However, the user should note that GSE data files do not always get stored in an easy to parse format.
#'  Although this function has worked succesfully with some GSE data files, we cannot make any
#'  guarantee that it will work for all GSE data files.
#'
#' To acquire GSE data files, use the function \code{getGEO} from the GEOquery package (e.g.,
#'  \code{getGEO("GSExxxxx", GSEMatrix = FALSE)}). For more information, see the GEOquery package.
#'
#' @param gse A GSE data object retrieved using GEOquery.
#' @param colID A character string. The GSE column name that contains the feature identity.
#' @param colVal A character string. The GSE column name that contains the feature value.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{arrayRead}}, \code{\link{arrayEset}}
#' @import GEOquery Biobase plyr
#' @export
GSE2eSet <- function(gse, colID = "ID_REF", colVal = "VALUE"){

  # Check for non-unique platforms
  gsms <- unlist(lapply(GEOquery::GSMList(gse), function(g){ GEOquery::Meta(g)$platform}))
  if(length(unique(gsms)) > 1) stop("GSE contains non-unique platforms!")

  # Provide an opportunity for user to select a new platform ID column
  if(is.null(colID)){

    cat("The columns available for platform ID include:\n")
    print(Columns(GEOquery::GSMList(gse)[[1]]))
    cat("\n")
    colID <- readline(prompt = "Which column will you use for platform ID?\n")
  }

  # Provide an opportunity for user to select a new platform VALUE column
  if(is.null(colVal)){

    cat("The columns available for platform VALUE include:\n")
    print(Columns(GEOquery::GSMList(gse)[[1]]))
    cat("\n")
    colVal <- readline(prompt = "Which column will you use for platform VALUE?\n")
  }

  # Establish a template for all probes
  probesets <- Table(GEOquery::GPLList(gse)[[1]])$ID

  # Retrieve probe values for each sample
  vals <- lapply(GEOquery::GSMList(gse),
                 function(g){

                   as.numeric(GEOquery::Table(g)[match(probesets, GEOquery::Table(g)[, colID]), colVal])
                 }
  )

  # Combine samples into data.frame
  exprs <- data.frame(vals, row.names = probesets)

  # Retrieve annotations for each sample
  pdata <- lapply(GEOquery::GSMList(gse),
                  function(g){

                    # Retrieve all columns containing "characteristics_ch" in name
                    annotIndex <- grepl("characteristics_ch", names(GEOquery::Meta(g)))
                    annotSlots <- names(GEOquery::Meta(g))[annotIndex]

                    # If there is comma delimitation and multiple annotations
                    if(any(grepl(",", unlist(GEOquery::Meta(g)[annotSlots])) &
                           unlist(lapply(GEOquery::Meta(g)[annotSlots], function(x) length(x) == 1)))){

                      # Split comma containing annotations
                      characteristics <- unlist(strsplit(unlist(GEOquery::Meta(g)[annotSlots]), ",\\s*"))

                    }else{

                      # Pass along annotations unsplit
                      characteristics <- unlist(GEOquery::Meta(g)[annotSlots])
                    }

                    # Retrieve substance of "characteristics" without characteristic name
                    annots <- lapply(strsplit(characteristics, ":\\s*"), function(split) split[2])

                    # Prepare annotations in AnnotatedDataFrame format
                    annots <- data.frame(annots, row.names = GEOquery::Meta(g)$geo_accession)

                    # Rename columns based on "characteristics" name
                    colnames(annots) <- unlist(lapply(strsplit(characteristics, ":\\s*"), function(x) x[1]))

                    # Add a "sample" column
                    cbind(data.frame("sample" = GEOquery::Meta(g)$geo_accession), annots)
                  }
  )

  # Combine annotations into data.frame
  phenoData <- do.call(plyr::rbind.fill, pdata)
  rownames(phenoData) <- phenoData$sample

  # Build eSet object
  eset <- new("ExpressionSet", exprs = as.matrix(exprs), phenoData = as(phenoData, "AnnotatedDataFrame"))

  return(eset)
}
