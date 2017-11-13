###########################################################
### Import data from data.frame, eSet, or file

#' Import Data as ExprsArray
#'
#' A convenience function that builds an \code{ExprsArray} object.
#'
#' Importing a \code{data.frame} object:
#'
#' This function expects that the imported \code{data.frame} has the following format:
#'  rows indicate subject entries while columns indicate measured variables.
#'  The first several columns should contain annotation information (e.g., age, sex, diagnosis).
#'  The remaining columns should contain feature data (e.g., expression values).
#'  The argument \code{begin} defines the j-th column at which the feature
#'  data starts. This function automatically removes any features with \code{NA} values.
#'  Take care to remove any \code{factor} columns before importing.
#'
#' Importing an \code{ExpressionSet} object:
#'
#' The package Biobase maintains a popular class object called \code{ExpressionSet} that
#'  often gets used to store expression data. This function converts this \code{eSet}
#'  object into an \code{ExprsArray} object. This function automatically removes any
#'  features with \code{NA} values.
#'
#' Importing a \code{file}:
#'
#' \code{arrayExprs} can also build an \code{ExprsArray} object from a tab-delimited
#'  data file, passing along the \code{file} and \code{...} argument(s) to
#'  \code{\link{read.delim}}. All rules for \code{data.frame} import also apply here.
#'  By default, \code{arrayExprs} forces \code{stringsAsFactors = FASE}.
#'
#' @param object What to import as an \code{ExprsArray} object. See Details.
#' @param colBy A numeric or character index. The column that contains group annotations.
#' @param include A list of character vectors. Specifies which annotations in \code{colBy}
#'  to include in which groups. Each element of the list specifies a unique group while
#'  each element of the character vector specifies which annotations define that group. For
#'  binary classification, the first list element defines the negative, or control, group.
#' @param colID A numeric or character index. The column used to name subjects.
#'  For \code{data.frame} or file import only.
#' @param begin A numeric scalar. The j-th column at which feature data starts.
#'  For \code{data.frame} or file import only.
#' @param ... Additional arguments passed along to \code{read.delim}.
#'  For file import only.
#'
#' @return An \code{ExprsArray} object.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{GSE2eSet}}
#' @importFrom utils read.delim
#' @export
arrayExprs <- function(object, colBy, include, colID, begin, ...){

  if(missing(colBy)) stop("Uh oh! User must specify 'colBy' argument!")
  if(missing(include)) stop("Uh oh! User must specify 'include' argument!")
  if(!class(include) == "list") stop("Uh oh! User must provide 'include' argument as list!")
  if(length(include) < 2) stop("Uh oh! User must provide at least two classes!")

  if(class(object) == "data.frame"){

    if(missing(colID)) stop("Uh oh! User must specify 'colID' argument!")
    if(missing(begin)) stop("Uh oh! User must specify 'begin' argument!")
    rownames(object) <- object[, colID]
    exprs <- exprs <- t(object[, begin:ncol(object)])
    annot <- object[, 1:(begin-1)]

  }else if(class(object) == "ExpressionSet"){

    if(!requireNamespace("Biobase", quietly = TRUE)){
      stop("Biobase needed for this function to work. Please install it.",
           call. = FALSE)
    }

    exprs <- Biobase::exprs(object)
    annot <- object@phenoData@data

  }else if(class(object) == "character" & file.exists(object)){

    if(missing(colID)) stop("Uh oh! User must specify 'colID' argument!")
    if(missing(begin)) stop("Uh oh! User must specify 'begin' argument!")
    args <- getArgs(...)
    args <- forceArg("stringsAsFactors", FALSE, args)
    args <- append(args, list("file" = object))
    object <- do.call("read.delim", args)
    rownames(object) <- object[, colID]
    exprs <- t(object[, begin:ncol(object)])
    annot <- object[, 1:(begin-1)]

  }else{

    stop("Uh oh! No default method for importing this object as an ExprsArray.")
  }

  newClass <- ifelse(length(include) == 2, "ExprsBinary", "ExprsMulti")
  array <- new(newClass, exprs = exprs, annot = annot, preFilter = NULL, reductionModel = NULL)

  # Force @annot rownames to mirror proper @exprs colnames
  colnames(array@exprs) <- make.names(colnames(array@exprs), unique = TRUE)
  rownames(array@annot) <- colnames(array@exprs)

  # Use 'include' to filter subjects and label classes
  array@annot <- array@annot[array@annot[, colBy] %in% unlist(include), ]
  array@exprs <- array@exprs[, rownames(array@annot)]
  for(i in 1:length(include)){

    array@annot[array@annot[, colBy] %in% include[[i]], "defineCase"] <- i
  }

  if(newClass == "ExprsBinary"){

    array@annot$defineCase <- ifelse(array@annot$defineCase == 1, "Control", "Case")

  }else if(newClass == "ExprsMulti"){

    array@annot$defineCase <- factor(array@annot$defineCase)

  } # if regrso, do nothing here

  # Remove features with missing values
  if(any(is.na(array@exprs))){

    cat("Removing features with missing values...\n")
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
#' @param colBy A character string. The GSE column name that contains the feature value.
#'  If missing, function will prompt user for a column name after previewing options.
#' @param colID A character string. The GSE column name that contains the feature identity.
#'  If missing, function will prompt user for a column name after previewing options.
#'
#' @return An \code{ExpressionSet} object.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{arrayExprs}}
#' @export
GSE2eSet <- function(gse, colBy, colID){

  if(!requireNamespace("Biobase", quietly = TRUE)){
    stop("Biobase needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(!requireNamespace("GEOquery", quietly = TRUE)){
    stop("GEOquery needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Check for non-unique platforms
  gsms <- unlist(lapply(GEOquery::GSMList(gse), function(g){ GEOquery::Meta(g)$platform}))
  if(length(unique(gsms)) > 1) stop("GSE contains non-unique platforms!")

  # Provide an opportunity for user to select a new platform ID column
  if(missing(colID)){

    cat("The columns available for platform ID include:\n")
    print(GEOquery::Columns(GEOquery::GSMList(gse)[[1]]))
    cat("\n")
    colID <- readline(prompt = "Which column (by name) will you use for platform ID? \n")
  }

  # Provide an opportunity for user to select a new platform VALUE column
  if(missing(colBy)){

    cat("The columns available for platform VALUE include:\n")
    print(GEOquery::Columns(GEOquery::GSMList(gse)[[1]]))
    cat("\n")
    colBy <- readline(prompt = "Which column (by name) will you use for platform VALUE? \n")
  }

  # Establish a template for all features
  featsets <- GEOquery::Table(GEOquery::GPLList(gse)[[1]])$ID

  # Retrieve feature values for each sample
  vals <- lapply(GEOquery::GSMList(gse),
                 function(g){

                   as.numeric(GEOquery::Table(g)[match(featsets, GEOquery::Table(g)[, colID]), colBy])
                 }
  )

  # Combine samples into data.frame
  exprs <- data.frame(vals, row.names = featsets)

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
