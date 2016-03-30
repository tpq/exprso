###########################################################
### Import data from file

# NOTE: 'include' expects a list with the first element set to CONTROL
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
  if(class(array) == "ExprsBinary") array@annot$defineCase <- ifelse(array@annot$defineCase == 1, "Control", "Case")
  
  # Set factor if multi
  if(class(array) == "ExprsMulti") array@annot$defineCase <- factor(array@annot$defineCase)
  
  # Filter @exprs
  array@exprs <- array@exprs[, rownames(array@annot)]
  
  return(array)
}

###########################################################
### Import data from GSE

# NOTE: retrieve GSE object with getGEO("GSExxxxx", GSEMatrix = FALSE)
GSE2eSet <- function(gse, colID = "ID_REF", colVal = "VALUE"){
  
  require(GEOquery)
  require(Biobase)
  require(plyr)
  
  # Check for non-unique platforms
  gsms <- unlist(lapply(GSMList(gse), function(g){ Meta(g)$platform}))
  
  # Check for non-unique platforms
  if(length(unique(gsms)) > 1) stop("GSE contains non-unique platforms!")
  
  # Provide an opportunity for user to select a new platform ID column
  if(is.null(colID)){
    
    cat("The columns available for platform ID include:\n")
    print(Columns(GSMList(gse)[[1]]))
    cat("\n")
    colID <- readline(prompt = "Which column will you use for platform ID?\n")
  }
  
  # Provide an opportunity for user to select a new platform VALUE column
  if(is.null(colVal)){
    
    cat("The columns available for platform VALUE include:\n")
    print(Columns(GSMList(gse)[[1]]))
    cat("\n")
    colVal <- readline(prompt = "Which column will you use for platform VALUE?\n")
  }
  
  # Establish a template for all probes
  probesets <- Table(GPLList(gse)[[1]])$ID
  
  # Retrieve probe values for each sample
  vals <- lapply(GSMList(gse), function(g){ as.numeric(Table(g)[match(probesets, Table(g)[, colID]), colVal])})
  
  # Combine samples into data.frame
  exprs <- data.frame(vals, row.names = probesets)
  
  # Retrieve annotations for each sample
  pdata <- lapply(GSMList(gse),
                  function(g){
                    
                    # Retrieve all columns containing "characteristics_ch" in name
                    annotSlots <- names(Meta(g))[grepl("characteristics_ch", names(Meta(g)))]
                    
                    # If there is any comma delimitation
                    grepl(",", unlist(Meta(g)[annotSlots])) &
                      unlist(lapply(Meta(g)[annotSlots], function(x) length(x) == 1))
                    
                    # If there is any comma delimitation
                    if(any(grepl(",", unlist(Meta(g)[annotSlots])) & # IF ANY commas AND...
                           unlist(lapply(Meta(g)[annotSlots], function(x) length(x) == 1)))){ # ...!length > 1
                      
                      # Split comma containing annotations
                      characteristics <- unlist(strsplit(unlist(Meta(g)[annotSlots]), ",\\s*"))
                      
                    }else{
                      
                      # Pass along annotations unsplit
                      characteristics <- unlist(Meta(g)[annotSlots])
                    }
                    
                    # Retrieve substance of "characteristics" without characteristic name
                    annots <- lapply(strsplit(characteristics, ":\\s*"), function(split) split[2])
                    
                    # Prepare annotations in AnnotatedDataFrame format
                    annots <- data.frame(annots, row.names = Meta(g)$geo_accession)
                    
                    # Rename columns based on "characteristics" name
                    colnames(annots) <- unlist(lapply(strsplit(characteristics, ":\\s*"), function(x) x[1]))
                    
                    # Add a "sample" column
                    cbind(data.frame("sample" = Meta(g)$geo_accession), annots)
                  })
  
  # Combine annotations into data.frame
  phenoData <- do.call(rbind.fill, pdata)
  
  # Set row.names to "sample" column
  rownames(phenoData) <- phenoData$sample
  
  # Build eSet object
  eset <- new("ExpressionSet", exprs = as.matrix(exprs), phenoData = as(phenoData, "AnnotatedDataFrame"))
  
  return(eset)
}

# NOTE: 'include' expects a list with the first element set to CONTROL
arrayEset <- function(eSet, colBy, include){
  
  require(Biobase)
  
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
  if(class(array) == "ExprsBinary") array@annot$defineCase <- ifelse(array@annot$defineCase == 1, "Control", "Case")
  
  # Set factor if multi
  if(class(array) == "ExprsMulti") array@annot$defineCase <- factor(array@annot$defineCase)
  
  # Filter @exprs
  array@exprs <- array@exprs[, rownames(array@annot)]
  
  # Remove probes with missing values
  if(any(is.na(array@exprs))){
    
    cat("Removing probes with missing values...\n")
    array@exprs <- array@exprs[apply(array@exprs, 1, function(x) !any(is.na(x))), ]
  }
  
  return(array)
}
