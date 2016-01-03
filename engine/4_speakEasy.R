###########################################################
### speakEasy

# Converts expression values as measured by one PROBEID to expression values of another PROBEID
# NOTE: Function only works on ExprsArray objects that have not undergone feature selection
# NOTE: 'intermediate' argument determines how AnnotationDbi converts between PROBEIDs
# NOTE: IF platformTo = NULL, function returns expression values as 'intermediate' itself
# NOTE: IF fillBlanks = TRUE, this function will set all unused keys to the global median
# NOTE: 'fillBlanks' should equal TRUE when assembling an external validation set
setGeneric("speakEasy",
           function(object, ...){
             standardGeneric("speakEasy")
           }
)

setMethod("speakEasy", "ExprsArray",
          function(object, platformFrom, platformTo = NULL, intermediate = "ENTREZID", fillBlanks = TRUE){
            
            require(AnnotationDbi)
            require(plyr)
            
            if(!is.null(object@preFilter) | !is.null(object@reductionModel)){
              
              stop("This function is not equipped to handle ExprsArray objects that have undergone feature selection!")
            }
            
            if(is.null(platformTo) & fillBlanks){
              
              warning("Without a selected 'platformTo', 'fill in the blanks' depends on 'platformFrom'!")
            }
            
            # Prepare data for PROBEID platform conversion
            data <- data.frame("id" = rownames(object@exprs), object@exprs, stringsAsFactors = FALSE)
            
            # Convert probe vector to intermediate metric
            require(platformFrom, character.only = TRUE)
            key <- select(get(platformFrom), keys = data[, "id"], columns = intermediate, keytype = "PROBEID")
            
            # Merge expression data with intermediate conversion key
            mergedInter <- merge(x = key, y = data, by.x = "PROBEID", by.y = "id")
            
            # Remove NA intermediates
            mergedInter <- mergedInter[!is.na(mergedInter[, intermediate]), ]
            
            # Correlate intermediate metric with ~median~ expression of all corresponding probes
            intermeds <- ldply(unique(mergedInter[, intermediate]),
                               function(id){
                                 
                                 # Subset expression values for single intermediate metric
                                 df <- mergedInter[mergedInter[, intermediate] %in% id, ]
                                 
                                 # Calculate all median expression values
                                 final <- apply(df[, !colnames(df) %in% c("PROBEID", intermediate)], 2, median)
                                 
                                 # Clean up data
                                 final <- data.frame(id, data.frame(as.list(final)), stringsAsFactors = FALSE)
                                 
                                 return(final)
                               })
            
            # convert Entrez IDs to NEW probe vector
            if(!is.null(platformTo)){
              
              # Convert probe vector to intermediate metric
              require(platformTo, character.only = TRUE)
              key <- select(get(platformTo), keys = intermeds[, "id"], columns = "PROBEID", keytype = intermediate)
              
              # Merge expression data with final PROBEID conversion key
              mergedTo <- merge(x = key, y = intermeds, by.x = intermediate, by.y = "id")
              
              # Remove NA PROBEIDs
              mergedTo <- mergedTo[!is.na(mergedTo[, "PROBEID"]), ]
              
              # Correlate final PROBEID with ~median~ expression of all corresponding intermediate metrics
              finals <- ldply(unique(mergedTo[, "PROBEID"]),
                              function(id){
                                
                                # Subset expression values for single intermediate metric
                                df <- mergedTo[mergedTo[, "PROBEID"] %in% id, ]
                                
                                # Calculate all median expression values
                                final <- apply(df[, !colnames(df) %in% c("PROBEID", intermediate)], 2, median)
                                
                                # Clean up data
                                final <- data.frame(id, data.frame(as.list(final)), stringsAsFactors = FALSE)
                                
                                return(final)
                              })
              
              if(fillBlanks){
                
                # Find missing PROBEIDs
                missing <- data.frame("missing" = keys(get(platformTo), keytype = "PROBEID"), stringsAsFactors = FALSE)
                
                # Merge missing PROBEIDs with median transformed results
                mergeBlanks <- merge(x = finals, y = missing, by.x = "id", by.y = "missing", all.y = TRUE)
                
                # Clean up data
                mergeBlanks <- data.frame(mergeBlanks[, !colnames(mergeBlanks) %in% "id"], row.names = mergeBlanks[, "id"])
                
                # Impute missing values with global median
                global <- median(as.matrix(finals[, !colnames(finals) %in% "id"]))
                mergeBlanks[is.na(mergeBlanks[, 1]), ] <- global
                
                # Prepare @exprs output
                exprs <- as.matrix(mergeBlanks)
                
                
              }else{
                
                # Clean up data
                finals <- data.frame(finals[, !colnames(finals) %in% "id"], row.names = finals[, "id"])
                
                # Prepare @exprs output
                exprs <- as.matrix(finals)
              }
              
            }else{
              
              if(fillBlanks){
                
                # Find missing PROBEIDs
                missing <- data.frame("missing" = keys(get(platformFrom), keytype = intermediate), stringsAsFactors = FALSE)
                
                # Merge missing PROBEIDs with median transformed results
                mergeBlanks <- merge(x = intermeds, y = missing, by.x = "id", by.y = "missing", all.y = TRUE)
                
                # Clean up data
                mergeBlanks <- data.frame(mergeBlanks[, !colnames(mergeBlanks) %in% "id"], row.names = mergeBlanks[, "id"])
                
                # Impute missing values with global median
                global <- median(as.matrix(intermeds[, !colnames(intermeds) %in% "id"]))
                mergeBlanks[is.na(mergeBlanks[, 1]), ] <- global
                
                # Prepare @exprs output
                exprs <- as.matrix(mergeBlanks)
                
              }else{
                
                # Clean up data
                intermeds <- data.frame(intermeds[, !colnames(intermeds) %in% "id"], row.names = intermeds[, "id"])
                
                # Prepare @exprs output
                exprs <- as.matrix(intermeds)
                
              }
            }
            
            # For NULL platformTo, record intermediate
            if(is.null(platformTo)) platformTo <- intermediate
            
            # Prepare speakEasy terms for @annot
            terms <- list("platformFrom" = platformFrom, "intermediate" = intermediate, "platformTo" = platformTo, "fillBlanks" = fillBlanks)
            
            # Prepare @annot output
            annot <- cbind(object@annot, terms, stringsAsFactors = FALSE)
            
            # Combine results into ExprsArray object
            array <- new("ExprsArray", exprs = exprs, annot = annot, preFilter = object@preFilter, reductionModel = object@reductionModel)
            
            return(array)
          }
)

# Prunes multiple speakEasy results so that they all have the same probeset
# Necessary because each speakEasy run may use non-overlapping intermediates
# Use when platformTo = NULL or fillBlanks = FALSE during speakEasy call
# NOTE: Function only works on ExprsArray objects that have not undergone feature selection
setGeneric("abridge",
           function(object, ...){
             standardGeneric("abridge")
           }
)

setMethod("abridge", "ExprsArray",
          function(object, ...){ # args to include additional ExprsArray objects
            
            # Prepare list of ExprsArray objects
            args <- list(...)
            index <- unlist(lapply(args, function(arg) inherits(arg, "ExprsArray")))
            args <- append(list(object), args[index])
            
            if(!length(args) > 1 | any(!sapply(args, function(e) identical(class(e), class(object))))){
              
              stop("User must provide additional ExprsArray objects of the same class!")
            }
            
            if(any(!sapply(args, function(e) is.null(e@preFilter) | is.null(e@reductionModel)))){
              
              stop("This function is not equipped to handle objects that have undergone feature selection!")
            }
            
            # IF there is any missing 'platformTo'
            if(any(unlist(lapply(args, function(e) is.null(e@annot$platformTo))))){
              
              stop("All ExprsArray objects must have same 'platformTo'!")
            }
            
            # IF there is not one 'platformTo'
            if(!length(unique(unlist(lapply(args, function(e) e@annot$platformTo)))) == 1){
              
              stop("All ExprsArray objects must have same 'platformTo'!")
            }
            
            # IF any 'fillBlanks' == FALSE
            if(any(!unlist(lapply(args, function(e) e@annot$fillBlanks)))){
              
              warning("User did not 'fill in the blanks' during speakEasy! Expect fewer probes in results.")
            }
            
            # Initialize while-loop probe counter
            i <- 1
            
            # Extract those probes present in all ExprsArray objects
            while(i < length(args)){
              
              if(i == 1){
                
                # Intersect probes from the first with the second
                probes <- intersect(rownames(args[[i]]@exprs), rownames(args[[i + 1]]@exprs))
                
              }else{
                
                # Intersect probes from the prior with the next
                probes <- intersect(probes, rownames(args[[i + 1]]@exprs))
              }
              
              # Increase counter
              i <- i + 1
            }
            
            # Subset each ExprsArray object with intersected probeset
            arrays <- lapply(args,
                             function(array){
                               
                               array@exprs <- array@exprs[match(probes, rownames(array@exprs)), ]
                               array@annot <- cbind(array@annot, list("abridged" = TRUE))
                               
                               return(array)
                             })
            
            return(arrays)
          }
)
