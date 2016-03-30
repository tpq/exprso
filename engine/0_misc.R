###########################################################
### Define generic functions

setGeneric("modFilter", function(object, ...) standardGeneric("modFilter"))
setGeneric("modTransform", function(object, ...) standardGeneric("modTransform"))
setGeneric("modNormalize", function(object, ...) standardGeneric("modNormalize"))

###########################################################
### Pre-process data

# Sets minimum and maximum thresholds, then removes genes with low MAX:MIN relationship
setMethod("modFilter", "ExprsArray",
          function(object, threshold, maximum, beta1, beta2, displayAll = FALSE){
            
            # Prepare layout for multiple plots
            layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))
            
            # For large datasets, plot only a sub-sample
            v <- as.vector(object@exprs)
            if(!displayAll & length(v) > 1000) v <- sample(v, 1000)
            
            qqnorm(v, main = "Normal Q-Q Plot (Before)")
            plot(density(v), main = "Density Plot (Before)")
            boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
            
            # Set minimum and maximum threshold
            object@exprs[object@exprs < threshold] <- threshold
            object@exprs[object@exprs > maximum] <- maximum
            
            # Calculate RANGE and MAX:MIN
            index.beta1 <- apply(object@exprs, MARGIN = 1, function(gene) max(gene) - min(gene))
            index.beta2 <- apply(object@exprs, MARGIN = 1, function(gene) max(gene) / min(gene))
            
            # INCLUDE those with variance GREATER THAN beta1 AND beta2
            object@exprs <- object@exprs[index.beta1 > beta1 & index.beta2 > beta2, ]
            
            # For large datasets, plot only a sub-sample
            v <- as.vector(object@exprs)
            if(!displayAll & length(v) > 1000) v <- sample(v, 1000)
            
            qqnorm(v, main = "Normal Q-Q Plot (After)")
            plot(density(v), main = "Density Plot (After)")
            boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")
            
            return(object)
          }
)

setMethod("modTransform", "ExprsArray",
          function(object, displayAll = FALSE){
            
            # Prepare layout for multiple plots
            layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))
            
            # For large datasets, plot only a sub-sample
            v <- as.vector(object@exprs)
            if(!displayAll & length(v) > 1000) v <- sample(v, 1000)
            
            qqnorm(v, main = "Normal Q-Q Plot (Before)")
            plot(density(v), main = "Density Plot (Before)")
            boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
            
            # Perform log transformation
            object@exprs <- log(object@exprs, base = 2)
            
            # For large datasets, plot only a sub-sample
            v <- as.vector(object@exprs)
            if(!displayAll & length(v) > 1000) v <- sample(v, 1000)
            
            qqnorm(v, main = "Normal Q-Q Plot (After)")
            plot(density(v), main = "Density Plot (After)")
            boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")
            
            return(object)
          }
)

# MARGIN = 1 to normalize by gene vector
# MARGIN = 2 to normalize by subject vector
# MARGIN = c(1, 2) for both
setMethod("modNormalize", "ExprsArray",
          function(object, MARGIN = c(1, 2), displayAll = FALSE){
            
            # Prepare layout for multiple plots
            layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))
            
            # For large datasets, plot only a sub-sample
            v <- as.vector(object@exprs)
            if(!displayAll & length(v) > 1000) v <- sample(v, 1000)
            
            qqnorm(v, main = "Normal Q-Q Plot (Before)")
            plot(density(v), main = "Density Plot (Before)")
            boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
            
            # Perform normalization
            if(2 %in% MARGIN) object@exprs <- apply(object@exprs, MARGIN = 2, function(samp) (samp - mean(samp)) / sd(samp))
            if(1 %in% MARGIN) object@exprs <- t(apply(object@exprs, MARGIN = 1, function(gene) (gene - mean(gene)) / sd(gene)))
            
            # For large datasets, plot only a sub-sample
            v <- as.vector(object@exprs)
            if(!displayAll & length(v) > 1000) v <- sample(v, 1000)
            
            qqnorm(v, main = "Normal Q-Q Plot (After)")
            plot(density(v), main = "Density Plot (After)")
            boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")
            
            return(object)
          }
)

###########################################################
### Compare datasets

# NOTE: Use arraySubset to define new 'array.train' and 'array.valid' objects for pairwise comparisons
# NOTE: Provided 'colBy' determines categorical variable used for "internal" comparisons
# NOTE: This function will always perform "internal" comparisons whether or not 'array.valid' is provided
# NOTE: IF 'array.valid' is provided, function will also compare 'array.train' against 'array.valid'
compare <- function(array.train, array.valid = NULL, colBy = "defineCase", include = c("Case", "Control"), cutoff = .05){
  
  if(!inherits(array.train, "ExprsArray")) stop("Uh oh! You can only compare ExprsArray objects.")
  
  # Subset the first ExprsArray object
  arrays <- list(arraySubset(array.train, colBy = colBy, include = include))
  
  # Subset second object if provided
  if(!is.null(array.valid)){
    
    if(!inherits(array.train, "ExprsArray")) stop("Uh oh! You can only compare ExprsArray objects.")
    arrays <- append(arrays, arraySubset(array.valid, colBy = colBy, include = include))
  }
  
  # Perform ANOVA or chi-square across each provided array
  annots.each <- lapply(arrays,
                        function(array){
                          
                          cat("\nPerforming ANOVA or chi-square analysis for each annotation...\n")
                          
                          # Make sure colBy contains categorical data
                          if(!class(array@annot[, colBy]) %in% c("character", "factor")){
                            
                            cat("Expect strange results when 'colBy' does not contain categorical data!\n")
                          }
                          
                          # Prepare the annotations to test as independent variables
                          annots <- colnames(array@annot)[!colnames(array@annot) %in% colBy]
                          
                          # For each independent variable...
                          index.each <- lapply(annots,
                                               function(annot){
                                                 
                                                 if(class(array@annot[, annot]) %in% c("character", "factor")){
                                                   
                                                   # Build contingency table for one independent variable
                                                   compare <- table(array@annot[, c(annot, colBy)])
                                                   
                                                   # Perform chi-square test and concatenate results
                                                   chisq.p <- chisq.test(compare)$p.value
                                                   cat("\t", "Annotation ", annot, ": ", chisq.p, "\n", sep = "")
                                                   
                                                   # Return statistical significance as boolean
                                                   ifelse(is.na(chisq.p), return(FALSE), return(chisq.p < cutoff))
                                                   
                                                 }else{
                                                   
                                                   # Prepare data for ANOVA
                                                   df <- data.frame("y" = array@annot[, annot], "x" = as.factor(array@annot[, colBy]))
                                                   
                                                   # Check if ANOVA will work after removing NA values
                                                   if(length(unique(na.omit(df)$x)) > 1){
                                                     
                                                     # Fit linear model to independent variable
                                                     fit <- lm(y ~ x, data = df, na.action = na.omit)
                                                     
                                                     # Perform ANOVA test and concatenate results
                                                     anova.p <- anova(fit)$Pr[1]
                                                     cat("\t", "Annotation ", annot, ": ", anova.p, "\n", sep = "")
                                                     
                                                     # Return statistical significance as boolean
                                                     return(anova.p < cutoff)
                                                     
                                                   }else{
                                                     
                                                     # ANOVA fails
                                                     cat("\t", "Annotation ", annot, ": insufficient measurements \n", sep = "")
                                                     return(FALSE)
                                                   }
                                                 }
                                               })
                          
                          # Return statistically significant annotations
                          return(annots[unlist(index.each)])
                        })
  
  if(!is.null(array.valid)){
    
    cat("\nPerforming ANOVA or chi-square analysis between ExprsArray objects for each annotation...\n")
    
    # Prepare the annotations to test as independent variables
    annots <- colnames(array.train@annot)[colnames(array.train@annot) %in% colnames(array.valid@annot)]
    
    # For each independent variable...
    index.both <- lapply(annots,
                         function(annot){
                           
                           # Build a data.frame to use in ExprsArray object comparisons
                           df <- rbind(data.frame("array" = "array.train", "annot" = array.train@annot[, annot]),
                                       data.frame("array" = "array.valid", "annot" = array.valid@annot[, annot]))
                           
                           if(class(array.train@annot[, annot]) %in% c("character", "factor")){
                             
                             # Build contingency table for one independent variable
                             compare <- table(df)
                             
                             # Perform chi-square test and concatenate results
                             chisq.p <- chisq.test(compare)$p.value
                             cat("\t", "Annotation ", annot, ": ", chisq.p, "\n", sep = "")
                             
                             # Return statistical significance as boolean
                             ifelse(is.na(chisq.p), return(FALSE), return(chisq.p < cutoff))
                             
                           }else{
                             
                             # Check if ANOVA will work after removing NA values
                             if(length(unique(na.omit(df)$array)) > 1){
                               
                               # Fit linear model to independent variable
                               fit <- lm(annot ~ array, data = df, na.action = na.omit)
                               
                               # Perform ANOVA test and concatenate results
                               anova.p <- anova(fit)$Pr[1]
                               cat("\t", "Annotation ", annot, ": ", anova.p, "\n", sep = "")
                               
                               # Return statistical significance as boolean
                               return(anova.p < cutoff)
                               
                             }else{
                               
                               # ANOVA fails
                               cat("\t", "Annotation ", annot, ": insufficient measurements \n", sep = "")
                               return(FALSE)
                             }
                           }
                         })
    
    # Return statistically significant annotations
    annots.both <- annots[unlist(index.both)]
    
    # Prepare result
    result <- append(annots.each, list(annots.both))
    names(result) <- c("differences.within.train", "differences.within.valid", "differences.between")
    
  }else{
    
    # Set index.both = NULL if array.valid = NULL
    annots.both <- NULL
    
    # Prepare result
    result <- annots.each
    names(result) <- c("differences.train")
  }
  
  return(result)
}
