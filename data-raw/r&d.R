# Use this script to help make new fs and build methods
data(array)
data(arrayMulti)

array2data <- function(object, top){

  if(class(top) == "numeric"){

    if(length(top) == 1){

      if(top > nrow(object@exprs)) top <- 0
      if(top == 0) top <- nrow(object@exprs)
      top <- rownames(object@exprs[1:top, ])

    }else{

      top <- rownames(object@exprs[top, ])
    }
  }

  t(object@exprs[top, ])
}

labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
data <- array2data(array)
dataMulti <- array2data(arrayMulti)

# Use for new build method:
uniqueFx <- function(data, labels, ...){


}

# Use for new fs method:
uniqueFx <- function(data, top, ...){


}
