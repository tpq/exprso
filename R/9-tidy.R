#' Build args List
getArgs <- function(...){

  args <- as.list(substitute(list(...)))[-1]
  return(args)
}

#' Set Default args Value
defaultArg <- function(what, as, args){

  if(!what %in% names(args)){

    cat("Setting", what, "to", as.character(as), "(default behavior, override explicitly)...\n")
    names(as) <- what
    args <- append(args, as.list(as))
  }

  return(args)
}

#' Force args Value
forceArg <- function(what, as, args){

  if(!what %in% names(args)){

    cat("Setting", what, "to", as.character(as), "(default behavior, override explicitly)...\n")
    names(as) <- what
    args <- append(args, as.list(as))

  }else{

    if(args[[what]] == as){

      cat(paste0("Uh oh! This function requires ", what, " = ", as,
                 ". Setting ", what, " to ", as, "...\n"))
      args[[what]] <- as
    }
  }

  return(args)
}
