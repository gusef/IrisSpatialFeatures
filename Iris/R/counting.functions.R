##################################
####### Get all count data


setGeneric("get.counts", function(object, ...) standardGeneric("get.counts"))
setMethod("get.counts",
          signature = "Iris",
          definition = function(object){
              return(object@counts)
})


setGeneric("extract.counts", function(object, ...) standardGeneric("extract.counts"))
setMethod("extract.counts",
          signature = "Iris",
          definition = function(object){
              counts <- lapply(object@samples,extract.counts.sample)
              nams <- unique(unlist(lapply(counts,colnames)))
              for (i in 1:length(counts)){
                  counts[[i]] <- counts[[i]][,match(nams,colnames(counts[[i]]))]
              }
              combined <- sapply(counts,colSums,na.rm=T)
              
              if (class(combined) != 'matrix'){
                  nams <- sort(unique(unlist(lapply(combined,names))))
                  counter <- rep(0,length(nams))
                  names(counter) <- nams
                  object@counts <- sapply(combined,extractCountsF,counter)
              }else{
                  object@counts <- combined
                  rownames(object@counts) <- nams
              }
              object@counts <- object@counts[order(rownames(object@counts)),]
              return(object)
})

setGeneric("extract.counts.sample", function(object, ...) standardGeneric("extract.counts.sample"))
setMethod("extract.counts.sample",
          signature = "Sample",
          definition = function(object){
              counts <- lapply(object@coordinates,function(x)table(x@ppp$marks))
              nams <- unique(unlist(lapply(counts,names)))
              counter <- rep(0,length(nams))
              names(counter) <- nams
              counts <- t(sapply(counts,extractCountsF,counter))
              return(counts)
})

#helperfunction to count the features making sure missing celltypes don't cause problems
extractCountsF <- function(x,counter){
    counter[match(names(x),names(counter))] <- x
    return(counter)
}

setGeneric("extract.counts.noncollapsed", function(object, ...) standardGeneric("extract.counts.noncollapsed"))
setMethod("extract.counts.noncollapsed",
          signature = "Sample",
          definition = function(object){
              counts <- lapply(object@samples,extract.counts.sample)
              nams <- unique(unlist(lapply(counts,colnames)))
    
              standardize <- function(x,nams){
                  y <- matrix(0,nrow=nrow(x),ncol=length(nams))
                  colnames(y) <- nams
                  rownames(y) <- rownames(x)
                  y[,colnames(x)] <- x
                  return(y)
              }    
              counts <- lapply(counts,standardize,nams)
    return(counts)
})

