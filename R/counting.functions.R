##################################
####### Get all count data


#' Get all the counts on a per sample basis, collapsing the single coordinates
#' 
#' 
setGeneric("get.counts.collapsed", function(object, ...) standardGeneric("get.counts.collapsed"))
setMethod("get.counts.collapsed",
          signature = "Iris",
          definition = function(object){
              combined <- sapply(object@counts,colSums,na.rm=T)
              nams <- object@markers
              if (class(combined) != 'matrix'){
                  counter <- rep(0,length(nams))
                  names(counter) <- nams
                  counts <- sapply(combined,extractCountsF,counter)
              }else{
                  counts <- combined
                  rownames(counts) <- nams
              }
              counts <- counts[order(rownames(counts)),]
              
              return(counts)
})


#' Get all the counts on a per mm2 basis non-collapsed
#' @export
#' 
#' 
setGeneric("get.counts.per.mm2.noncollapsed", function(object, ...) standardGeneric("get.counts.per.mm2.noncollapsed"))
setMethod("get.counts.per.mm2.noncollapsed",
          signature = "Iris",
          definition = function(object){
              
              sizes <- sapply(object@samples,function(x)sapply(x@coordinates,function(y)y@size_in_px))
              
              #counts per mm2
              counts <- lapply(object@counts,function(x,y) x/(y@microns_per_pixel^2),object)
              samps <- colnames(sizes)
              counts <- lapply(samps,function(x,counts,sizes)1000000*sweep(counts[[x]],1,sizes[,x],'/'),counts,sizes)
              names(counts) <- samps
              
              return(counts)
          })


#' Get all the counts on a per mm2 basis
#' @export
#' 
#' 
setGeneric("get.counts.per.mm2", function(object, ...) standardGeneric("get.counts.per.mm2"))
setMethod("get.counts.per.mm2",
          signature = "Iris",
          definition = function(object, digits=2){
              counts <- get.counts.per.mm2.noncollapsed(object)
              if (length(object@counts)>1){
                  means <- sapply(counts,colMeans,na.rm=T)
                  se <- sapply(counts,function(x)apply(x,2,function(y)sd(y,na.rm = T)/sqrt(length(y[!is.na(y)]))))
                  res <- means
                  for (i in 1:ncol(means)){
                     res [,i] <- paste(format(means[,i],digits=digits),'+/-',format(se[,i],digits=digits))   
                  }
              }else{
                  res <- format(counts,digits=digits)
              }
              return(res)
          })


#' Get ratio of counts between two markers
#' @export
#' 
#' 
setGeneric("get.count.ratios", function(object, ...) standardGeneric("get.count.ratios"))
setMethod("get.count.ratios",
          signature = "Iris",
          definition = function(object, marker1, marker2, digits=2){
              ratios <- sapply(object@counts,function(x,m1,m2)x[,m1]/x[,m2],marker1,marker2)
              ratios[is.infinite(ratios)] <- NA
              means <- apply(ratios,2,mean,na.rm=T)
              se <- apply(ratios,2,function(x)sd(x,na.rm = T)/sqrt(length(x[!is.na(x)])))
              res <- paste(format(means,digits=digits),'+/-',format(se,digits=digits))   
              names(res) <- sapply(object@samples,function(x)x@sample_name)
              return(res)
 })


setGeneric("extract.counts", function(object, ...) standardGeneric("extract.counts"))
setMethod("extract.counts",
          signature = "Iris",
          definition = function(object){
              counts <- lapply(object@samples,extract.counts.sample)
              nams <- sort(unique(unlist(lapply(counts,colnames))))
              for (i in 1:length(counts)){
                  counts[[i]] <- counts[[i]][,match(nams,colnames(counts[[i]]))]
                  colnames(counts[[i]]) <- nams
                  counts[[i]][is.na(counts[[i]])] <- 0
              }
              object@counts <- counts
              object@markers <- nams
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

setGeneric("get.counts.noncollapsed", function(object, ...) standardGeneric("get.counts.noncollapsed"))
setMethod("get.counts.noncollapsed",
          signature = "Iris",
          definition = function(object){
              counts <- object@counts
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

