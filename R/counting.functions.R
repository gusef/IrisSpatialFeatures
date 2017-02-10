##################################
####### Get all count data


#' Get all the counts on a per sample basis, collapsing the single coordinates
#' 
#' 
setGeneric("get.counts.collapsed", function(x, ...) standardGeneric("get.counts.collapsed"))
setMethod("get.counts.collapsed",
          signature = "Iris",
          definition = function(x){
              combined <- sapply(x@counts,colSums,na.rm=T)
              nams <- x@markers
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
setGeneric("get.counts.per.mm2.noncollapsed", function(x, ...) standardGeneric("get.counts.per.mm2.noncollapsed"))
setMethod("get.counts.per.mm2.noncollapsed",
          signature = "Iris",
          definition = function(x){
              
              sizes <- lapply(x@samples,function(y)sapply(y@coordinates,function(z)z@size_in_px))
              
              #counts per mm2
              counts <- lapply(x@counts,function(y,z) y/(z@microns_per_pixel^2),x)
              samps <- names(sizes)
              counts <- lapply(samps,function(y,counts,sizes)1000000*sweep(counts[[y]],1,sizes[[y]],'/'),counts,sizes)
              names(counts) <- samps
              
              return(counts)
          })


#' Get all the counts on a per mm2 basis
#' @importFrom stats sd
#' @export
#' 
#' 
setGeneric("get.counts.per.mm2", function(x, ...) standardGeneric("get.counts.per.mm2"))
setMethod("get.counts.per.mm2",
          signature = "Iris",
          definition = function(x, digits=2){
              counts <- get.counts.per.mm2.noncollapsed(x)
              if (length(x@counts)>1){
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
#' @importFrom stats sd
#' @export
#' 
#' 
setGeneric("get.count.ratios", function(x, ...) standardGeneric("get.count.ratios"))
setMethod("get.count.ratios",
          signature = "Iris",
          definition = function(x, marker1, marker2, digits=2){
              ratios <- sapply(x@counts,function(x,m1,m2)x[,m1]/x[,m2],marker1,marker2)
              for (idx in 1:length(ratios)){
                 ratios[[idx]][is.infinite(ratios[[idx]])] <- NA
              }
              means <- sapply(ratios,mean,na.rm=T)
              se <- sapply(ratios,function(x)sd(x,na.rm = T)/sqrt(length(x[!is.na(x)])))
              res <- paste(format(means,digits=digits),'+/-',format(se,digits=digits))   
              names(res) <- sapply(x@samples,function(x)x@sample_name)
              return(res)
 })


setGeneric("extract.counts", function(x, ...) standardGeneric("extract.counts"))
setMethod("extract.counts",
          signature = "Iris",
          definition = function(x){
              counts <- lapply(x@samples,extract.counts.sample)
              nams <- sort(unique(unlist(lapply(counts,colnames))))
              for (i in 1:length(counts)){
                  counts[[i]] <- counts[[i]][,match(nams,colnames(counts[[i]]))]
                  colnames(counts[[i]]) <- nams
                  counts[[i]][is.na(counts[[i]])] <- 0
              }
              x@counts <- counts
              x@markers <- nams
              return(x)
})

setGeneric("extract.counts.sample", function(x, ...) standardGeneric("extract.counts.sample"))
setMethod("extract.counts.sample",
          signature = "Sample",
          definition = function(x){
              counts <- lapply(x@coordinates,function(x)table(x@ppp$marks))
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

setGeneric("get.counts.noncollapsed", function(x, ...) standardGeneric("get.counts.noncollapsed"))
setMethod("get.counts.noncollapsed",
          signature = "Iris",
          definition = function(x){
              counts <- x@counts
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

