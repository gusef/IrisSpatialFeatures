##################################
####### Get all count data


setGeneric("get_counts_collapsed",
            function(x, ...) standardGeneric("get_counts_collapsed"))
setMethod("get_counts_collapsed",
          signature = "ImageSet",
          definition = function(x){
              combined <- sapply(x@counts,colSums,na.rm=TRUE)
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
#' 
#' @param x Iris ImageSet object.
#' @param ... Additional arguments  
#' @examples
#' get_counts_per_mm2_noncollapsed(new("ImageSet"))
#' 
#' @return Iris ImageSet object.
#' @docType methods
#' @export
#' @rdname get_counts_per_mm2_noncollapsed
setGeneric("get_counts_per_mm2_noncollapsed", 
            function(x, ...) standardGeneric("get_counts_per_mm2_noncollapsed"))

#' @rdname get_counts_per_mm2_noncollapsed
#' @aliases get_counts_per_mm2_noncollapsed,ANY,ANY-method
setMethod("get_counts_per_mm2_noncollapsed",
          signature = "ImageSet",
          definition = function(x){
              
              sizes <- lapply(x@samples,
                              function(y)sapply(y@coordinates,
                                                function(z)z@size_in_px))
              
              #counts per mm2
              counts <- lapply(x@counts,
                               function(y,z) y/(z@microns_per_pixel^2),
                               x)
              samps <- names(sizes)
              counts <- lapply(samps,
                               function(y,counts,sizes)
                                  1000000*sweep(counts[[y]],1,sizes[[y]],'/'),
                               counts,
                               sizes)
              names(counts) <- samps
              
              return(counts)
          })


#' Get all the counts on a per mm2 basis
#' @param x An Iris ImageSet object
#' @param digits Number of digits that are shown in the output (default: 2)
#' @param blank (default: FALSE)
#' @param ... Additional arguments
#' @return counts per mm2 per sample, collapsing each coordinate and returning mean and standard error
#' 
#' @examples
#' get_counts_per_mm2(new("ImageSet"))
#' 
#' @docType methods
#' @export
#' @importFrom stats sd
#' @rdname get_counts_per_mm2
setGeneric("get_counts_per_mm2", function(x, ...) standardGeneric("get_counts_per_mm2"))

#' @rdname get_counts_per_mm2
setMethod("get_counts_per_mm2",
          signature = "ImageSet",
          definition = function(x, digits=2, blank=FALSE){
              counts <- get_counts_per_mm2_noncollapsed(x)
              if (length(x@counts)>1){
                  means <- sapply(counts,colMeans,na.rm=TRUE)
                  se <- sapply(counts,function(x)apply(x,2,
                     function(y)sd(y,na.rm = TRUE)/sqrt(length(y[!is.na(y)]))))
                  res <- means
                  if (!blank){
                      for (i in 1:ncol(means)){
                          res [,i] <- paste(format(means[,i],digits=digits),
                                            '+/-',
                                            format(se[,i],digits=digits))   
                      }
                  }
              }else if (!blank){
                  res <- format(counts,digits=digits)
              }else{
                  res <- counts
              }
              return(res)
          })


#' Get ratio of counts between two markers
#' 
#' @param x An Iris object
#' @param marker1 First cell-type.
#' @param marker2 Second cell-type.
#' @param digits Number of digits that should be shown in the the results. (Default: 2)
#' @param ... Additional arguments.
#' 
#' @docType methods
#' @importFrom stats sd
#' @rdname get_count_ratios
#' @export
#' @examples 
#' raw_data <- new("ImageSet")
#' raw_data <- read_raw(raw_data,
#'                      raw_dir_name=system.file("extdata", package = "Iris"),
#'                      format='Mantra')
#' dataset <- threshold_dataset(raw_data,
#'     marker='PD-Ligand-1 (Opal 690)',
#'     marker_name='PDL1',
#'     base=c('SOX10+'))
#' dataset <- threshold_dataset(dataset,
#'     marker='PD-1 (Opal 540)',
#'     marker_name='PD1',
#'     base=c('CD8+','OTHER'))
#'     get_count_ratios(dataset,'SOX10+ PDL1-','SOX10+ PDL1+')
setGeneric("get_count_ratios", 
           function(x, ...) standardGeneric("get_count_ratios"))

#' @rdname get_count_ratios
#' @aliases get_count_ratios,ANY,ANY-method
setMethod("get_count_ratios",
          signature = "ImageSet",
          definition = function(x, marker1, marker2, digits=2){
              ratios <- sapply(x@counts,
                               function(x,m1,m2)x[,m1]/x[,m2],marker1,marker2)
              for (idx in 1:length(ratios)){
                 ratios[[idx]][is.infinite(ratios[[idx]])] <- NA
              }
              means <- sapply(ratios,mean,na.rm=TRUE)
              se <- sapply(ratios,
                 function(x)sd(x,na.rm = TRUE)/sqrt(length(x[!is.na(x)])))
              res <- paste(format(means,digits=digits),
                           '+/-',
                           format(se,digits=digits))   
              names(res) <- sapply(x@samples,function(x)x@sample_name)
              return(res)
 })

setGeneric("extract_counts", function(x, ...) standardGeneric("extract_counts"))
setMethod("extract_counts",
          signature = "ImageSet",
          definition = function(x){
              counts <- lapply(x@samples,extract_counts_sample)
              nams <- sort(unique(unlist(lapply(counts,colnames))))
              for (i in 1:length(counts)){
                  if (nrow(counts[[i]])==1){
                      temp <- t(as.matrix(
                         counts[[i]][,match(nams,colnames(counts[[i]]))]))
                      rownames(temp) <- rownames(counts[[i]])
                      counts[[i]] <- temp
                  }else{
                      counts[[i]] <- counts[[i]][,match(nams,
                                                        colnames(counts[[i]]))]
                  }
                  
                  colnames(counts[[i]]) <- nams
                  counts[[i]][is.na(counts[[i]])] <- 0
              }
              x@counts <- counts
              x@markers <- nams
              return(x)
})

setGeneric("extract_counts_sample", 
           function(x, ...) standardGeneric("extract_counts_sample"))
setMethod("extract_counts_sample",
          signature = "Sample",
          definition = function(x){
              counts <- lapply(x@coordinates,function(x)table(x@ppp$marks))
              nams <- unique(unlist(lapply(counts,names)))
              counter <- rep(0,length(nams))
              names(counter) <- nams
              counts <- t(sapply(counts,extractCountsF,counter))
              return(counts)
})

#helperfunction to count the features making sure missing celltypes don't 
#  cause problems
extractCountsF <- function(x,counter){
    counter[match(names(x),names(counter))] <- x
    return(counter)
}

setGeneric("get_counts_noncollapsed", 
           function(x, ...) standardGeneric("get_counts_noncollapsed"))
setMethod("get_counts_noncollapsed",
          signature = "ImageSet",
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

