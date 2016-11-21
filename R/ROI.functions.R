
setGeneric("extract.ROI.Coordinate", function(object,...) standardGeneric("extract.ROI.Coordinate"))
setMethod("extract.ROI.Coordinate",
          signature = "Coordinate",
          definition = function(object, ROI, normalize, all_levels){
              #reduce to the filter
              mask <- object@mask[[ROI]]
              filter <- sapply(1:length(object@ppp$x),function(i,dat,mask)mask[dat$x[i],dat$y[i]]==1,object@ppp,mask)
              object@ppp <- object@ppp[filter,]
              object@raw@data <- object@raw@data[filter,]
              
              #get the ROI counts
              ROI.counts <- table(object@ppp$marks)
              ROI.counts <- ROI.counts[all_levels]
              ROI.counts[is.na(ROI.counts)] <- 0
              
              #normalize the counts by the ROI size to make them comparable across images
              if (normalize){
                  size <- sum(mask)/length(mask)
                  ROI.counts <- ROI.counts/size
              }
              return(list(object=object,ROI.counts=ROI.counts))
          })

setGeneric("extract.ROI", function(object, ...) standardGeneric("extract.ROI"))
setMethod("extract.ROI",
          signature = "Iris",
          definition = function(object, ROI='filled_margin', normalize=T){
              all_levels <- sort(unique(unlist(lapply(object@counts,colnames))))
              ret <- lapply(object@samples, extract.ROI.sample, ROI, normalize, all_levels)
              object@samples <- lapply(ret,function(x)x$object)
              
              #update the counts
              object@counts <- lapply(ret,function(x)x$ROI.counts)
              
              #reset all spatial stats
              object@nearest_neighbors <- list()
              object@interactions <- list()
              object@proximity <- list()
              
              return(object)
          })



setGeneric("extract.ROI.sample", function(object, ...) standardGeneric("extract.ROI.sample"))
setMethod("extract.ROI.sample",
          signature = "Sample",
          definition = function(object, ROI, normalize, all_levels){
              ret <-lapply(object@coordinates, extract.ROI.Coordinate, ROI, normalize, all_levels)
              object@coordinates <- lapply(ret,function(x)x$object)
              ROI.counts <- t(sapply(ret,function(x)x$ROI.counts))
              colnames(ROI.counts) <- all_levels
              return(list(object=object, ROI.counts=ROI.counts))
          })