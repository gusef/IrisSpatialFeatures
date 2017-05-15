#' Method that reduces the current dataset to a specific region of interest, discarding all cell coordinates outside of that region
#'  
#' @param x Iris ImageSet object 
#' @param ROI Region of interest (default: 'invasive_margin')
#' @param ... Additional arguments  
#' 
#' @return Iris ImageSet object
#' @examples
#' raw_data <- new("ImageSet")
#' raw_data <- read_raw(raw_data,
#'                      raw_dir_name=system.file("extdata", package = "Iris"),
#'                      format='Mantra')
#' dataset <- threshold_dataset(raw_data,
#'                              marker='PD-Ligand-1 (Opal 690)',
#'                              marker_name='PDL1',
#'                              base=c('SOX10+'))
#' dataset <- threshold_dataset(dataset,
#'                              marker='PD-1 (Opal 540)',
#'                              marker_name='PD1',
#'                              base=c('CD8+','OTHER'))                     
#' @docType methods
#' @export
#' @rdname extract_ROI
setGeneric("extract_ROI", function(x, ...) standardGeneric("extract_ROI"))

#' @rdname extract_ROI
#' @aliases extract_ROI,ANY,ANY-method
setMethod("extract_ROI",
          signature = "ImageSet",
          definition = function(x, ROI='invasive_margin'){
            if (length(x@samples[[1]]@coordinates[[1]]@mask[[ROI]])==0){
              stop('There is no mask for "',ROI,'"')
            }
            
            x@samples <- lapply(x@samples, extract_ROI_sample, ROI)
            
            #update the counts
            x <- extract_counts(x)
            
            #reset all spatial stats
            x@nearest_neighbors <- list()
            x@interactions <- list()
            x@proximity <- list()
            
            return(x)
          })

setGeneric("extract_ROI_sample", function(x, ...) standardGeneric("extract_ROI_sample"))
setMethod("extract_ROI_sample",
          signature = "Sample",
          definition = function(x, ROI){
            x@coordinates <- lapply(x@coordinates, extract_ROI_Coordinate, ROI)
            return(x)
          })


setGeneric("extract_ROI_Coordinate", function(x, ...) standardGeneric("extract_ROI_Coordinate"))
setMethod("extract_ROI_Coordinate",
          signature = "Coordinate",
          definition = function(x, ROI){
            #reduce to the filter
            mask <- x@mask[[ROI]]
            filter <- sapply(1:length(x@ppp$x),function(i,dat,mask)mask[dat$x[i],dat$y[i]]==1,x@ppp,mask)
            x@ppp <- x@ppp[filter,]
            x@raw@data <- x@raw@data[filter,]
            x@size_in_px <- sum(mask>0) 
            
            return(x)
          })
