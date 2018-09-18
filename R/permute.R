#####################################################################################################################################
################ shuffle_labels
#####################################################################################################################################
#
### randomly shuffle the cell labels on each frame
##
##
### shuffle_labels(dataset,"SOX10+ PDL1-","SOX10+ PDL1+", "CD8+ PD1+",TRUE)
###
#' shuffle labels of each frame
#'
#' @param x IrisSpatialFeatures ImageSet object that has had extract nearest neighbors run
#'
#' @return ImageSet with labels shuffled
#'
#' @docType methods
#' @export
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' shuffle_labels(dataset)
#'
#' @rdname shuffle_labels
setGeneric("shuffle_labels", function(x, ...)
    standardGeneric("shuffle_labels"))

#' @rdname shuffle_labels
#' @aliases extract.nearest.neighbor,ANY,ANY-method
setMethod(
    "shuffle_labels",
    signature = c("ImageSet"),
    definition = function(x,subset=NULL) {
        x_copy <- x
        samples <- names(x_copy@samples)
        for (sample in samples) {
            frames <- names(x_copy@samples[[sample]]@coordinates)
            for (frame in frames) {
                marks <- x_copy@samples[[sample]]@coordinates[[frame]]@ppp$marks
                if (is.null(subset)) {
                    x_copy@samples[[sample]]@coordinates[[frame]]@ppp$marks <- sample(marks)
                } else {
                    to_place = x@samples[[sample]]@coordinates[[frame]]@ppp$marks %in% subset
                    to_shuffle = which(to_place)
                    #print(to_shuffle)
                    shuffled = sample(x@samples[[sample]]@coordinates[[frame]]@ppp$marks[to_shuffle])
                    #print(shuffled)
                    for (i in seq(1,length(shuffled))) {
                        marks[to_shuffle[i]] <- shuffled[i]
                    }
                    x_copy@samples[[sample]]@coordinates[[frame]]@ppp$marks <- marks
                }
            }
        }
        return(x_copy)
    }
)