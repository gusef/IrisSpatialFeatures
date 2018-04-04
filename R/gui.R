#' Launch a GUI to look at the data
#'
#' @param x IrisSpatialFeatures ImageSet object
#'
#' @return launched shiny server
#'
#' @docType methods
#' @export
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' launch_shiny_gui(dataset)
#'
#' @rdname launch_shiny_gui
setGeneric("launch_shiny_gui", function(x)
    standardGeneric("launch_shiny_gui"))

#' @rdname launch_shiny_gui
#' @aliases ANY,ANY-method
setMethod(
    "launch_shiny_gui",
    signature = c("ImageSet"),
    definition = function(x) {
        print('haro thar')
        return(as.data.frame(x))
    }
)
