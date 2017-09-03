#' Read inForm output and store it in an IrisSpatialFeatures ImageSet object.
#'
#' @param x Iris ImageSet boject
#' @param ... Additional arguments
#'
#' @return A dataframe
#' @examples
#' dataset <- IrisSpatialFeatures_data
#' dataframe <- as.data.frame(dataset)
#' @docType methods
#' @export
#' @rdname as.data.frame
#' @importFrom data.table rbindlist
setMethod("as.data.frame",
          signature = c(x="ImageSet"),
          function(x) {
              sample_names <- names(x@samples)
              sample_dfs <- lapply(sample_names,function(sample_name){
                  sample <- x@samples[[sample_name]]
                  frame_names <- names(sample@coordinates)
                  frame_dfs <- lapply(frame_names,function(frame_name){
                      frame <- sample@coordinates[[frame_name]]
                      coords <- frame@ppp
                      frame_df <- as.data.frame(coords)
                      frame_df$frame <- rep(frame_name,dim(frame_df)[1])
                      return(frame_df)
                  })
                  sample_df <- rbindlist(frame_dfs)
                  sample_df$sample <- rep(sample_name,dim(sample_df)[1])
                  return(sample_df)
              })
              imageset_df <- rbindlist(sample_dfs)
              return(imageset_df)
})
