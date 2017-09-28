#' Extract the distance to each nearest neighbor for each cell-type
#'
#' @param x IrisSpatialFeatures ImageSet object
#' @param min_num_cells Minimum number of cell that a coordinate needs to have in order to calculate the statistics (Default: 10)
#' @param ... Additional arguments
#'
#' @return distance to nearest neighbor for each
#'
#' @docType methods
#' @export
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' extract_nearest_neighbor(dataset)
#'
#' @rdname extract_nearest_neighbor
setGeneric("extract_nearest_neighbor", function(x, ...)
    standardGeneric("extract_nearest_neighbor"))

#' @rdname extract_nearest_neighbor
#' @aliases extract.nearest.neighbor,ANY,ANY-method
setMethod(
    "extract_nearest_neighbor",
    signature = "ImageSet",
    definition = function(x, min_num_cells = 10) {
        all_levels <- x@markers
        x@nearest_neighbors <- lapply(x@samples,
                                      nearest_neighbor_sample,
                                      all_levels,
                                      min_num_cells)
        return(x)
    }
)

#' Compare nearest neighbors by a data.frame
#'
#' @param x IrisSpatialFeatures ImageSet object that has had extract nearest neighbors run
#' @param markerA First marker
#' @param markerB Second marker
#' @param reference Reference marker
#' @param from_reference If true calculate distance from the reference to the markers by NN
#'
#' @return data.frame of markers and distances
#'
#' @docType methods
#' @export
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' extract_nearest_neighbor(dataset)
#'
#' @importFrom data.table rbindlist
#' @rdname nn_comparison_dataframe
setGeneric("nn_comparison_dataframe", function(x,markerA,markerB,reference,from_reference=TRUE)
    standardGeneric("nn_comparison_dataframe"))

#' @rdname nn_comparison_dataframe
#' @aliases extract.nearest.neighbor,ANY,ANY-method
setMethod(
    "nn_comparison_dataframe",
    signature = c("ImageSet","character","character","character","logical"),
    definition = function(x, markerA, markerB, reference, from_reference=TRUE) {
        samples <- names(nn@nearest_neighbors)
        neighbors <- lapply(samples,function(sample) {
            means = nn@nearest_neighbors[sample][[1]]$means
            if (from_reference) {
                v1 = data.frame(sample=sample,markerA=markerA,markerB=markerB,reference=reference,from_reference=from_reference,distanceA=means[reference,markerA],distanceB=means[reference,markerB])
                return(v1)
            } else {
                v1 = data.frame(sample=sample,markerA=markerA,markerB=markerB,reference=reference,from_reference=from_reference,distanceA=means[markerA,reference],distanceB=means[markerB,reference])
                return(v1)
            }
        })
        return(rbindlist(neighbors))

    }
)

setGeneric("nearest_neighbor_sample", function(x, ...)
    standardGeneric("nearest_neighbor_sample"))
setMethod(
    "nearest_neighbor_sample",
    signature = "Sample",
    definition = function(x, all_levels, min_num_cells) {
        res <-
            lapply(x@coordinates,
                   nearest_neighbor_coord_raw,
                   all_levels,
                   min_num_cells)

        means <- lapply(res, function(x)
            x$means)
        vars <- lapply(res, function(x)
            x$vars)
        nums <- lapply(res, function(x)
            x$nums)

        #collapse the different coordinates
        means <- collapseMatrices(means, rowMeans)
        vars <- collapseMatrices(vars, rowMeans)
        nums <- collapseMatrices(nums, rowSums)

        #there is a special case where there is only 1 cell of a type
        #this leads to an NA variance. In this case the means are set to NA as well
        means[is.na(vars)] <- NA

        #calculate the standard error on the combined coordinates
        ses <- sqrt(vars) / sqrt(nums)
        return(list(means = means, SE = ses))
    }
)

setGeneric("nearest_neighbor_coord_raw", function(x, ...)
    standardGeneric("nearest_neighbor_coord_raw"))
setMethod(
    "nearest_neighbor_coord_raw",
    signature = "Coordinate",
    definition = function(x, all_levels, min_num_cells) {
        ppp <- x@ppp
        res <-
            lapply(all_levels, getToNeighbors, all_levels, ppp, min_num_cells)
        means <- t(sapply(res, function(x)
            x['means', ]))
        vars <- t(sapply(res, function(x)
            x['vars', ]))
        nums <- t(sapply(res, function(x)
            x['nums', ]))

        rownames(nums) <-
            rownames(means) <- rownames(vars) <- colnames(means)

        #in some cases there are not all samples represented
        means[!is.finite(means)] <- NA
        vars[!is.finite(vars)] <- NA
        nums[!is.finite(nums)] <- NA

        return(list(
            means = means,
            vars = vars,
            nums = nums
        ))
    }
)


#' @importFrom spatstat nncross
#' @importFrom stats var
getNearestNeighbor <- function(from, to, ppp, min_num_cells) {
    if (sum(ppp$marks == from) < min_num_cells |
        sum(ppp$marks == to) < min_num_cells) {
        dis <- NA
    } else{
        dis <- nncross(ppp[ppp$marks == from, ], ppp[ppp$marks == to, ])[, 1]
    }
    means <- mean(dis)
    vars <- var(dis)
    nums <- sum(ppp$marks == from)
    return(c(
        means = means,
        vars = vars,
        nums = nums
    ))
}

getToNeighbors <- function(to, classes, ppp, min_num_cells) {
    res <- sapply(classes, getNearestNeighbor, to, ppp, min_num_cells)
    return(res)
}

################################################################
##### Nearest Neighbors getters


#' Get the nearest neighbor for each cell-type
#'
#' @param x An IrisSpatialFeatures ImageSet object
#' @param ... Additional arguments
#' @return Nearest neighbor for each cell-type
#'
#' @docType methods
#' @export
#'
#' @examples
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_nearest_neighbor(dataset)
#' get_all_nearest_neighbors(dataset)
#'
#' @rdname get_all_nearest_neighbors
setGeneric("get_all_nearest_neighbors", function(x, ...)
    standardGeneric("get_all_nearest_neighbors"))

#' @rdname get_all_nearest_neighbors
#' @aliases get_all_nearest_neighbors,ANY,ANY-method
setMethod(
    "get_all_nearest_neighbors",
    signature = "ImageSet",
    definition = function(x) {
        return(x@nearest_neighbors)
    }
)

#' Get the nearest neighbor for a specified cell-type
#'
#' @param x An IrisSpatialFeatures ImageSet object
#' @param marker Cell type for which the nearest neighbor should be calculated
#' @param ... Additional arguments
#' @return nearest neighbors for the specified cell-type
#'
#' @export
#' @rdname get_nearest_neighbors
#' @examples
#' #' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_nearest_neighbor(dataset,min_num_cells=2)
#' get_nearest_neighbors(dataset,"SOX10+ PDL1+")
#'
setGeneric("get_nearest_neighbors", function(x, ...)
    standardGeneric("get_nearest_neighbors"))

#' @rdname get_nearest_neighbors
#' @aliases get_nearest_neighbors,ANY,ANY-method
setMethod(
    "get_nearest_neighbors",
    signature = "ImageSet",
    definition = function(x, marker) {
        if (!marker %in% x@markers) {
            stop(paste('There is no celltype: ', marker))
        }

        nn <-
            sapply(x@nearest_neighbors, function(x, y)
                x$means[, y], marker)
        se <-
            sapply(x@nearest_neighbors, function(x, y)
                x$SE[, y], marker)
        return(list(mean = nn, SE = se))
    }
)

###############################################################################################################
##################################### Nearest neighbor plotting functions
###############################################################################################################

#' Plot average nearest neighbor barplots for two cell types. This measurement is not symmetric, so if 'from' and 'to' are switched it will result in different results.
#' For the 'to' parameter this function allows a cell-type without '+' or '-' in the end. Indicating that the distances from the first cell-type should be calculated
#' against both '+/-' and a paired t-test should be calculated. For example we want to calculate the average distance between SOX10 PDL1+ melanoma cells against
#' both CD8 PD1+ and CD8 PD1- cells, the 'CD8 PD1' would be speficified as 'to' parameter, 2 distances would be calculated for each sample and a two-sided paired t-test calculated
#' to test for significant differences.
#'
#' @param x IrisSpatialFeatures ImageSet object.
#' @param from Cell-type from which the nearest neighbor is calculated.
#' @param to Cell-type to which the nearest neighbor is calculated.
#' @param ttest Flag indicating whether a paired t-test should be calculated. (default: TRUE)
#' @param transposed Switches 'from' and 'to' cell-type. This way the (default: FALSE)
#' @param remove_NAs dont plot samples with less than min cells
#' @param ... Additional arguments.
#' @return plot average nearest neighbor barplots for two cell types
#'
#' @importFrom graphics barplot
#' @importFrom graphics mtext
#' @importFrom stats t.test
#' @docType methods
#' @export
#' @rdname plot_nearest_neighbor
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_nearest_neighbor(dataset)
#' p <- plot_nearest_neighbor(dataset,'CD8+ PD1+','SOX10+ PDL1')
#'
setGeneric("plot_nearest_neighbor", function(x, ...)
    standardGeneric("plot_nearest_neighbor"))

#' @rdname plot_nearest_neighbor
#' @aliases plot_nearest_neighbor,ANY,ANY-method
setMethod(
    "plot_nearest_neighbor",
    signature = "ImageSet",
    definition = function(x,
                          from,
                          to,
                          ttest = TRUE,
                          transposed = FALSE,
                          remove_NAs = FALSE) {
        marker_names <- x@markers

        #grab the relevant markers
        comp <- grep(to, marker_names, fixed = TRUE)

        if (length(comp) == 0 || !from %in% marker_names) {
            stop('One or both selected markers are not included in the dataset')
        }

        #drop distances to self
        comp <- comp[!marker_names[comp] %in% from]

        x.mean <-
            extractNNVals(x@nearest_neighbors, 'means', from, comp[1], transposed)
        x.se <-
            extractNNVals(x@nearest_neighbors, 'SE', from, comp[1], transposed)

        if (length(comp) > 1) {
            y.mean <-
                extractNNVals(x@nearest_neighbors, 'means', from, comp[2], transposed)
            y.se <-
                extractNNVals(x@nearest_neighbors, 'SE', from, comp[2], transposed)
            current.mean <- (t(cbind(x.mean, y.mean)))
            #print(current.mean)
            current.se <- (t(cbind(x.se, y.se)))
        } else{
            current.mean <- t(x.mean)
            current.se <- t(x.se)
        }

        #remove NA values
        if (remove_NAs){
            drop <- colSums(is.na(current.mean)) > 0
            current.mean <- current.mean[,!drop]
            current.se <- current.se[,!drop]
        }else{
            current.se[is.na(current.mean)] <- 0
            current.mean[is.na(current.mean)] <- 0
        }

        max_idx <- which.max(rowSums(current.mean))
        ord <- order(current.mean[max_idx, ], decreasing = TRUE)
        current.mean <- current.mean[, ord]
        current.se <- current.se[, ord]

        #Fixing the legend for the plot
        if (length(comp) == 1) {
            ext <- ''
            leg <- marker_names[comp]
            COLS <- c("lightgrey")
        } else{
            ext <- '+/-'
            leg <- c(marker_names[comp[1]],
                     marker_names[comp[2]])
            COLS <- c("lightgrey", "black")
        }
        label <- buildLabel(from, to, ext, transposed)
        bp <- barplot(
            current.mean,
            main = label,
            xlab = "",
            ylab = "Avg. distance to NN",
            col = COLS,
            legend = leg,
            ylim = c(0, max(current.mean) + max(current.se)),
            las = 2,
            beside = TRUE
        )

        if (length(comp) > 1) {
            # Add a check for NA values if there is an NA value in one row, we need it in the other row.
            #print(current.mean)
            plotSE(bp, current.mean, current.se, 1)
            plotSE(bp, current.mean, current.se, 2)
        } else{
            bp <- t(bp)
            current.mean <- t(current.mean)
            current.se <- t(current.se)
            plotSE(bp, current.mean, current.se, 1)
        }
        #paired t test to test for significance
        if (ttest & length(comp) > 1) {
            current.mean <- current.mean[,colSums(current.mean!=0)==2]
            current.se <- current.se[,colSums(current.se!=0)==2]
            pval <-
                t.test(current.mean[1, ], current.mean[2, ], paired = TRUE)$p.value
            mtext(paste('Paired t-test:', format(pval, digits = 4)), 3)
        } else{
            pval <- NA
        }
        return(list(
            means = current.mean,
            ses = current.se,
            pval = pval,
            label = label
        ))
    }
)

#' @importFrom graphics arrows
plotSE <- function(bp, current.mean, current.se, idx) {
    arrows(
        bp[idx, ],
        current.mean[idx, ] + current.se[idx, ],
        bp[idx, ],
        current.mean[idx, ],
        angle = 90,
        code = 3,
        length = 0.02,
        col = 'black'
    )
}

#extract the marker pos/neg values across the samples
extractNNVals <- function(mat, val, from, to, transposed) {
    if (transposed) {
        dest <- from
        src <- to
    } else{
        dest <- to
        src <- from
    }
    return(sapply(mat, function(x, y, z)
        x[[val]][z, y], src, dest))
}

#generate the label of the plot
buildLabel <- function(from, to, ext, transposed) {
    if (transposed) {
        label <- paste('Distance from', to, ext, 'to', from)
    } else{
        label <- paste('Distance from', from, 'to', to, ext)
    }
    return(label)
}

###############################################################################################################
######################################ray plots
###############################################################################################################

#' Plot nearest neighbor ray plots for each samples
#'
#' @param x An IrisSpatialFeatures ImageSet object
#' @param from_type Cell type from which the rays are drawn
#' @param to_type Cell type to which the rays are drawn
#' @param from_col Color for the 'from' cell-type (Default: '#EE7600')
#' @param to_col Color for the 'to' cell-type  (Default: '#028482')
#' @param format Format of the output file, can be '.pdf' or '.png' (Default: '.pdf')
#' @param plot_dir Directory in which the images are written (Default: './')
#' @param lineColor Color of the rays (Default: '#666666')
#' @param height Height of the pdf. (Default: 7)
#' @param width Width of the pdf. (Default: 10)
#' @param ... Additional arguments.
#' @return nearest neighbor ray plots
#'
#' @docType methods
#' @export
#'
#' @importFrom spatstat superimpose
#' @rdname neighbor_ray_plot
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_nearest_neighbor(dataset)
#' get_nearest_neighbors(dataset,"SOX10+ PDL1+")
#' plot_dir <- file.path('./ray_plots')
#' if (!file.exists(plot_dir)){
#'  dir.create(file.path(plot_dir))
#' }
setGeneric("neighbor_ray_plot", function(x, ...)
    standardGeneric("neighbor_ray_plot"))

#' @rdname neighbor_ray_plot
#' @aliases neighbor.ray.plot,ANY,ANY-method
setMethod(
    "neighbor_ray_plot",
    signature = "ImageSet",
    definition = function(x,
                          from_type,
                          to_type,
                          from_col = '#EE7600',
                          to_col = '#028482',
                          format = '.pdf',
                          plot_dir = './',
                          lineColor = '#666666',
                          height = 7,
                          width = 10) {
        #generate the mapping directory
        #out_dir <- file.path(getwd(), plot_dir)
        out_dir <- plot_dir
        if (!file.exists(out_dir)) {
            dir.create(out_dir, showWarnings = FALSE)
        }

        #generate ray plots for each sample
        lapply(
            x@samples,
            neighbor_ray_plot_sample,
            from_type,
            to_type,
            from_col,
            to_col,
            out_dir,
            format,
            lineColor,
            height,
            width
        )
    }
)

setGeneric("neighbor_ray_plot_sample", function(x, ...)
    standardGeneric("neighbor_ray_plot_sample"))
setMethod(
    "neighbor_ray_plot_sample",
    signature = "Sample",
    definition = function(x,
                          from_type,
                          to_type,
                          from_col,
                          to_col,
                          out_dir,
                          format,
                          lineColor,
                          height,
                          width) {
        lapply(
            x@coordinates,
            neighbor_ray_plot_coord,
            x@sample_name,
            from_type,
            to_type,
            from_col,
            to_col,
            out_dir,
            format,
            lineColor,
            height,
            width
        )
    }
)


#' @importFrom spatstat nncross
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom grDevices png
setGeneric("neighbor_ray_plot_coord", function(x, ...)
    standardGeneric("neighbor_ray_plot_coord"))
setMethod(
    "neighbor_ray_plot_coord",
    signature = "Coordinate",
    definition = function(x,
                          samp_name,
                          from_type,
                          to_type,
                          from_col,
                          to_col,
                          out_dir,
                          format,
                          lineColor,
                          height,
                          width) {


        #extract the relevant cells
        if (sum(x@ppp$marks == from_type) > 0 &&
            sum(x@ppp$marks == to_type) > 0) {

            file_stub <- paste0(samp_name, '_', x@coordinate_name)
            if (format == '.pdf') {
                pdf(
                    file = file.path(out_dir, paste0(file_stub, '.pdf')),
                    width = width,
                    height = height
                )
            } else if (format == '.png') {
                png(file.path(out_dir, paste0(file_stub, '.png')),
                    width = 800,
                    height = 600)
            }

            rayplot_single_coordinate(x,
                                      from,
                                      to,
                                      samp_name,
                                      from_col = '#EE7600',
                                      to_col = '#028482',
                                      lineColor = '#666666')

            dev.off()
        }
    }
)


#' Plot nearest neighbor ray plots for a single coordinate
#'
#' @param x An IrisSpatialFeatures ImageSet object
#' @param from_col Cell type from which the rays are drawn
#' @param to_type Cell type to which the rays are drawn
#' @param samp_name Name of the sample
#' @param from_col Color for the 'from' cell-type (Default: '#EE7600')
#' @param to_col Color for the 'to' cell-type  (Default: '#028482')
#'
#' @export
#'
#' @importFrom spatstat nncross
#' @importFrom spatstat superimpose
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @rdname rayplot_single_coordinate
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_nearest_neighbor(dataset)
#' rayplot_single_coordinate(x = dataset@samples[[1]]@coordinates[[1]],
#'                           samp_name = dataset@samples[[1]]@sample_name,
#'                           from_type = "SOX10+ PDL1+",
#'                           to_type = "CD8+ PD1+")

rayplot_single_coordinate <- function(x,
                                      from_type,
                                      to_type,
                                      samp_name = '',
                                      from_col = '#EE7600',
                                      to_col = '#028482',
                                      lineColor = '#666666'){

    from <- x@ppp[x@ppp$marks == from_type, ]
    to <- x@ppp[x@ppp$marks == to_type, ]
    #get limits
    overlap <- superimpose(from, to)
    xlim <- max(data.frame(overlap)$x)
    ylim <- max(data.frame(overlap)$y)

    #get distances
    dist <- nncross(from, to)
    nearest <- data.frame(cbind(data.frame(from),
                                dist$dist,
                                data.frame(to[dist$which, ])))
    names(nearest) <-
        c('from_x',
          'from_y',
          'from_type',
          'dist',
          'to_x',
          'to_y',
          'to_type')

    #get the colors right
    df <- data.frame(overlap)
    df$cols <- as.character(df$marks)
    df$cols[df$cols == from_type] <- from_col
    df$cols[df$cols == to_type] <- to_col

    par(mar = c(4, 4, 4, 1))
    plot(
        df$x,
        df$y,
        col = as.character(df$cols),
        pch = 18,
        ylab = 'y (pixels)',
        xlab = 'x (pixels)',
        main = paste(samp_name, '-', x@coordinate_name)
    )
    segments(nearest$from_x,
             nearest$from_y,
             nearest$to_x,
             nearest$to_y,
             col = lineColor)
    legend(
        'bottomleft',
        col = c(from_col, to_col),
        legend = c(from_type, to_type),
        pch = 18,
        cex = 0.8
    )
}


#####################################################################################################################################
################ Normalized nearest neighbor functions
#####################################################################################################################################




#' Extract the distance to each nearest neighbor for specified
#' cell-types, normalized by downsampling each cell-type to the
#' same size, for a single sample, with no resampling
#'
#' @param sample_name sample_name sample name string
#' @param data IrisSpatialFeatures ImageSet object
#' @param markers vector of marker names to use
#' @param minimum_cells the smallest number of cells (default:50)
#' @param grouped_sample TRUE/FALSE if we want to group samples together and
#'                       thus normalize the frames to the smallest frame
#'                       count (Default: TRUE)
#'
#' @return data.frame
#'
#' @importFrom spatstat nncross
setGeneric("normal_nearest_neighbor_sample_once", function(sample_name,
                                                           data,
                                                           markers,
                                                           minimum_cells=50,
                                                           grouped_sample=TRUE)
    standardGeneric("normal_nearest_neighbor_sample_once"))

#' @rdname normal_nearest_neighbor_sample_once
#' @aliases normal.nearest.neighbor.sample.once,ANY,ANY-method
setMethod(
    "normal_nearest_neighbor_sample_once",
    signature(sample_name="character",data="ImageSet"),
    definition <- function(sample_name,data,markers,minimum_cells,grouped_sample) {
    # For a single sample designated by sample_name get a dataframe
    contains_markers <- data@markers[data@markers %in% markers]
    if(length(contains_markers)!=length(markers)) {
        stop("marker name problem")
    }
    sample <- data@samples[sample_name][[1]]
    frame_names <- names(sample@coordinates)
    # First lets get the smallest cell count
    functional_frame_names <- lapply(frame_names,function(frame_name){
        #get the smallest cell counts from the frames that have enough cells
        dat <- sample@coordinates[frame_name][[1]]
        mcnt <- min(sapply(markers,function(x){sum(dat@ppp$marks==x)}))
        if (mcnt >= minimum_cells) { return(TRUE)}
        return(FALSE)
    })
    true_minimum <- minimum_cells
    min_counts <- sapply(frame_names,function(frame_name){
        #get the smallest cell counts from the frames that have enough cells
        dat <- sample@coordinates[frame_name][[1]]
        mcnt <- min(sapply(markers,function(x){sum(dat@ppp$marks==x)}))
        if (mcnt < minimum_cells) { return(NA)}
        return(mcnt)
    })
    if(!all(is.na(min_counts))) {
        #if there is real number in there
        true_minimum <- min(min_counts,na.rm=TRUE)
        #print(true_minimum)
    }
    names(functional_frame_names) <- frame_names
    #print(functional_frame_names)
    #print(smallest_cell_count)
    smallest_cell_count <- true_minimum
    frame_df_list <- lapply(frame_names,function(frame_name){
        dat <- sample@coordinates[frame_name][[1]]
        # filter down to just the markers we're interested in
        tot <- sapply(markers,function(x){sum(dat@ppp$marks==x)})
        # get the number of cells in each of the categories of interest
        if (grouped_sample==FALSE) {
            smallest_cell_count <- min(tot)
        }
        # get the number to downsample to
        parr <- lapply(markers,function(x){
            mppp<-dat@ppp[dat@ppp$marks==x,]
            if (grouped_sample==TRUE) {
                if(functional_frame_names[frame_name][[1]]==TRUE) {
                    # if its grouped and it is going to get used, do it right
                    mppp<-mppp[sample(1:length(mppp$marks),true_minimum),]
                    return(mppp)
                } else {
                    return(mppp)
                }
            } else {
                # otherwise use its own cell count
                mppp<-mppp[sample(1:length(mppp$marks),smallest_cell_count),]
                return(mppp)
            }
        })
        names(parr) <- markers
        # exectue downsampling
        nn_df_list <- lapply(markers,function(marki){
            # Get the mean aand variance between all markers a list of lists
            pi <- parr[marki][[1]]
            outs <- lapply(markers,function(markj){
                # Get the mean and variance for nnearest distances bettween markj and marki
                pj <- parr[markj][[1]]
                #print(functional_frame_names[frame_name][[1]])
                if(smallest_cell_count < minimum_cells || functional_frame_names[frame_name][[1]]==FALSE) {
                    #print("return bad")
                    return(list(mean_dist=NA,
                                var_dist=NA))
                }
                #print(pj)
                #print(pi)
                #print('----')
                dis<-spatstat::nncross(pi,pj)[,1]
                res <- list(mean_dist=mean(dis),
                            var_dist=var(dis))
                return(res)
            })
            names(outs) <- markers
            mean_dist_arr <- sapply(markers,function(x){outs[[x]]$mean_dist})
            names(mean_dist_arr) <- markers
            var_dist_arr <- sapply(markers,function(x){outs[[x]]$var_dist})
            names(var_dist_arr) <- markers
            # Begin forming our dataframes early
            df <- cbind(as.data.frame(rep(marki,length(mean_dist_arr))),
                        as.data.frame(markers),
                        as.data.frame(mean_dist_arr),
                        as.data.frame(var_dist_arr))
            rownames(df) <- NULL
            colnames(df) <- c('marker_i','marker_j','mean','var')
            return(df)
        })
        nn_df <- do.call("rbind",nn_df_list)
        # concatonate the data frames
        nn_df$smallest_cell_count <- rep(smallest_cell_count,dim(nn_df)[1])
        # include our smallest cell count
        return(nn_df)
    })
    names(frame_df_list) <- frame_names
    if (grouped_sample==TRUE) {
        ### Case 1: We are grouping sample frames together
        #now we can get the mean and variance matrix from the mean of the frames
        # Remove frames that did not have enough data for the calculation

        #notna <- sapply(frame_names,function(x){
        #    if(is.na(frame_df_list[x][[1]]$mean[1])) { return(FALSE);}
        #    return(TRUE)
        #})
        #new_frame_names <- frame_names
        #if (length(notna[notna==TRUE])>0) { new_frame_names <- frame_names[notna];}
        usednames <- frame_names[unlist(functional_frame_names)]
        mean_data <- lapply(usednames,function(x){
            return(frame_df_list[x][[1]]$mean)
        })
        var_data <- lapply(usednames,function(x){
            return(frame_df_list[x][[1]]$var)
        })
        populations <- lapply(usednames,function(x){
            return(frame_df_list[x][[1]]$smallest_cell_count)
        })
        #Combine the frames to get aggrogate statistics of all the frames
        mean_combined = NA
        if (length(mean_data)>0) {
            mean_combined <- Reduce("+",mean_data)/length(mean_data)
        }
        #print(mean_combined)
        var_combined = NA
        if (length(var_data)>0) {
            var_combined <- Reduce("+",var_data)/length(var_data)
        }
        #print(var_combined)
        min_pop <- Reduce("min",populations)
        #print(min_pop)
        if(is.null(min_pop)) {min_pop=NA}
        else{ min_pop=min_pop[1]}
        #print(min_pop)
        #print(frame_df_list)
        #max_pop <- Reduce("max",populations)
        # Build a data frame with our data
        template <- frame_df_list[frame_names[1]][[1]]
        #print(template)
        df <- data.frame(marker_i=template$marker_i,
                marker_j=template$marker_j,
                mean=mean_combined,
                var=var_combined,
                original_frame_count=rep(length(frame_names),dim(template)[1]),
                useful_frame_count=rep(length(usednames),dim(template)[1]),
                #min_frame_cells=rep(min_pop,dim(template)[1]),
                #max_frame_cells=rep(max_pop,dim(template)[1]),
                smallest_cell_count=rep(min_pop,dim(template)[1]),
                sample=rep(sample_name,dim(template)[1])
                )
        return(df)
    } else if (grouped_sample==FALSE) {
        ### Case 2: We are leaving frames separate
        named_frames <-lapply(frame_names,function(x){
            # name the dataframes
            framedf <- frame_df_list[x][[1]]
            framedf$frame <- rep(x,dim(framedf)[1])
            #framedf$frame_cells <- rep(framedf$smallest_cell_count
            return(framedf)
        })
        nf_df <- do.call("rbind",named_frames)
        nf_df <- data.frame(nf_df,
                            sample=rep(sample_name,dim(nf_df)[1]),
                            check.names = FALSE)
        return(nf_df)
    }
})

#' Extract the distance to each nearest neighbor for specified
#' cell-types, normalized by downsampling each cell-type to the
#' same size (the smallest population from among the specified
#' markers), calculates for a single specified sample
#'
#' @param sample_name string name of the sample
#' @param data IrisSpatialFeatures ImageSet object
#' @param markers vector of marker names to use
#' @param n_resamples number of times to resample each frame (default:500)
#' @param minimum_cells smallest number of cells to consider a frame (default:50)
#' @param quantiles vector of numeric fractions to include in vector
#'        to show the mean distance calculated across resamplings
#'        (default:c(0.05,0.25,0.5,0.75,0.95))
#' @param grouped_sample TRUE/FALSE group samples together (default:TRUE)
#'
#' @return data.frame
#'
#' @docType methods
#' @export
#'
#' @examples

#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' normal_nearest_neighbor_sample("MEL2",dataset,c("SOX10+ PDL1+","SOX10+ PDL1-"),10)
#'
#' @rdname normal_nearest_neighbor_sample
#' @importFrom spatstat nncross
#' @importFrom matrixStats rowMedians
#' @importFrom matrixStats rowQuantiles
setGeneric("normal_nearest_neighbor_sample", function(sample_name,
                                                      data,markers,
                                                      n_resamples=500,
                                                      minimum_cells=50,
                                                      quantiles=c(0.05,0.25,0.5,0.75,0.95),
                                                      grouped_sample=TRUE)
    standardGeneric("normal_nearest_neighbor_sample"))

#' @rdname normal_nearest_neighbor_sample
#' @aliases normal.nearest.neighbor.sample,ANY,ANY-method
setMethod(
    "normal_nearest_neighbor_sample",
    signature(sample_name="character",data="ImageSet"),
    definition <- function(sample_name,data,markers,n_resamples,minimum_cells,quantiles,grouped_sample) {
    totals<-lapply(rep(sample_name,n_resamples),
                   normal_nearest_neighbor_sample_once,
                   data=data,
                   markers=markers,
                   minimum_cells=minimum_cells,
                   grouped_sample = grouped_sample)
    combine_mean <- sapply(totals,function(x){x$mean})
    combine_var <- sapply(totals,function(x){x$var})
    template <- totals[[1]]
    #build the dataframe
    if (grouped_sample==TRUE) {
        ### Case 1: We are putting samples frames together
        df <- data.frame(sample=template$sample,
                marker_i=template$marker_i,
                marker_j=template$marker_j,
                original_frame_count = template$original_frame_count,
                useful_frame_count = template$useful_frame_count,
                smallest_cell_count = template$smallest_cell_count,
                var=rowQuantiles(combine_var,probs=0.5),
                mean=rowQuantiles(combine_mean,probs=0.5),
                rowQuantiles(combine_mean,probs=quantiles),
                n_resamples = rep(n_resamples,dim(template)[1]),
                check.names=FALSE
                )
        return(df)
    } else if (grouped_sample==FALSE) {
        ### Case 2: We are leaving frames seperate
        df <- data.frame(sample=template$sample,
                         frame=template$frame,
                         marker_i=template$marker_i,
                         marker_j=template$marker_j,
                         smallest_cell_count=template$smallest_cell_count,
                         var=rowQuantiles(combine_var,probs=0.5),
                         mean=rowQuantiles(combine_mean,probs=0.5),
                         rowQuantiles(combine_mean,probs=quantiles),
                         n_resamples=rep(n_resamples,dim(template)[1]),
                         check.names=FALSE
                         )
        return(df)
    }
})

#' Extract the distance to each nearest neighbor for specified
#' cell-types, normalized by downsampling each cell-type to the
#' same size (the smallest population from among the specified
#' markers), calculates across all samples
#'
#' @param data IrisSpatialFeatures ImageSet object
#' @param markers vector of marker names to use
#' @param n_resamples number of times to resample each frame (default:500)
#' @param minimum_cells the smallest number of cells to consider a frame (default:50)
#' @param quantiles vector of numeric fractions to include in vector
#'        to show the mean distance calculated across resamplings
#'        (default:c(0.05,0.25,0.5,0.75,0.95))
#' @param grouped_sample TRUE/FALSE group samples together (default:TRUE)
#'
#' @return data.frame
#'
#' @docType methods
#' @export
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_nearest_neighbor(dataset)
#' normal_nearest_neighbor(dataset,c("SOX10+ PDL1+","SOX10+ PDL1-"),10)
#'
#' @rdname normal_nearest_neighbor
setGeneric("normal_nearest_neighbor", function(data, markers, n_resamples=500,minimum_cells=50,quantiles=c(0.05,0.25,0.5,0.75,0.95),grouped_sample=TRUE)
    standardGeneric("normal_nearest_neighbor"))

#' @rdname normal_nearest_neighbor
#' @aliases normal.nearest.neighbor,ANY,ANY-method
setMethod(
    "normal_nearest_neighbor",
    signature(data="ImageSet"),
    definition <- function(data,markers,n_resamples,minimum_cells,quantiles,grouped_sample) {
    sample_names <- names(data@samples)
    v<-lapply(sample_names,
              normal_nearest_neighbor_sample,
              data=data,
              markers=markers,
              n_resamples=n_resamples,
              quantiles=quantiles,
              minimum_cells=minimum_cells,
              grouped_sample=grouped_sample)
    names(v)<-sample_names
    df <- do.call("rbind",v)
    nnn <- new("NNN")
    nnn@df <- df
    nnn@microns_per_pixel <- data@microns_per_pixel
    return(nnn)
})

#' Class to represent a normalized nearest neighbor.
#'
#' @slot df A dataframe of marker labels and nearest neighbor distances.
#' @slot micros_per_pixel numeric for plotting scale
NNN <- setClass(
    "NNN",
    slots = c(
        df = "data.frame",
        microns_per_pixel = "numeric"
    )
)

#' Get the dataframe from a normalized nearest neighbor object.
#'
#' @param x Normalized nearest neighbor object
#' @param ... Additional arguments
#'
#' @return A dataframe
#' @examples
#' dataset <- IrisSpatialFeatures_data
#' dataframe <- as.data.frame(dataset)
#' @docType methods
#' @export
#' @rdname as.data.frame
setMethod("as.data.frame",
          signature = c(x="NNN"),
          function(x) {
              return(x@df)
          })

setOldClass("htest")
setOldClass("gg")
#' Class to represent a comparison of two markers to a reference.
#'
#' @import ggplot2
#' @slot to_reference_plot A plot comparing markerA and markerB's distance to the reference.
#' @slot to_reference_ttest A paired ttest comparing markerA and markerB's distance to the reference.
#' @slot from_reference_plot A plot comparing markerA and markerB's distance from the reference.
#' @slot from_reference_ttest A ttest comparing markerA and markerB's distance from the reference.
#' @slot to_reference_order For reordering other plots
#' @slot from_reference_order For reordering other plots
#' @slot to_reference_df The raw data
#' @slot from_reference_df The raw data
NNN_compare <- setClass(
    "NNN_Compare",
    slots = c(
        to_reference_plot = "gg",
        to_reference_ttest = "htest",
        from_reference_plot = "gg",
        from_reference_ttest = "htest",
        to_reference_order = "character",
        from_reference_order = "character",
        to_reference_df = "data.frame",
        from_reference_df = "data.frame"
    )
)

#' Compare distances between two markers to a reference marker.
#'
#' @param NNN Normalized nearest neighbor object
#' @param markerA Additional arguments
#' @param markerB Additional arguments
#' @param reference Additional arguments
#' @param order Optional character vector with sample names in an order for plotting
#'
#' @return Analysis data
#' @examples
#' dataset <- IrisSpatialFeatures_data
#' dataframe <- as.data.frame(dataset)
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @import RColorBrewer
#' @docType methods
#' @export
#' @rdname compare_normalized_nearest_neighbor
setGeneric("compare_normalized_nearest_neighbor", function(NNN, markerA, markerB,reference,order=FALSE)
    standardGeneric("compare_normalized_nearest_neighbor"))

#' @rdname compare_normalized_nearest_neighbor
setMethod("compare_normalized_nearest_neighbor",
          signature = c(NNN="NNN",markerA="character",markerB="character",reference="character"),
          function(NNN,markerA,markerB,reference,order=FALSE) {
    t <- as_tibble(as.data.frame(NNN))
    output <- new("NNN_Compare")
    #do_analysis <- function(t,markerA,markerB,reference) {
    font1 <- 8
    font2 <- 8
    pos1 <- t %>% filter(marker_i == markerA & marker_j == reference) %>% select(sample,mean)
    neg1 <- t %>% filter(marker_i == markerB & marker_j == reference) %>% select(sample,mean)
    pos2 <- t %>% filter(marker_j == markerA & marker_i == reference) %>% select(sample,mean)
    neg2 <- t %>% filter(marker_j == markerB & marker_i == reference) %>% select(sample,mean)


    output@to_reference_ttest<-t.test(pos1$mean,neg1$mean,paired=TRUE)
    output@from_reference_ttest<-t.test(pos2$mean,neg2$mean,paired=TRUE)

    # Plot them
    # First reorder factors by distance
    if(class(order)!="logical") {
        ordered_samples = order
    } else {
        ordered_samples <- t %>% filter(marker_i!=reference & marker_j==reference) %>% group_by(sample) %>% summarize(max_val=max(mean)) %>% arrange(desc(max_val))
        ordered_samples = ordered_samples$sample
    }
    #t$sample <- factor(t$sample,levels=ordered_samples$sample)

    # Now plot
    sub = t %>% filter(marker_i == markerA | marker_i == markerB) %>% filter(marker_j==reference)
    output@to_reference_df = sub
    output@to_reference_order = as.vector(ordered_samples)

    output@to_reference_plot <- plot_nnn(sub,markerA,markerB,reference,
                                         paste('From',markerA,'and',markerB,'to',reference),
                                         paste('p =',signif(output@to_reference_ttest$p.value,3)),
                                         paste('to',reference,'from'),
                                         order = as.vector(ordered_samples),
                                         font1,
                                         font2,
                                         NNN@microns_per_pixel)


    if (class(order)!="logical") {
        ordered_samples = order
    } else {
        ordered_samples <- t %>% filter(marker_i==reference & marker_j!=reference) %>% group_by(sample) %>% summarize(max_val=max(mean)) %>% arrange(desc(max_val))
        ordered_samples = ordered_samples$sample
    }
    #switch for other direction
    sub = t %>% filter(marker_j == markerA | marker_j == markerB) %>% filter(marker_i==reference)

    sub2 = sub
    sub2$marker_j<-sub$marker_i
    sub2$marker_i<-sub$marker_j
    output@from_reference_df = sub
    output@from_reference_order = as.vector(ordered_samples)
    output@from_reference_plot <- plot_nnn(sub2,markerA,markerB,reference,
                                         paste('From',reference,'to',markerA,'and',markerB),
                                         paste('p =',signif(output@from_reference_ttest$p.value,3)),
                                         paste('from',reference,'to'),
                                         order = as.vector(ordered_samples),
                                         font1,
                                         font2,
                                         NNN@microns_per_pixel)

    return(output)
})


#' Plot the normalized nearest neighbor (called internally)
#'
#' @param sub subset of data in a tibble
#' @param markerA for factors
#' @param markerB for factors
#' @param reference for factors
#' @param title top of plot
#' @param subtitle on plot
#' @param legendtitle title to put over legend
#' @param order a vector of sample names to display from left to right,
#'              if not specified will order by factor
#' @param font1 font size 1
#' @param font2 font size 2
#' @param microns_per_pixel scale data by this
#'
#' @return gg
#'
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @docType methods
#' @rdname plot_nnn
setGeneric("plot_nnn", function(sub,markerA,markerB,reference,title,subtitle,legendtitle,order,font1,font2,microns_per_pixel)
    standardGeneric("plot_nnn"))

#' @rdname plot_nnn
setMethod(
    "plot_nnn",
    signature(sub="tbl_df"),
    definition <- function(sub,markerA,markerB,reference,title,subtitle,legendtitle,order=c(),font1=6,font2=4,microns_per_pixel=1) {
        #sub = sub %>% filter(marker_i == markerA | marker_i == markerB) %>% filter(marker_j==reference)
        # Get factors in the order we want
        sub$marker_i<-factor(sub$marker_i,levels=c(markerA,markerB))
        sub$marker_j<-factor(sub$marker_j,levels=c(markerA,markerB))
        if(length(order)>0) {
            sub$sample <- factor(sub$sample,levels=order)
        }
        # use our scaling
        sub$mean <- sub$mean * microns_per_pixel
        sub$`75%` <- sub$`75%` * microns_per_pixel
        sub$`25%` <- sub$`25%` * microns_per_pixel
        sub$`95%` <- sub$`95%` * microns_per_pixel
        sub$`5%` <- sub$`5%` * microns_per_pixel
        output<-ggplot(sub,aes(x=factor(marker_i),color=marker_i))+
            geom_boxplot(aes(middle=mean,upper=`75%`,lower=`25%`,ymin=`5%`,ymax=`95%`),
                         stat="identity")+
            facet_wrap(~sample,ncol=26,strip.position="bottom")+
            theme_minimal()+
            theme(axis.text.y = element_text(size=font2),
                  strip.text = element_text(angle=90,size=font1),
                  strip.text.x = element_text(vjust=1),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank())+
            labs(title=title,
                 subtitle=subtitle,
                 color=legendtitle)+
            ylab("microns")
        return(output)
    })
