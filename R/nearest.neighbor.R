


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
#' extract_nearest_neighbor(new("ImageSet"))
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
#' get_all_nearest_neighbors(new("ImageSet"))
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
#' raw_data <- new("ImageSet")
#' raw_data <- read_raw(raw_data,
#'                      raw_dir_name=system.file("extdata", package = "IrisSpatialFeatures"),
#'                      format='Mantra')
#' dataset <- threshold_dataset(raw_data,
#'     marker='PD-Ligand-1 (Opal 690)',
#'     marker_name='PDL1',
#'     base=c('SOX10+'))
#' dataset <- threshold_dataset(dataset,
#'     marker='PD-1 (Opal 540)',
#'     marker_name='PD1',
#'     base=c('CD8+','OTHER'))
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
#' raw_data <- new("ImageSet")
#' raw_data <- read_raw(raw_data,
#'                      raw_dir_name=system.file("extdata", package = "IrisSpatialFeatures"),
#'                      format='Mantra')
#' dataset <- threshold_dataset(raw_data,
#'     marker='PD-Ligand-1 (Opal 690)',
#'     marker_name='PDL1',
#'     base=c('SOX10+'))
#' dataset <- threshold_dataset(dataset,
#'     marker='PD-1 (Opal 540)',
#'     marker_name='PD1',
#'     base=c('CD8+','OTHER'))
#' dataset <- extract_nearest_neighbor(dataset,min_num_cells=2)
#' get_nearest_neighbors(dataset,"SOX10+ PDL1+")
#' p <- plot_nearest_neighbor(dataset,'CD8+ PD1+','SOX10+ PDL1')
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
                          transposed = FALSE) {
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
            current.se <- (t(cbind(x.se, y.se)))
        } else{
            current.mean <- t(x.mean)
            current.se <- t(x.se)
        }

        #remove NA values
        current.se[is.na(current.mean)] <- 0
        current.mean[is.na(current.mean)] <- 0

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

        bp <- barplot(
            current.mean,
            main = buildLabel(from, to, ext, transposed),
            xlab = "Samples",
            ylab = "Avg. distance to NN",
            col = COLS,
            legend = leg,
            ylim = c(0, max(current.mean) + max(current.se)),
            las = 2,
            beside = TRUE
        )

        if (length(comp) > 1) {
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
            pval <-
                t.test(current.mean[1, ], current.mean[2, ], paired = TRUE)$p.value
            mtext(paste('Paired t-test:', format(pval, digits = 4)), 3)
        } else{
            pval <- NA
        }
        return(list(
            means = current.mean,
            ses = current.se,
            pval = pval
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

#####################################
#ray plots


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
#' raw_data <- new("ImageSet")
#' raw_data <- read_raw(raw_data,
#'                      raw_dir_name=system.file("extdata", package = "IrisSpatialFeatures"),
#'                      format='Mantra')
#' dataset <- threshold_dataset(raw_data,
#'     marker='PD-Ligand-1 (Opal 690)',
#'     marker_name='PDL1',
#'     base=c('SOX10+'))
#' dataset <- threshold_dataset(dataset,
#'     marker='PD-1 (Opal 540)',
#'     marker_name='PD1',
#'     base=c('CD8+','OTHER'))
#' dataset <- extract_nearest_neighbor(dataset,min_num_cells=2)
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
            dev.off()
        }
    }
)
normal_nearest_neighbor_sample_once <- function(sample_name,data,markers) {
    # For a single sample designated by sample_name get a dataframe
    markers <- data@markers[data@markers %in% markers]
    sample <- data@samples[sample_name][[1]]
    frame_names <- names(sample@coordinates)
    frame_df_list <- lapply(frame_names,function(frame_name){
        dat <- sample@coordinates[frame_name][[1]]
        # filter down to just the markers we're interested in
        tot <- sapply(markers,function(x){sum(dat@ppp$marks==x)})
        # get the number of cells in each of the categories of interest
        smallest_cell_count <- min(tot)
        # get the number to downsample to
        parr <- lapply(markers,function(x){
            mppp<-dat@ppp[dat@ppp$marks==x,]
            mppp<-mppp[sample(1:length(mppp$marks),smallest_cell_count),]
        })
        names(parr) <- markers
        # exectue downsampling
        nn_df_list <- lapply(markers,function(marki){
            # Get the mean aand variance between all markers a list of lists
            pi <- parr[marki][[1]]
            outs <- lapply(markers,function(markj){
                # Get the mean and variance for nnearest distances bettween markj and marki
                pj <- parr[markj][[1]]
                dis<-spatstat::nncross(pi,pj)[,1]
                list(mean_dist=mean(dis),
                            var_dist=var(dis))
            })
            names(outs) <- markers
            # we want one list of means an one list of means
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
    #now we can get the mean and variance matrix from the mean of the frames
    mean_data <- lapply(frame_names,function(x){
      return(frame_df_list[x][[1]]$mean)
    })
    var_data <- lapply(frame_names,function(x){
      return(frame_df_list[x][[1]]$var)
    })
    populations <- lapply(frame_names,function(x){
      return(frame_df_list[x][[1]]$smallest_cell_count)
    })
    #Combine the frames to get aggrogate statistics of all the frames
    mean_combined <- Reduce("+",mean_data)/length(mean_data)
    var_combined <- Reduce("+",var_data)/length(var_data)
    min_pop <- Reduce("min",populations)
    max_pop <- Reduce("max",populations)
    # Build a data frame with our data
    template <- frame_df_list[1][[1]]
    df <- data.frame(marker_i=template$marker_i,
                marker_j=template$marker_j,
                mean=mean_combined,
                var=var_combined,
                frame_count=rep(length(frame_names),dim(template)[1]),
                min_frame_cells=rep(min_pop,dim(template)[1]),
                max_frame_cells=rep(max_pop,dim(template)[1]),
                sample=rep(sample_name,dim(template)[1])
                )
    return(df)
}
normal_nearest_neighbor_n <- function(sample_name,data,markers,nresample) {
    # For a single sample name, resample it nresample timesand return dataframes
    totals<-lapply(rep(sample_name,nresample),
                   sample_normal_nearest_neighbor_once,
                   data=data,markers=markers)
    combine_mean <- sapply(totals,function(x){x$mean$mean})
    combine_var <- sapply(totals,function(x){x$var$variance})
    template <- totals[[1]]$mean
    mean_df<-cbind(as.data.frame(template$sample_name),
         template$marker_j,
         template$marker_i,
         rowQuantiles(combine_mean,probs=c(0.05,0.5,0.95)))
    var_df<-cbind(as.data.frame(template$sample_name),
         template$marker_j,
         template$marker_i,
         rowQuantiles(combine_var,probs=c(0.05,0.5,0.95)))
    names(mean_df)<-c("sample_name","marker_j","marker_i","5%","50%","95%")
    names(var_df)<-c("sample_name","marker_j","marker_i","5%","50%","95%")  
    return(list(mean=mean_df,var=var_df,nresample=nresample,sample=sample_name))
}
normal_nearest_neighbor <- function(data,markers,nresample) {
    # Pre: take a list of markers and return median nearest neighbor based 
    #  on n resamplings of the the same number of cells for the smallest 
    #  number of cells in the marker list (per-sample)
    # Post: Return a data frame with sample, the i an j markers used in
    #  mean distance, and the 5%, median (50%) and 95% distances
    sample_names <- names(data@samples)
    v<-lapply(sample_names,
            normal_nearest_neighbor_n,
            data=data,
            markers=markers,
            nresample=nresample)
    names(v)<-sample_names
    mean_frame_list<-lapply(sample_names,function(x){
        v[[x]]$mean
    })
    mean_df <- do.call("rbind",mean_frame_list)
    return(mean_df)
}


