
#' DataFrame from all the counts on a per mm2 basis per sample
#'
#' @param x IrisSpatialFeatures ImageSet object.
#' @param ... Additional arguments
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' interactions_sample_data_frame(dataset)
#'
#' @importFrom reshape2 melt
#' @return data frame
#' @docType methods
#' @export
#' @rdname interactions_sample_data_frame
setGeneric("interactions_sample_data_frame",
           function(x, ...)
               standardGeneric("interactions_sample_data_frame"))

#' @rdname interactions_sample_data_frame
#' @aliases  interactions_sample_data_frame,ANY,ANY-method
setMethod(
    "interactions_sample_data_frame",
    signature = "ImageSet",
    definition = function(x) {
        # make sure interactions have been extracted in beforehand
        stop("Use 'interaction_proportion_data_frame' function to access this")
        if (length(x@interactions) == 0) {
            stop(paste(
                'Please run extract.interactions before plotting the interactions.'
            ))
        }
        v <- get_all_interactions(x)
        dfs <- lapply(names(v),function(sample){
            mat <- v[[sample]]
            mat <- melt(mat)
            colnames(mat) <- c('neighbor','reference','proportion')
            mat$sample <- sample
            return(mat)
        })
        dfs <- do.call(rbind,dfs)[,c('sample','reference','neighbor','proportion')]
        return(as.tibble(dfs))
    }
)


#' Extract interactions between all cell-types
#'
#' @param x IrisSpatialFeatures ImageSet object
#' @param membrane_border_width_px the width of the membrane border (default 2)
#' @param interaction_distance_px distance beyond the membrane 
#' @param ... Additional arguments
#' @return list of interactions
#'
#'
#' @docType methods
#' @export
#' @rdname extract_interactions
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_interactions(dataset)
#'
setGeneric("extract_interactions", function(x,
                                            membrane_border_width_px=2,
                                            interaction_distance_px=2, 
                                            ...)
    standardGeneric("extract_interactions"))

#' @rdname extract_interactions
#' @aliases extract_interactions,ANY,ANY-method
setMethod(
    "extract_interactions",
    signature = "ImageSet",
    definition = function(x, 
                          membrane_border_width_px, 
                          interaction_distance_px) {
        x@interactions <-
            lapply(x@samples, 
                          interactions_per_sample, 
                          membrane_border_width_px, 
                          interaction_distance_px, 
                          x@markers)
        names(x@interactions) <- names(x@samples)
        return(x)
    }
)

setGeneric("interactions_per_sample", function(x,
                                            membrane_border_width_px,
                                            interaction_distance_px, 
                                            all_levels, 
                                            ...)
    standardGeneric("interactions_per_sample"))
setMethod(
    "interactions_per_sample",
    signature = "Sample",
    definition = function(x, 
                          membrane_border_width_px, 
                          interaction_distance_px, 
                          all_levels) {
        message(paste(x@sample_name, ' ... processing...'))
        interactions <-
            lapply(x@coordinates, 
                          interaction_events, 
                          membrane_border_width_px, 
                          interaction_distance_px,
                          all_levels)

        areas_with_counts <- sapply(interactions, length) > 0
        if (sum(areas_with_counts) == 0) {
            stop(
                'There has to be at least one coordinate that has cells in it to calculate interactions'
            )
        }

        interactions <- interactions[areas_with_counts]

        ppps <- lapply(interactions, function(x)
            x$ppp)
        ints <- lapply(interactions, function(x)
            x$ints)

        #reshape so we get lists of matrices
        #per_means <-
        #    lapply(interactions, function(x)
        #        x$stats$per$mean)
        #per_vars <-
        #    lapply(interactions, function(x)
        #        x$stats$per$var)
        #avg_means <-
        #    lapply(interactions, function(x)
        #        x$stats$avg$mean)
        #avg_vars <-
        #    lapply(interactions, function(x)
        #        x$stats$avg$var)
        #nums <-
        #    rowSums(sapply(interactions, function(x)
        #        x$stats$nums))

        ##get the total number of interactions
        #total_ints <-
        #    lapply(1:length(per_means), function(x, int)
        #        sweep(int[[x]]$stats$avg$mean, 2, int[[x]]$stats$nums, '*'), interactions)

        #collapse the different coordinates
        #per_means <- collapseMatrices(per_means, rowMeans)
        #per_vars <- collapseMatrices(per_vars, rowMeans)
        #avg_means <- collapseMatrices(avg_means, rowMeans)
        #avg_vars <- collapseMatrices(avg_vars, rowMeans)
        #total_ints <- collapseMatrices(total_ints, rowSums)

        ##calculate the standard error on the combined coordinates
        #per_SE <- sweep(sqrt(per_vars), 2, sqrt(nums), '/')
        #avg_SE <- sweep(sqrt(avg_vars), 2, sqrt(nums), '/')

        return(list(
            #per = list(mean = per_means,
            #           SE = per_SE),
            #avg = list(mean = avg_means,
            #           SE = avg_SE),
            #total = total_ints,
            ppp = ppps,
            ints = ints
            #nums = nums
        ))
    }
)

#' DataFrame from all the interaction counts on a per mm2 basis per frame grouped by sample
#' This function is more like a density of edges between nodes on the interaction graph
#' than an average neighbor count
#'
#' @param x IrisSpatialFeatures ImageSet object.
#' @param ... Additional arguments
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' interactions_sample_data_frame(dataset)
#'
#' @return data frame
#' @import magrittr dplyr tibble
#' @docType methods
#' @export
#' @rdname interaction_counts_sample_data_frame
setGeneric("interaction_counts_sample_data_frame",
           function(x, 
                    verbose = FALSE,
                    ...)
               standardGeneric("interaction_counts_sample_data_frame"))

#' @rdname interaction_counts_sample_data_frame
#' @aliases  interaction_counts_sample_data_frame,ANY,ANY-method
setMethod(
    "interaction_counts_sample_data_frame",
    signature = "ImageSet",
    definition = function(x, 
                    verbose) {
        if (length(x@interactions) == 0) {
            stop(paste(
                'Please run extract.interactions before plotting the interactions.'
            ))
        }
        sresults <- interaction_counts_data_frame(x, 
                                                  verbose)
        sresults <- sresults %>% group_by(sample,reference_phenotype,neighbor_phenotype) %>% summarize(total_interactions=sum(interaction_count),mean_interactions_per_mm2=mean(interactions_per_mm2),stderr_interactions_per_mm2=sd(interactions_per_mm2)/sqrt(n()),frame_count=n())
        return(sresults)
    }
)

#' DataFrame from all the interaction counts on a per mm2 basis per frame
#' This function is more like a density of edges between nodes on the interaction graph
#' than an average neighbor count
#'
#' @param x IrisSpatialFeatures ImageSet object.
#' @param ... Additional arguments
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' interactions_sample_data_frame(dataset)
#'
#' @return data frame
#' @import magrittr dplyr tibble
#' @docType methods
#' @export
#' @rdname interaction_counts_data_frame
setGeneric("interaction_counts_data_frame",
           function(x, 
                    verbose=FALSE, 
                    ...)
               standardGeneric("interaction_counts_data_frame"))

#' @rdname interaction_counts_data_frame
#' @aliases  interaction_counts_data_frame,ANY,ANY-method
setMethod(
    "interaction_counts_data_frame",
    signature = "ImageSet",
    definition = function(x, 
                    verbose) {
        if (length(x@interactions) == 0) {
            stop(paste(
                'Please run extract.interactions before plotting the interactions.'
            ))
        }

        sresults <- lapply(names(x@samples),function(sample_name){
            if (verbose) { print(sample_name) }
            sample <- x@samples[[sample_name]]
            fresults <- lapply(names(sample@coordinates),function(frame_name){
                frame <- sample@coordinates[[frame_name]]
                if (verbose) { print(frame_name) }

                marks = as.character(x@interactions[[sample_name]]$ppp[[frame_name]]$marks)
                if (length(x@interactions[[sample_name]]$ints[[frame_name]])==0) {
                    print("zero interactions case. whole frame will not contribute. if you want to force zero counts for these, add this.")
                    return(NULL)
                }
                ## Lets set up a structure to hold the converted data
                phenotypes <- levels(x@samples[[sample_name]]@coordinates[[frame_name]]@ppp$marks)
                #print(phenotypes)
                cnts = list()
                for (phenotype_name1 in phenotypes) {
                    cnts[[phenotype_name1]] == list()
                    for (phenotype_name2 in phenotypes) {
                        cnts[[phenotype_name1]][[phenotype_name2]] = 0
                    }
                }
                for (i in seq(1,length(x@interactions[[sample_name]]$ints[[frame_name]]),1)) {
                    reference_phenotype = marks[i]
                    if (is.null(x@interactions[[sample_name]]$ints[[frame_name]][[i]])) { next }
                    count_table <- table(sapply(x@interactions[[sample_name]]$ints[[frame_name]][[i]],function(x){return(marks[x])}))
                    for (neighbor_phenotype in names(count_table)) {
                        cnts[[reference_phenotype]][[neighbor_phenotype]] = count_table[[neighbor_phenotype]] + cnts[[reference_phenotype]][[neighbor_phenotype]]
                    }

                }
                v1 <- lapply(phenotypes,function(reference_phenotype) {
                    v2 <- lapply(phenotypes,function(neighbor_phenotype) {
                        df <- data.frame(reference_phenotype = reference_phenotype,
                            neighbor_phenotype=neighbor_phenotype,
                            interaction_count=cnts[[reference_phenotype]][[neighbor_phenotype]])
                        colnames(df) <- c('reference_phenotype','neighbor_phenotype','interaction_count')
                        return(df)
                    })
                    v2 <- do.call(rbind,v2)
                    return(v2)
                })
                cnts <- do.call(rbind,v1)

                # Check for zero remainder
                values = list()
                #print("cycle markers")
                z = 0
                for (marki in x@markers) {
                    for (markj in x@markers) {
                        if (dim(cnts %>% filter(neighbor_phenotype==markj & reference_phenotype==marki))[1]==0) {
                            z = z + 1
                            values[[z]] <- c(markj,marki)
                        }
                    }
                }
                #print(cnts)
                #print("rbind")
                if(length(values)>0) {
                    values <- as.data.frame(do.call(rbind,values))
                    colnames(values) <- c("neighbor_phenotype","reference_phenotype")
                    values$interaction_count = 0
                    #print(values)
                    #print(cnts)
                    #cnts <- cnts %>% bind_rows(values)
                    cnts <- rbind(cnts,values)
                }
                cnts$size_in_px <- frame@size_in_px
                cnts$frame <- frame_name
                cnts$frame <- as.character(cnts$frame)
                cnts$sample <- sample_name
                #print(cnts)
                return(cnts)
            })
            #print(fresults)
            fresults <- do.call(rbind,fresults)
            return(fresults)
        })
        sresults <- do.call(rbind,sresults)
        sresults <- sresults[,c('sample','frame','size_in_px','reference_phenotype','neighbor_phenotype','interaction_count')]
        sresults <- sresults %>% mutate(interactions_per_mm2=1000000*interaction_count/(size_in_px*(x@microns_per_pixel ^2)))
        return(sresults)
    }
)

#' DataFrame from all the interaction proportions freshly calculated from our current extracted counts
#'
#' @param x IrisSpatialFeatures ImageSet object.
#' @param ... Additional arguments
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' interactions_sample_data_frame(dataset)
#'
#' @return data frame
#' @import magrittr dplyr tibble
#' @docType methods
#' @export
#' @rdname interaction_proportion_data_frame
setGeneric("interaction_proportion_data_frame",
           function(x, 
                    verbose=FALSE, 
                    ...)
               standardGeneric("interaction_proportion_data_frame"))

#' @rdname interaction_proportion_data_frame
#' @aliases  interaction_proportion_data_frame,ANY,ANY-method
setMethod(
    "interaction_proportion_data_frame",
    signature = "ImageSet",
    definition = function(x, 
                    verbose) {
        if (length(x@interactions) == 0) {
            stop(paste(
                'Please run extract.interactions before plotting the interactions.'
            ))
        }
        
        sresults <- lapply(names(x@samples),function(sample_name){
            if (verbose) { print(sample_name) }
            ## Initialize the place to store data
            phenotypes <- levels(x@samples[[sample_name]]@coordinates[[1]]@ppp$marks)
            cnts <- list()
            totals <- list()
            for (phenotype_name1 in phenotypes) {
                cnts[[phenotype_name1]] <- list()
                totals[[phenotype_name1]] <- 0
                for (phenotype_name2 in phenotypes) {
                    cnts[[phenotype_name1]][[phenotype_name2]] <- 0
                }
            }
            sample <- x@samples[[sample_name]]
            for (frame_name in names(sample@coordinates)) {
                frame <- sample@coordinates[[frame_name]]
                if (verbose) { print(frame_name) }
                marks = as.character(x@interactions[[sample_name]]$ppp[[frame_name]]$marks)
                if (length(x@interactions[[sample_name]]$ints[[frame_name]])==0) {
                    print("zero interactions case. whole frame will not contribute. if you want to force zero counts for these, add this.")
                    next
                }
                cell_numbers = seq(1,length(x@interactions[[sample_name]]$ints[[frame_name]]),1)
                for (i in cell_numbers) {
                    reference_phenotype = marks[i]
                    if (is.null(x@interactions[[sample_name]]$ints[[frame_name]][[i]])) { next } # skip cells with no neighbors
                    count_table <- table(sapply(x@interactions[[sample_name]]$ints[[frame_name]][[i]],function(x){return(marks[x])})) # the number of each neighbor
                    totals[[reference_phenotype]] = totals[[reference_phenotype]] + length(x@interactions[[sample_name]]$ints[[frame_name]][[i]])
                    for (neighbor_phenotype in names(count_table)) {
                        cnts[[reference_phenotype]][[neighbor_phenotype]] = cnts[[reference_phenotype]][[neighbor_phenotype]] + count_table[[neighbor_phenotype]]
                    }
                }
            }
            # Now assemble a dataframe for the sample
            v1 <- lapply(phenotypes,function(phenotype1) {
                v2 <- lapply(phenotypes,function(phenotype2) {
                    row <- list(sample=sample_name,reference_phenotype=phenotype1,neighbor_phenotype=phenotype2,interactions=0,count=0,proportion=NULL)
                    if (totals[[phenotype1]] > 0) {
                        row[['interactions']] = totals[[phenotype1]]
                        row[['count']] = cnts[[phenotype1]][[phenotype2]]
                        row[['proportion']] = cnts[[phenotype1]][[phenotype2]]/row[['interactions']]
                    } else { 
                        return(NULL)
                    }
                    return(data.frame(row))
                })
                rows <- do.call(rbind,v2)
            })
            return(do.call(rbind,v1))
        })
        sresults <- do.call(rbind,sresults)
        return(as.tibble(sresults))
    }
)

setGeneric("interaction_events", function(x,
                                          membrane_border_width_px,
                                          interaction_distance_px, 
                                          all_levels, ...)
    standardGeneric("interaction_events"))
setMethod(
    "interaction_events",
    signature = "Coordinate",
    definition = function(x,
                          membrane_border_width_px,
                          interaction_distance_px, 
                          all_levels) {
        #extract membrane map and set membranes to -1
        if (is.null(x@raw@mem_seg_map) ||
            any(dim(x@raw@mem_seg_map) == 0)) {
            stop(
                'The interaction analysis can only be run on datasets that include the membrane maps. Try the proximity analysis instead.'
            )
        }

        if (length(x@ppp$x) == 0) {
            return(list())
        }

        #fill in all of the cells in the membrane map using the cell ID
        ret <- watershed(x)

        #update the values
        filled_map <- ret$map
        #saveRDS(filled_map,'test.RDS')
        #writeTIFF(t(as.matrix(filled_map)),'test.tiff',compression='LZW')
        #print(paste0(membrane_border_width_px," ",interaction_distance_px))
        x <- ret$x

        #extract the interactions
        interactions <- getNeighbors(filled_map,
                                     membrane_border_width_px,
                                     interaction_distance_px)
        #print(interactions)
        #stop()

        #extract the means and variances
        #inter_stats <-
        #    extract_interaction_stats(x, interactions, all_levels)

        return(list(
            #stats = inter_stats,
            ppp = x@ppp,
            ints = interactions
        ))
    }
)


setGeneric("watershed", function(x, ...)
    standardGeneric("watershed"))
setMethod(
    "watershed",
    signature = "Coordinate",
    definition = function(x) {
        mem_map <- t(x@raw@mem_seg_map)
        mem_map[mem_map > 0] <- -1
        #watershed filling in all cells with their cell ID
        padded_map <- rbind(-1, cbind(-1, mem_map, -1), -1)

        #need to offset the coordinates because of the padding
        cell_coords <- cbind(1:length(x@ppp$x),
                             x@ppp$x,
                             x@ppp$y)

        #run the watershed algorithm and fill up all cells
        ret <- watershedC(padded_map, cell_coords)

        #filled in cells
        padded_map <- ret[[1]]

        #updated coordinates
        x@ppp$x <- ret[[2]][, 2]
        x@ppp$y <- ret[[2]][, 3]

        #remove the padding
        padded_map <-
            padded_map[-c(1, nrow(padded_map)), -c(1, ncol(padded_map))]

        return(list(map = padded_map, x = x))
    }
)

#' @importFrom stats var
get_single_int <- function(lvl, int, labels, all_levels) {
    #generate a per cell summary
    ints <-
        lapply(int[labels == lvl], function(i, lev)
            factor(i, levels = lev), all_levels)

    if (length(ints) == 0) {
        num_cells <- 0
        per_means <- rep(NA, length(all_levels))
        names(per_means) <- all_levels
        avg_means <- per_vars <- avg_vars <- per_means
    } else{
        per_cell_summary <- t(sapply(ints, table))

        num_cells <- nrow(per_cell_summary)

        #calculate the average interaction measurements
        avg_means <- colMeans(per_cell_summary, na.rm = TRUE)
        avg_vars <- apply(per_cell_summary, 2, var, na.rm = TRUE)


        #calculate the percent interaction measurements
        per_cell_summary <- per_cell_summary > 0
        per_means <- colMeans(per_cell_summary, na.rm = TRUE)
        per_vars <- apply(per_cell_summary, 2, var, na.rm = TRUE)
    }


    return(list(
        per = list(mean = per_means, var = per_vars),
        avg = list(mean = avg_means, var = avg_vars),
        nums = num_cells
    ))
}

setGeneric("extract_interaction_stats", function(x, ...)
    standardGeneric("extract_interaction_stats"))
setMethod(
    "extract_interaction_stats",
    signature = "Coordinate",
    definition = function(x, interactions, all_levels) {
        print("WARNING: This function 'extract_interaction_stats' is depricated")
        labels <- as.character(x@ppp$marks)

        #translate the coordinates in the lists to labels
        int <- lapply(interactions, function(x, lab)
            lab[x], labels)

        #get all the stats
        stats <-
            lapply(all_levels, get_single_int, int, labels, all_levels)

        #reshape the data so we get 4 matrices and 1 vector
        per_means <- sapply(stats, function(x)
            x$per$mean)
        per_vars <- sapply(stats, function(x)
            x$per$var)
        avg_means <- sapply(stats, function(x)
            x$avg$mean)
        avg_vars <- sapply(stats, function(x)
            x$avg$var)
        nums <- sapply(stats, function(x)
            x$nums)

        #make sure we get the right column names
        colnames(per_means) <- colnames(per_vars) <-  all_levels
        colnames(avg_means) <- colnames(avg_vars) <-  all_levels
        names(nums) <- all_levels

        return(list(
            per = list(mean = per_means, var = per_vars),
            avg = list(mean = avg_means, var = avg_vars),
            nums = nums
        ))
    }
)


getNeighbors <- function(filled_map,
                         membrane_border_width_px,
                         interaction_distance_px) {
    step <- membrane_border_width_px+interaction_distance_px
    interactions <- getInteractionsC(filled_map,step)[[1]]
    interactions <- interactions[-(1:2)]
    #transform the list so the indices correspond to the names
    interactions <-
        lapply(as.character(1:max(filled_map)), function(x, int)
            int[[x]], interactions)
    return(interactions)
}

collapseMatrices <- function(mat, fun) {
    #  Make a 3D array from list of matrices
    arr <-
        array(unlist(mat) , c(nrow(mat[[1]]), nrow(mat[[1]]), length(mat)))
    #  Get mean of third dimension
    collapsed <- fun(arr , dims = 2 , na.rm = TRUE)
    colnames(collapsed) <- rownames(collapsed) <- colnames(mat[[1]])
    return(collapsed)
}

################################################################
##### Interaction getters

#' Get all interactions between all cell-types
#'
#' @param x An IrisSpatialFeatures ImageSet object.
#' @param ... Additional arguments.
#' @return For each cell-type return interactions
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_interactions(dataset)
#' get_all_interactions(dataset)
#'
#' @docType methods
#' @export
#' @rdname get_all_interactions
setGeneric("get_all_interactions", function(x, ...)
    standardGeneric("get_all_interactions"))

#' @rdname get_all_interactions
#' @aliases get_all_interactions,ANY,ANY-method
setMethod(
    "get_all_interactions",
    signature = "ImageSet",
    definition = function(x) {
        ints <- lapply(x@interactions, function(y)
            y$avg$mean)
        int_norm <-
            lapply(ints, function(y)
                sweep(y, 2, colSums(y), '/'))
        return(int_norm)
    }
)

#' Get interactions for a specific marker
#'
#' @param x An IrisSpatialFeatures ImageSet object
#' @param marker Cell-type for which the interactions should be pulled
#' @param normalize Flag to indicated whether to normalize each sample so all interactions sum up to 1 (Default: 1)
#' @param ... Additional arguments.
#' @return interactions for a specific marker
#'
#' @docType methods
#' @export
#' @rdname get_interactions
#' @examples
#'
#' #' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_interactions(dataset)
#' get_interactions(dataset,'CD8+ PD1+')
#'
setGeneric("get_interactions", function(x, ...)
    standardGeneric("get_interactions"))

#' @rdname get_interactions
#' @aliases get_interactions,ANY,ANY-method
setMethod(
    "get_interactions",
    signature = "ImageSet",
    definition = function(x, marker, normalize = TRUE) {
        if (!marker %in% x@markers) {
            stop(paste('There is no celltype: ', marker))
        }

        int <- lapply(x@interactions, function(x)
            x$avg$mean)
        marker_int <- sapply(int, function(x)
            x[, marker])

        if (normalize) {
            marker_int <- sweep(marker_int, 2, colSums(marker_int), '/')
        }
        return(marker_int)
    }
)

################################################################
##### Interaction summary plotting functions

#' Interaction summary plot for all cell-types and all samples in a dataset
#'
#' @param x IrisSpatialFeatures ImageSet object to be plotted
#' @param label The cell type the interaction profile should be plotted for
#' @param ordering Ordering of the samples (default: NULL)
#' @param normalize Normalize the interactions with a given cell-type, so they sum up to 1 (default: TRUE)
#' @param palette Color palette for all the cell-types (default: Spectral scheme from RColorbrewer)
#' @param celltype_order Order in which the cell-types are displayed. (default: Alphabethically)
#' @param xlim_fix Whitespace on the right side so the legend can be displayed clearly. (default: 13)
#' @param topbar_cols Color of the barplots that are shown on top. (default: 'darkgrey')
#' @param ... Additional arguments
#' @return plot of all cell-types and samples interactions
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_interactions(dataset)
#' plot_interactions(dataset,'SOX10+ PDL1+',xlim_fix=3)
#'
#' @importFrom graphics axis
#' @importFrom graphics layout
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics text
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#'
#' @docType methods
#' @export
#' @rdname plot_interactions
setGeneric("plot_interactions", function(x, ...)
    standardGeneric("plot_interactions"))

#' @rdname plot_interactions
#' @aliases plot_interactions,ANY,ANY-method
setMethod(
    "plot_interactions",
    signature = "ImageSet",
    definition = function(x,
                          label,
                          ordering = NULL,
                          normalize = TRUE,
                          palette = NULL,
                          celltype_order = NULL,
                          xlim_fix = 13,
                          topbar_cols = 'darkgrey') {
        if (length(x@interactions) == 0) {
            stop(paste(
                'Please run extract.interactions before plotting the interactions.'
            ))
        }

        int <- lapply(x@interactions, function(x)
            x$avg$mean)
        dat <- sapply(int, function(x)
            x[, label])
        count <- get_counts_per_mm2(x, blank = TRUE)[label, ]
        labels <- rownames(dat)

        if (normalize) {
            dat <- sweep(dat, 2, colSums(dat), '/')
            ylab <- 'Proportions of interactions'

        } else{
            ylab <- 'Average interactions'
        }

        if (!is.null(ordering)) {
            if (length(ordering) == 1) {
                ordering <-
                    order(colSums(dat[grep(ordering,
                                           rownames(dat),
                                           fixed = TRUE), ]),
                          decreasing = TRUE)
            }
            dat <- dat[, ordering]
            count <- count[ordering]
        }

        if (!is.null(celltype_order)) {
            dat <- dat[celltype_order, ]
        } else{
            celltype_order <- rownames(dat)
        }

        if (is.null(palette)) {
            palette <- brewer.pal(length(labels), "Spectral")
        }

        #generate the plots
        op <- par(no.readonly = TRUE)
        layout(
            as.matrix(2:1),
            widths = c(1),
            heights = c(0.4, 0.8),
            respect = FALSE
        )

        par(mar = c(6, 4, 0, 0))
        bp <- barplot(
            dat,
            cex.names = 1,
            # makes x-axis labels small enough to show all
            col = palette,
            # colors
            xlab = '',
            ylab = ylab,
            las = 2,
            xaxt = "n",
            xlim = c(0, ncol(dat) + xlim_fix),
            width = 1
        )
        text(
            cex = 1,
            x = bp + 0.8,
            y = -0.05,
            colnames(dat),
            xpd = TRUE,
            srt = 45,
            pos = 2
        )

        legend("right",
               legend = celltype_order,
               fill = palette)

        par(mar = c(0.5, 4, 4, 0))
        p <- barplot(
            count,
            xlim = c(0, ncol(dat) + xlim_fix),
            col = topbar_cols,
            axisnames = FALSE,
            axes = FALSE,
            cex.names = 0.5,
            main = paste('Interactions with', label)
        )
        mtext('Counts / mm2', side = 2, line = 2)
        axis(
            side = 2,
            tick = TRUE,
            labels = TRUE,
            line = -1,
            las = 1,
            cex.axis = 0.5
        )
        par(op)
        return(dat)
    }
)

################################################################
##### Interaction maps


#' Plot interaction maps for all samples
#' @param x An IrisSpatialFeatures ImageSet object
#' @param int_markers Cell-types that should be considered. If two cells
#'                    from different cell-types interact they are filled in,
#'                    if a cell is not interacting it is just outlined.
#' @param int_marker_cols Colors for the cell-types
#' @param silent_markers Cell-types that should only be outlined (Default: c())
#' @param silent_col Colors for silent markers (Default: c())
#' @param outline_transparency Dimming factor for the outlines cells(Default: 0.9)
#' @param use_dapi Use the DAPI channel as a background (Default: FALSE)
#' @param outdir Output directory (Default: './interaction_maps')
#' @param useMask (Default: NULL)
#' @param format Output format of the images. Can be '.png' or '.tiff' (Default: '.png')
#' @param ... Additional arguments.
#' @return plot of interactions for all samples
#'
#' @docType methods
#' @export
#' @rdname interaction_maps
#' @examples
#'
#' #' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_interactions(dataset)
#' get_interactions(dataset,'CD8+ PD1+')
#' int_markers <- c('CD8+ PD1+','SOX10+ PDL1+')
#' int_marker_cols <- c('#dd1c77','#99d8c9')
#' silent_markers <- c('CD8+ PD1-')
#' silent_col=c('yellow')
#' p <- interaction_maps(dataset,int_markers,int_marker_cols,silent_markers,
#'                       silent_col)
setGeneric("interaction_maps", function(x, ...)
    standardGeneric("interaction_maps"))

#' @rdname interaction_maps
#' @aliases interaction_maps,ANY,ANY-method
setMethod(
    "interaction_maps",
    signature = "ImageSet",
    definition = function(x,
                          int_markers,
                          int_marker_cols,
                          silent_markers = c(),
                          silent_col = c(),
                          outline_transparency = 0.9,
                          use_dapi = FALSE,
                          outdir = 'interaction_maps',
                          useMask = NULL,
                          format = '.png') {
        if (length(x@interactions) == 0) {
            stop('Please run "extract_interactions" before plotting the interaction maps.')
        }

        #generate the mapping directory
        #map_dir <- file.path(getwd(), outdir)
        map_dir <-outdir
        if (!file.exists(map_dir)) {
            dir.create(map_dir, showWarnings = FALSE)
        }

        #generate a map for each sample
        lapply(
            x@samples,
            interaction_map_sample,
            x@interactions,
            int_markers,
            int_marker_cols,
            silent_markers,
            silent_col,
            map_dir,
            outline_transparency,
            use_dapi,
            useMask,
            format
        )
        return('Done!')

    }
)

setGeneric("interaction_map_sample", function(x, ...)
    standardGeneric("interaction_map_sample"))
setMethod(
    "interaction_map_sample",
    signature = "Sample",
    definition = function(x,
                          interactions,
                          int_markers,
                          int_marker_cols,
                          silent_markers,
                          silent_col,
                          map_dir,
                          outline_transparency,
                          use_dapi,
                          useMask,
                          format) {
        message("Working on sample: ", x@sample_name)
        lapply(
            x@coordinates,
            generate_interaction_map,
            x@sample_name,
            interactions[[x@sample_name]],
            int_markers,
            int_marker_cols,
            silent_markers,
            silent_col,
            map_dir,
            outline_transparency,
            use_dapi,
            useMask,
            format
        )
    }
)


#' @importFrom gplots colorpanel
#' @importFrom graphics image
#' @importFrom graphics legend
#' @importFrom grDevices col2rgb
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom spatstat rgb2hex
#' @importFrom tiff writeTIFF
setGeneric("generate_interaction_map", function(x, ...)
    standardGeneric("generate_interaction_map"))
setMethod(
    "generate_interaction_map",
    signature = "Coordinate",
    definition = function(x,
                          samp_name,
                          interactions,
                          int_markers,
                          int_marker_cols,
                          silent_markers,
                          silent_col,
                          map_dir,
                          outline_transparency,
                          use_dapi,
                          useMask,
                          format) {
        nams <- paste(samp_name, x@coordinate_name, sep = '_')
        #extract all data
        int <- interactions$ints[[x@coordinate_name]]
        ppp <-
            interactions$ppp[[x@coordinate_name]]

        #get the marker prefix
        marker_prefix <-
            paste0(c(int_markers, silent_markers), collapse = '__')
        #remove spaces from prefix so R CMD check does not mope
        marker_prefix <- gsub(' ', '_', marker_prefix)

        #extract membrane map and set membranes to -1
        if (is.null(x@raw@mem_seg_map)) {
            stop(
                'The interaction maps can only be created on datasets that include the membrane maps.'
            )
        }
        mem <- t(x@raw@mem_seg_map)
        mem[mem > 0] <- -1

        if (use_dapi) {
            #extract membrane map and set membranes to -1
            if (is.null(x@raw@dapi_map)) {
                stop('No DAPI map available, please set dapi_map flag to FALSE')
            }
            dapi_map <- t(as.matrix(x@raw@dapi_map))
            dapi_map <- dapi_map / max(dapi_map)
        } else{
            dapi_map <- mem
            dapi_map[dapi_map != 0] <- 0
        }

        if (all(int_markers %in% ppp$marks)) {
            #generate the masks
            #first the interaction ones
            int_marker_masks <-
                lapply(int_markers, generate_mask, mem, ppp)
            names(int_marker_masks) <- int_markers

            #then the ones we don't fill out
            sil_marker_masks <-
                lapply(silent_markers, generate_mask, mem, ppp)
            names(sil_marker_masks) <- silent_markers

            #fill in the marker masks that are relevant for interactions
            int_marker_masks <-
                lapply(int_markers,
                       fill_in_maps,
                       int_markers,
                       int_marker_masks,
                       mem,
                       int,
                       ppp)

            #make the outlines transparent to increase the contrast
            int_marker_cols2 <-
                unlist(lapply(int_marker_cols, function(x)
                    c(
                        rgb2hex(col2rgb(x)[, 1] * outline_transparency), x
                    )))
            silent_col <-
                unlist(lapply(silent_col, function(x)
                    rgb2hex(col2rgb(x)[, 1] * outline_transparency)))

            #color panel
            cols <-
                c(colorpanel(10, 'black', 'blue'),
                  silent_col,
                  int_marker_cols2)
            breaks <-
                c(seq(0, 1, 0.1), (1:(
                    length(silent_col) + length(int_marker_cols2)
                )) + 1.5)

            #add up all the outlines
            if (length(sil_marker_masks) > 0) {
                for (i in 1:length(sil_marker_masks)) {
                    dapi_map[sil_marker_masks[[i]] == 1] <- i + 1

                }
            }

            col_count <- length(sil_marker_masks)
            #add up all the masks that have interactions
            for (i in 1:length(int_marker_masks)) {
                dapi_map[int_marker_masks[[i]] == 1] <- col_count + i + 1

                dapi_map[int_marker_masks[[i]] == 2] <-
                    col_count + i + 2

                col_count <- col_count + 1
            }

            if (format == '.png') {
                png(
                    file.path(map_dir, paste0(
                        nams, '_', marker_prefix, '.png'
                    )),
                    width = nrow(dapi_map),
                    height = ncol(dapi_map)
                )
                image(
                    dapi_map[, ncol(dapi_map):1],
                    col = cols,
                    breaks = breaks,
                    yaxt = 'n',
                    xaxt = 'n'
                )
                legend(
                    'bottomleft',
                    c(int_markers, silent_markers),
                    col = c(int_marker_cols, silent_col),
                    cex = 1.5,
                    pch = 18
                )
                dev.off()
            } else if (format == '.tiff') {
                img <- t(int_marker_masks[[1]])
                tif <- array(0, dim = c(nrow(img), ncol(img), 4))
                tif[, , 4] <- 0

                #actual markers
                for (i in 1:length(int_marker_masks)) {
                    #if a mask was specified show the
                    #interactions only within those masks
                    if (!is.null(useMask)) {
                        int_marker_masks[[i]][x@mask[[useMask]] == 0] <- 0
                    }

                    current_col <-
                        col2rgb(int_marker_cols[i])[, 1] / 255
                    for (j in 1:3) {
                        tmp <- tif[, , j]
                        tmp[t(int_marker_masks[[i]]) > 0] <-
                            current_col[j]
                        tif[, , j] <- tmp
                    }
                    tmp <- tif[, , 4]
                    tmp[t(int_marker_masks[[i]]) == 1] <-
                        outline_transparency
                    tmp[t(int_marker_masks[[i]]) == 2] <- 1
                    tif[, , 4] <- tmp
                }

                #silent markers
                if (length(sil_marker_masks) > 0) {
                    for (i in 1:length(sil_marker_masks)) {
                        current_col <- col2rgb(silent_col[i])[, 1] / 255
                        #if a mask was specified show the
                        #interactions only within those masks
                        if (!is.null(useMask)) {
                            sil_marker_masks[[i]][x@mask[[useMask]] == 0] <- 0
                        }
                        for (j in 1:3) {
                            tmp <- tif[, , j]
                            tmp[t(sil_marker_masks[[i]]) > 0] <-
                                current_col[j]
                            tif[, , j] <- tmp
                        }
                        tmp <- tif[, , 4]
                        tmp[t(sil_marker_masks[[i]]) == 1] <-
                            outline_transparency
                        tif[, , 4] <- tmp
                    }
                }
                writeTIFF(tif, file.path(map_dir, paste0(
                    nams, '_', marker_prefix, '.tiff'
                )))

            }
        }
    }
)

fill_in_maps <- function(marker,
                         markers,
                         marker_masks,
                         mem,
                         int,
                         ppp) {
    mask <- marker_masks[[marker]]
    others <- markers[markers != marker]
    labels <- as.character(ppp$marks)

    #translate the interactions into actual labels
    inter <- lapply(int, function(x, lab)
        lab[x], labels)

    #look only at the cells from the current marker
    indicator <- labels == marker

    #focus only on the cells with the current marker
    inter <- inter[indicator]
    coords <- ppp[indicator, ]

    #coordinates to fill
    inter <-
        sapply(inter, function(x, others)
            sum(x %in% others) > 0, others)
    coords <- data.frame(coords)[inter, 1:2]
    coords <- as.matrix(coords)

    if (nrow(coords) > 0) {
        #current mask
        mask <- rbind(-1, cbind(-1, mask, -1), -1)
        padded <- rbind(-1, cbind(-1, mem, -1), -1)
        #c code
        mask <- fillMaskC(mask, padded, coords)
        #remove the padding
        mask <- mask[-c(1, nrow(mask)), -c(1, ncol(mask))]
    }
    return(mask)
}

#generates a mask for a single marker
generate_mask <- function(lvl, mem, ppp) {
    #add padding to the membrane map so the watershed
    #cannot go outside of the boundaries
    padded_map <- rbind(-1, cbind(-1, mem, -1), -1)

    #extract the coordinates for the current marker
    marker_coords <- ppp[ppp$marks == lvl, ]
    #adust the cell coordinates to take the padding into account
    cell_coords <- cbind(marker_coords$x, marker_coords$y)

    #make a blank mask
    marker_map <- padded_map
    marker_map[marker_map != 0] <- 0
    marker_map <-
        generate_maskC(marker_map, padded_map, cell_coords)

    #remove the padding
    marker_map <-
        marker_map[-c(1, nrow(marker_map)), -c(1, ncol(marker_map))]

    return(marker_map)

}

############ Permutation Test Interactions ##############
#' Calculate a permutation test result for nearest neighbors to say for each sample to see 
#' if the neighbor distance something seen under the null assumption
#'
#' @param x IrisSpatialFeatures ImageSet object that has had extract nearest neighbors run
#' @param permutations Set to 100 by default
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
#' interaction_permutation_test(dataset)
#'
#' @importFrom data.table rbindlist
#' @rdname interaction_permutation_test
setGeneric("interaction_permutation_test", function(x,permutations=20,subset=NULL,verbose=FALSE,...)
    standardGeneric("interaction_permutation_test"))

#' @rdname interaction_permutation_test
#' @aliases interaction_permutation_test,ANY,ANY-method
setMethod(
    "interaction_permutation_test",
    signature = c("ImageSet"),
    definition = function(x, permutations,subset,vebose) {
        obs <- as.tibble(interaction_counts_sample_data_frame(x))
        if (verbose) { print("finished observed") }
        expected <- lapply(seq(1,permutations),function(i){
            if (verbose) { print(i) }
            #datar1 <- extract_nearest_neighbor(s)
            vr <- as.tibble(interaction_counts_sample_data_frame(shuffle_labels(x,subset=subset)))
            vr$iter <- i
            return(vr)
        })
        expected <- rbindlist(expected)

        ci <- expected %>% group_by(sample, reference_phenotype, neighbor_phenotype) %>% summarize(`5%`=quantile(total_interactions,probs=0.05),
                                      `95%`=quantile(total_interactions,probs=0.95), expected_total_mean=mean(total_interactions),expected_total_sd=sd(total_interactions))
        annot <- obs %>% full_join(ci,by=c('sample','reference_phenotype','neighbor_phenotype')) 
        annot$z_score <- (annot$total_interactions-annot$expected_total_mean)/(annot$expected_total_sd/sqrt(permutations))
        annot$permutations = permutations

        # now get the p value
        subset = annot %>% rename(observed = total_interactions) %>% select(sample,reference_phenotype, neighbor_phenotype, observed,z_score,permutations)
        r2 = as.tibble(expected) %>% full_join(subset,by=c('sample','reference_phenotype','neighbor_phenotype'))
        low_values = annot %>% filter(z_score <=0)
        high_values = annot %>% filter(z_score > 0)
        low = r2 %>% filter(z_score <=0) %>% filter(total_interactions <= observed)
        high = r2 %>% filter(z_score >0) %>% filter(total_interactions >= observed)
        hcnt = high %>% group_by(sample, reference_phenotype, neighbor_phenotype) %>% summarize(count=n())
        high_values = high_values %>% full_join(hcnt,by=c('sample','reference_phenotype','neighbor_phenotype')) %>% mutate(p_value=count/permutations)  
        high_values$p_value <- replace_na(high_values$p_value,0)

        lcnt = low %>% group_by(sample, reference_phenotype, neighbor_phenotype) %>% summarize(count=n())
        low_values = low_values %>% full_join(lcnt,by=c('sample','reference_phenotype','neighbor_phenotype')) %>% mutate(p_value=count/permutations)
        low_values$p_value <- replace_na(low_values$p_value,0)
        #return(list(low=low_values,high=high_values,annot=annot))
        myna = annot %>% filter(is.na(z_score))
        myna = annot %>% filter(is.na(z_score))
        if (length(myna$sample) > 0) {
           myna$count = NA
           myna$p_value = NA
        }
        #myna$count = NA
        #myna$p_value = NA
        output = rbind(high_values,low_values,myna) %>% arrange(sample,reference_phenotype,neighbor_phenotype)

        return(list(result=output,expected=expected))
    }
)


############ Permutation Test Interaction Comparison ##############
#' Calculate a permutation test result for nearest neighbors to say for each sample to see 
#' if the neighbor distance something seen under the null assumption
#'
#' @param x IrisSpatialFeatures ImageSet object that has had extract nearest neighbors run
#' @param reference Reference phenotype to compare proportions around
#' @param phenotype_A First cell type to compare proportions of
#' @param phenotype_B Second cell type to compare proportions of
#' @param permutations Set to 100 by default
#' @param subset Limit permutations to these types
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
#' interaction_permutation_test(dataset)
#'
#' @importFrom data.table rbindlist
#' @rdname interaction_proportion_comparison_permutation_test
setGeneric("interaction_proportion_comparison_permutation_test", function(x,reference,phenotype_A,phenotype_B,permutations=20,subset=NULL,verbose=FALSE,...)
    standardGeneric("interaction_proportion_comparison_permutation_test"))

#' @rdname interaction_proportion_comparison_permutation_test
#' @aliases interaction_proportion_comparison_permutation_test,ANY,ANY-method
setMethod(
    "interaction_proportion_comparison_permutation_test",
    signature = c("ImageSet"),
    definition = function(x, reference,phenotype_A,phenotype_B,permutations,subset,verbose) {
        obs <- interaction_proportion_data_frame(x) %>% filter(reference_phenotype==reference) %>% filter(neighbor_phenotype %in% c(phenotype_A,phenotype_B))
        # Get our frame for A - B
        obsA = obs %>% filter(neighbor_phenotype == phenotype_A) %>% rename(proportion_A = proportion)
        obsB = obs %>% filter(neighbor_phenotype == phenotype_B) %>% rename(proportion_B = proportion)
        if (length(obsA$sample)==0) { return(NULL) }
        if (length(obsB$sample)==0) { return(NULL) }
        obs = obsA %>% inner_join(obsB,by=c('sample','reference_phenotype')) %>% select(sample,reference_phenotype,proportion_A,proportion_B) %>% mutate(delta = proportion_A-proportion_B)

        if (verbose) { print("finished observed") }
        expected <- lapply(seq(1,permutations),function(i){
            if (verbose) { print(i) }
            vr <- interaction_proportion_data_frame(shuffle_labels(x,subset=subset)) %>% filter(reference_phenotype==reference) %>% filter(neighbor_phenotype %in% c(phenotype_A,phenotype_B))
            expA = vr %>% filter(neighbor_phenotype == phenotype_A) %>% rename(proportion_A = proportion)
            expB = vr %>% filter(neighbor_phenotype == phenotype_B) %>% rename(proportion_B = proportion)
            if (length(expA$sample)==0) { return(NULL) }
            if (length(expB$sample)==0) { return(NULL) }
            vr = expA %>% inner_join(expB,by=c('sample','reference_phenotype')) %>% select(sample,reference_phenotype,proportion_A,proportion_B) %>% mutate(delta = proportion_A-proportion_B)
            vr$iter <- i
            return(vr)
        })
        expected <- rbindlist(expected)
        ci <- expected %>% group_by(sample) %>% summarize(`5%`=quantile(delta,probs=0.05),
                                           `95%`=quantile(delta,probs=0.95),
                                           expected_delta_mean=mean(delta),
                                           expected_delta_sd=sd(delta))

        annot <- obs %>% full_join(ci,by=c('sample')) 

        annot$z_score <- (annot$delta-annot$expected_delta_mean)/(annot$expected_delta_sd/sqrt(permutations))
        annot$permutations = permutations

        # now get the p value
        subset = annot %>% rename(observed = delta) %>% select(sample, observed,z_score,permutations)
        r2 = as.tibble(expected) %>% full_join(subset,by=c('sample'))
        low_values = annot %>% filter(z_score <=0)
        high_values = annot %>% filter(z_score > 0)
        low = r2 %>% filter(z_score <=0) %>% filter(delta <= observed)
        high = r2 %>% filter(z_score >0) %>% filter(delta >= observed)
        hcnt = high %>% group_by(sample) %>% summarize(count=n())
        high_values = high_values %>% full_join(hcnt,by=c('sample')) %>% mutate(p_value=count/permutations)  
        high_values$p_value <- replace_na(high_values$p_value,0)

        lcnt = low %>% group_by(sample) %>% summarize(count=n())
        low_values = low_values %>% full_join(lcnt,by=c('sample')) %>% mutate(p_value=count/permutations)
        low_values$p_value <- replace_na(low_values$p_value,0)
        #return(list(low=low_values,high=high_values,annot=annot))
        myna = annot %>% filter(is.na(z_score))
        myna = annot %>% filter(is.na(z_score))
        if (length(myna$sample) > 0) {
           myna$count = NA
           myna$p_value = NA
        }
        #myna$count = NA
        #myna$p_value = NA
        output = rbind(high_values,low_values,myna) %>% arrange(sample)

        return(list(result=output,expected=expected))
    }
)
