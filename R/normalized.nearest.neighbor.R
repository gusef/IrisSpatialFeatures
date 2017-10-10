####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

#Let's wait before before including this into the master branch



#####################################################################################################################################
################ nn_comparison_dataframe
#####################################################################################################################################
#
### I think the regular nearest neighbor function might return something very similar
##
##
### And I get an error when trying this piece of code:
### dataset <- IrisSpatialFeatures_data
### extract_nearest_neighbor(dataset)
### nn_comparison_dataframe(dataset,"SOX10+ PDL1-","SOX10+ PDL1+", "CD8+ PD1+")
### nn_comparison_dataframe(dataset,"SOX10+ PDL1-","SOX10+ PDL1+", "CD8+ PD1+",TRUE)
###
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
