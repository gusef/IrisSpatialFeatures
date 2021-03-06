

#' Extract proximity for a range of distances.
#' Pass an array of radii distances (in micrometers) to the function 
#' that would otherwise work like "proximity_sample_data_frame"
#'
#' @param x IrisSpatialFeatures ImageSet object
#' @param radii The array of distances to extract in microns.
#' @param ... Additional arguments
#' @return The proximial events for each cell
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' extract_proximity(dataset)
#'
#' @docType methods
#' @export
#' @rdname multiple_proximity_sample_data_frame
setGeneric("multiple_proximity_sample_data_frame", function(x, ...) standardGeneric("multiple_proximity_sample_data_frame"))

#' @rdname multiple_proximity_sample_data_frame
#' @aliases multiple_proximity_sample_data_frame,ANY,ANY-method
setMethod("multiple_proximity_sample_data_frame",
          signature = "ImageSet",
          definition = function(x, radii=c(10,20,40,100,200,400,1000)){
              all_levels <- x@markers
              v <- lapply(radii,function(distance){
                print(distance)
                x1 <- extract_proximity(x,distance/x@microns_per_pixel)
                df <- proximity_sample_data_frame(x1)
                df$distance <- distance
                return(df)
                })
              df <- do.call(rbind,v)
              return(df)
})

#' Return a data frame of the per-sample proximities
#' like those output by the plot_proximties
#'
#' @param x IrisSpatialFeatures ImageSet object
#' @param ... Additional arguments
#' @return The proximial events for each cell
#' @importFrom reshape2 melt
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' extract_proximity(dataset)
#'
#' @docType methods
#' @export
#' @rdname proximity_sample_data_frame
setGeneric("proximity_sample_data_frame", function(x, ...) standardGeneric("proximity_sample_data_frame"))

#' @rdname extract_proximity
#' @aliases run_proximity,ANY,ANY-method
setMethod("proximity_sample_data_frame",
          signature = "ImageSet",
          definition = function(x){
              if (length(x@proximity)==0){
                  stop(paste('Please run extract.proximity before plotting the interactions.'))
              }
              #return('output')

              # Do the avg_proximities like plot
              avg_proximities <- lapply(x@proximity,function(x)x$avg_proximities)
              df <- melt(avg_proximities)
              colnames(df) <- c('neighbor_phenotype','reference_phenotype','avg_proximity','sample')
              df <- as.tibble(df)
              df2 <- df %>% group_by(reference_phenotype,sample) %>% summarize(cummulative=sum(avg_proximity))
              prox_df <- df %>% inner_join(df2,by=c('reference_phenotype','sample')) %>% mutate(normalized_averaged_fraction=avg_proximity/cummulative)


              # Do the totals in a similar mannor
              totals <- lapply(x@proximity,function(x)x$total)
              df <- melt(totals)
              colnames(df) <- c('neighbor_phenotype','reference_phenotype','neighbor_count','sample')
              totals_df <- as.tibble(df)
              #df <- df %>% inner_join(df2,by=c('reference_phenotype','sample')) %>% mutate(normalized=avg_proximity/cummulative)

              # Get the reference counts
              v1 <- lapply(names(x@samples),function(sample) {
                #print(sample)
                v2 <- lapply(names(x@samples[[sample]]@coordinates),function(frame) {
                  #print(frame)
                  df <- melt(table(x@samples[[sample]]@coordinates[[frame]]@ppp$marks))
                  colnames(df) <- c('reference_phenotype','reference_count')
                  df$frame <- frame
                  return(df)
                })
                v2 <- do.call(rbind,v2)
                v2$sample <- sample
                return(v2)
              })
              cnts <- do.call(rbind,v1)
              # take them to sample level
              sample_counts <- as.tibble(cnts) %>% group_by(reference_phenotype,sample) %>% summarize(reference_count=sum(reference_count))

              # incorperate the sample counts into the totals_df
              totals_df <- totals_df %>% inner_join(sample_counts,by=c('reference_phenotype','sample')) %>% mutate(neighbors_per_reference_cell=neighbor_count/reference_count)
              df2 <- totals_df %>% group_by(sample,reference_phenotype) %>% summarize(neighborhood_total=sum(neighbor_count))
              totals_df <- totals_df %>% inner_join(df2,by=c('reference_phenotype','sample')) %>% mutate(normalized_total_fraction=neighbor_count/neighborhood_total)

              v1 <- prox_df %>% select(sample,reference_phenotype,neighbor_phenotype,avg_proximity,normalized_averaged_fraction)
              v2 <- totals_df %>% select(sample,reference_phenotype,neighbor_phenotype,neighbor_count,reference_count,normalized_total_fraction)
              df <- v1 %>% inner_join(v2,by=c('sample','reference_phenotype','neighbor_phenotype'))
              return(df)
})

#' Run a proximity analysis on all samples.
#' There are two modes this function can be run. In the first mode it uses the major and minor axes for each cell as provided by inform. It then averages
#' half of those two axes and adds an uncertainty margin, which are used to provide an estimate on whether two cells are touching. This mode can be used to approximate the interaction analysis.
#' The second mode uses a user specified distance to count the cells within the proximity of a given cell-type. With increasing distances usually cells fall into the proximity of multiple cells
#' of a given type so the function allows the restriction of only counting the cell only once.
#'
#' @param x IrisSpatialFeatures ImageSet object
#' @param radii The size of the radius in pixels or the default, which are the inForm output columns that indicated the minor and major axis of each cell.
#' @param uncertainty_margin Only for the approximation of the interaction analysis, where it indicates how many pixels further should be search to find a touching cell (deafult: 1).
#' @param only_closest For the proximity analysis a target cell can be in the vincinity of multiple source cells, so the counts are articicially inflated. E.g. a CD8 PD1+ T-cell is within <50 pixels of 30 HRS cells, this cell should only be counted for the closes HRS cell. (default: FALSE)
#' @param ... Additional arguments
#' @return The proximial events for each cell
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' extract_proximity(dataset)
#'
#' @docType methods
#' @export
#' @rdname extract_proximity
setGeneric("extract_proximity", function(x, ...) standardGeneric("extract_proximity"))

#' @rdname extract_proximity
#' @aliases run_proximity,ANY,ANY-method
setMethod("extract_proximity",
          signature = "ImageSet",
          definition = function(x, radii=c('Entire.Cell.Major.Axis', 'Entire.Cell.Minor.Axis'),
                                uncertainty_margin=1, only_closest=FALSE){
              all_levels <- x@markers
              x@proximity <- lapply(x@samples, touches_per_sample, radii,
                                             uncertainty_margin, all_levels, only_closest)
              return(x)
})

setGeneric("touches_per_sample", function(x, ...) standardGeneric("touches_per_sample"))
setMethod("touches_per_sample",
          signature = "Sample",
          definition = function(x, radii, uncertainty_margin, all_levels, only_closest){
              message(paste(x@sample_name,' ... processing...'))
              proximities <- lapply(x@coordinates, touching_events, all_levels,
                                       radii, uncertainty_margin, only_closest)
              avg_proxies <- collapseMatrices(lapply(proximities,function(x)x$avg_touching),rowMeans)
              total <- collapseMatrices(lapply(proximities,function(x)x$total_touching),rowSums)
              return(list(avg_proximities=avg_proxies, total=total))
})


#' @importFrom SpatialTools dist2
setGeneric("touching_events", function(x, ...) standardGeneric("touching_events"))
setMethod("touching_events",
          signature = "Coordinate",
          definition = function(x, all_levels, radii, uncertainty_margin, only_closest){
              gc(verbose = FALSE)
              total <- matrix(0,nrow=length(all_levels),ncol=length(all_levels))
              colnames(total) <- rownames(total) <- all_levels
              for (from in all_levels){
                  #extract the 'from' cells
                  f <- x@ppp[x@ppp$marks==from]
                  fr <- x@raw@data[x@ppp$marks==from,]

                  #count only if there actually cells of the 'from' type present. Otherwise just use the NA's
                  if (f$n > 0){
                      for (to in all_levels){
                          #extract the 'to' cells
                          t <- x@ppp[x@ppp$marks==to]
                          tr <- x@raw@data[x@ppp$marks==to,]

                          #if there are no to cells, just count 0
                          if (t$n >0){                              #get the distances between 'from' and 'to' cells
                              d <- dist2(cbind(t$x, t$y),
                                         cbind(f$x, f$y))
                              #in case we look for connections with same cell-type we need to remove the self distances
                              if (from == to){
                                  diag(d) <- 9999
                              }
                              total[to,from] <- extract_proximity_single(d, fr, tr, radii, uncertainty_margin, only_closest)
                          }
                      }
                  }
              }
              gc()
              #getting the average number of interactions
              counts <- table(x@ppp$marks)[all_levels]
              avg_touching_events <- sweep(total,2,counts,'/')

              return(list(avg_touching=avg_touching_events,
                          total_touching=total))
})


extract_proximity_single <- function(d, fr, tr, radii, uncertainty_margin, only_closest){
    if (class(radii) =='character' && c( c('Entire.Cell.Major.Axis', 'Entire.Cell.Minor.Axis') %in% radii)){
        #print("use uncertainty margin approach")
        #remove average ((minor+major)/2) radius from both cell.
        d <- sweep(d,2,fr[[radii[1]]]/4,FUN = '-')
        d <- sweep(d,2,fr[[radii[2]]]/4,FUN = '-')
        d <- sweep(d,1,tr[[radii[1]]]/4,FUN = '-')
        d <- sweep(d,1,tr[[radii[2]]]/4,FUN = '-')

        #remove a an additional amount for noise
        d <- d - uncertainty_margin

    }else{
        d <- d - radii
    }
    #print(d)
    if (only_closest){
        #when looking at proximity we don't want to double count cells so for each 'to' cell
        #we only count the distance to the closest 'from' cell
        d <- apply(d,1,min)
    }

    #now call all of the distances that are touching
    counts <- sum(d<0)
    return(counts)
}

################################################################
##### Getter functions

#' Get all proximity data for all cell-types in a sample
#'
#' @param x An IrisSpatialFeatures ImageSet object
#' @param ... Additional arguments
#' @return all proximity data for cell-types
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' get_all_proximities(dataset)
#'
#' @docType methods
#' @export
#' @rdname get_all_proximities
setGeneric("get_all_proximities", function(x, ...) standardGeneric("get_all_proximities"))

#' @rdname get_all_proximities
#' @aliases get_all_proximities,ANY,ANY-method
setMethod("get_all_proximities",
          signature = "ImageSet",
          definition = function(x){
              return(x@proximity)
          })

#' Get proximity data for a given cell-type
#'
#' @param x An IrisSpatialFeatures ImageSet object.
#' @param marker Cell type for which the proximity data should be extracted.
#' @param normalize Flag indicating whether the populations should be normalized so that the sum of all is 1 (default: TRUE).
#' @param ... Additional arguments.
#' @return proximities for a specific cell-type
#'
#' @docType methods
#' @rdname get_proximities
setGeneric("get_proximities", function(x, ...) standardGeneric("get_proximities"))

#' @rdname get_proximities
#' @aliases get_proximities,ANY,ANY-method
setMethod("get_proximities",
          signature = "ImageSet",
          definition = function(x,marker,normalize=TRUE){
              if (!marker %in% x@markers){
                  stop(paste('There is no celltype: ',marker))
              }

              int <- lapply(x@proximity,function(x)x$avg_proximities)
              marker_int <- sapply(int,function(x)x[,marker])

              if (normalize){
                  marker_int <- sweep(marker_int,2,colSums(marker_int),'/')
              }
              return(marker_int)
          })

################################################################
##### Interaction summary plotting functions

#' Plot proximity analysis data
#' @param x An IrisSpatialFeatures ImageSet object
#' @param label Cell-type for which the proximit profile is plotted
#' @param ordering Ordering of the samples (Default: NULL)
#' @param normalize Flag, should the populations of different cell-types sum up to one in each sample? (Default: TRUE)
#' @param palette Color palette, by default it uses Spectral from RColorbrewer (Default:NULL)
#' @param celltype_order Ordering of the cell-type. (Default: NULL)
#' @param xlim_fix Space on the right side to show the legend (Default: 13)
#' @param topbar_cols Color of the barplots on top (Default: 'darkgrey'
#' @param ... Additional parameters.
#' @return plot proximity analysis
#'
#' @importFrom graphics axis
#' @importFrom graphics barplot
#' @importFrom graphics layout
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics text
#' @docType methods
#' @rdname plot_proximities
#' @export
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_proximity(dataset,only_closest=TRUE,radii=25)
#' p <- plot_proximities(dataset,"SOX10+ PDL1-",xlim_fix=3)
#'
setGeneric("plot_proximities", function(x, ...) standardGeneric("plot_proximities"))

#' @rdname plot_proximities
#' @aliases plot_proximities,ANY,ANY-method
setMethod("plot_proximities",
          signature = "ImageSet",
          definition = function(x, label, ordering=NULL, normalize=TRUE, palette=NULL,
                                celltype_order=NULL, xlim_fix=13, topbar_cols='darkgrey'){
              if (length(x@proximity)==0){
                  stop(paste('Please run extract.proximity before plotting the interactions.'))
              }

              int <- lapply(x@proximity,function(x)x$avg_proximities)
              # For each of the samples get the avg_proximities
              dat <- sapply(int,function(x)x[,label]) # narrow down the avg_proximities to those referenced by the column 'label'
              count <- get_counts_collapsed(x)[label,] # This appears to be an error.
              labels <- rownames(dat)

              if (normalize){
                  dat <- sweep(dat,2,colSums(dat),'/')
                  ylab <- 'Proportions of cells in proximity'

              }else{
                  ylab <- 'Average number of cells in proximity'
              }

              if (!is.null(ordering)){
                  if(length(ordering)==1){
                      ordering <- order(colSums(dat[grep(ordering,rownames(dat),fixed=TRUE),]),decreasing=TRUE)
                  }
                  dat <- dat[,ordering]
                  count <- count[ordering]
              }

              if (!is.null(celltype_order)){
                  dat <- dat[celltype_order,]
              }

              if (is.null(palette)){
                  palette <- brewer.pal(length(labels),"Spectral")
              }

              #generate the plots
              op <- par(no.readonly = TRUE)
              layout(as.matrix(2:1), widths = c(1), heights = c(0.2,0.8), respect = FALSE)

              par(mar = c(6, 4, 0, 0))
              bp <- barplot(dat,
                            cex.names = 1, # makes x-axis labels small enough to show all
                            col = palette, # colors
                            xlab = 'Sample',
                            ylab = ylab,
                            las=2,
                            xaxt="n",
                            xlim = c(0,ncol(dat)+xlim_fix),
                            width = 1)
              text(cex=1, x=bp+0.8, y=-0.05, colnames(dat), xpd=TRUE, srt=45, pos=2)

              legend("right",
                     legend = labels,
                     fill = palette)

              par(mar = c(0.5, 4, 4, 0))
              p <- barplot(count,
                           xlim = c(0,ncol(dat)+xlim_fix),
                           col = topbar_cols,
                           axisnames = FALSE,
                           axes=FALSE,
                           cex.names = 0.5,
                           main=paste('In proximity to',label))
              axis(side = 2,tick = TRUE, labels = TRUE,line = -1,las=1,cex.axis=0.5)
              par(op)
              return(dat)
          })

