
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
#' extract_proximity(new("ImageSet"))
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
#' get_all_proximities(new("ImageSet"))
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
              dat <- sapply(int,function(x)x[,label])
              count <- get_counts_collapsed(x)[label,]
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

