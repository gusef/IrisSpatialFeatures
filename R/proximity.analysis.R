
#' Run a proximity analysis on all samples
#' @export
#' 
#' 
setGeneric("extract.proximity", function(object, ...) standardGeneric("extract.proximity"))
setMethod("extract.proximity",
          signature = "Iris",
          definition = function(object, radii=c('Entire.Cell.Major.Axis', 'Entire.Cell.Minor.Axis'),
                                uncertainty_margin=1, only_closest=F){
              all_levels <- object@markers
              object@proximity <- lapply(object@samples, touches.per.sample, radii, 
                                             uncertainty_margin, all_levels, only_closest)
              return(object)
})
              
setGeneric("touches.per.sample", function(object, ...) standardGeneric("touches.per.sample"))
setMethod("touches.per.sample",
          signature = "Sample",
          definition = function(object, radii, uncertainty_margin, all_levels, only_closest){
              message(paste(object@sample_name,' ... processing...'))           
              proximities <- lapply(object@coordinates, touching.events, all_levels, 
                                       radii, uncertainty_margin, only_closest)
              avg_proxies <- collapseMatrices(lapply(proximities,function(x)x$avg_touching),rowMeans)
              total <- collapseMatrices(lapply(proximities,function(x)x$total_touching),rowSums)
              return(list(avg_proximities=avg_proxies, total=total))
})
              
setGeneric("touching.events", function(object, ...) standardGeneric("touching.events"))
setMethod("touching.events",
          signature = "Coordinate",
          definition = function(object, all_levels, radii, uncertainty_margin, only_closest){
              gc(verbose = F)
              total <- matrix(0,nrow=length(all_levels),ncol=length(all_levels))
              colnames(total) <- rownames(total) <- all_levels
              for (from in all_levels){
                  #extract the 'from' cells
                  f <- object@ppp[object@ppp$marks==from]
                  fr <- object@raw@data[object@ppp$marks==from,]
            
                  #count only if there actually cells of the 'from' type present. Otherwise just use the NA's 
                  if (f$n > 0){
                      for (to in all_levels){
                          #extract the 'to' cells
                          t <- object@ppp[object@ppp$marks==to]
                          tr <- object@raw@data[object@ppp$marks==to,]
                    
                          #if there are no to cells, just count 0
                          if (t$n >0){                              #get the distances between 'from' and 'to' cells
                              d <- dist2(cbind(t$x, t$y),
                                         cbind(f$x, f$y))
                              #in case we look for connections with same cell-type we need to remove the self distances
                              if (from == to){
                                  diag(d) <- 9999
                              }
                              total[to,from] <- extract_proximity(d, fr, tr, radii, uncertainty_margin, only_closest)
                          }
                      }
                  }
              }
              gc()
              #getting the average number of interactions
              counts <- table(object@ppp$marks)[all_levels]
              avg_touching_events <- sweep(total,2,counts,'/')
            
              return(list(avg_touching=avg_touching_events,
                          total_touching=total))
})
             

extract_proximity <- function(d, fr, tr, radii, uncertainty_margin, only_closest){
    if (class(radii) =='character' && c( c('Entire.Cell.Major.Axis', 'Entire.Cell.Minor.Axis') %in% radii)){
        
        #remove radius from first cell
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
        #when looking at proximity we don't want to double count cells
        d <- apply(d,1,min)
    }
    
    #now call all of the distances that are touching
    counts <- sum(d<0)
    
    return(counts)
}

################################################################
##### Getter functions

#' Get all proximity data for all cell-types in a sample
#' @export
#' 
#' 
setGeneric("get.all.proximities", function(object, ...) standardGeneric("get.all.proximities"))
setMethod("get.all.proximities",
          signature = "Iris",
          definition = function(object){
              return(object@proximity)
          })

#' Get proximity data for a given cell-type
#' @export
#' 
#' 
setGeneric("get.proximities", function(object, ...) standardGeneric("get.proximities"))
setMethod("get.proximities",
          signature = "Iris",
          definition = function(object,marker,normalize=T){
              if (!marker %in% object@markers){
                  stop(paste('There is no celltype: ',marker))
              }
              
              int <- lapply(object@proximity,function(x)x$avg_proximities)
              marker_int <- sapply(int,function(x)x[,marker])
              
              if (normalize){
                  marker_int <- sweep(marker_int,2,colSums(marker_int),'/')
              }
              return(marker_int)
          })

################################################################
##### Interaction summary plotting functions

#' Plot proximity analysis data
#' @export
#' 
#' 
setGeneric("plot.proximities", function(object, ...) standardGeneric("plot.proximities"))
setMethod("plot.proximities",
          signature = "Iris",
          definition = function(object, label, ordering=NULL, normalize=T, palette=NULL,
                                celltype_order=NULL, xlim_fix=13, topbar_cols='darkgrey'){
              if (length(object@proximity)==0){
                  stop(paste('Please run extract.proximity before plotting the interactions.'))
              }
              
              int <- lapply(object@proximity,function(x)x$avg_proximities)
              dat <- sapply(int,function(x)x[,label])
              count <- get.counts.collapsed(object)[label,]
              labels <- rownames(dat)
              
              if (normalize){
                  dat <- sweep(dat,2,colSums(dat),'/')
                  ylab <- 'Proportions of cells in proximity'
                  
              }else{
                  ylab <- 'Average number of cells in proximity'
              }
              
              if (!is.null(ordering)){
                  if(length(ordering)==1){
                      ordering <- order(colSums(dat[grep(ordering,rownames(dat),fixed=T),]),decreasing=T)
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
                           axes=F,
                           cex.names = 0.5,
                           main=paste('In proximity to',label)) 
              axis(side = 2,tick = T, labels = T,line = -1,las=1,cex.axis=0.5)
              par(op)
              return(dat)
          })

