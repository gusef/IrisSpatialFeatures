setGeneric("extract.proximity", function(object, ...) standardGeneric("extract.proximity"))
setMethod("extract.proximity",
          signature = "Iris",
          definition = function(object, radii=c('Entire.Cell.Major.Axis', 'Entire.Cell.Minor.Axis'),
                                uncertainty_margin=1, only_closest=T){
              all_levels <- sort(rownames(object@counts))
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
              total <- matrix(NA,nrow=length(all_levels),ncol=length(all_levels))
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
                          if (t$n == 0){
                              total[to,from] <- 0
                          }else{
                              #get the distances between 'from' and 'to' cells
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





