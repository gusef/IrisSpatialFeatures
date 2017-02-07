


#' Extract the distance to each nearest neighbor for each cell-type
#' @export
#' 
#' 
setGeneric("extract.nearest.neighbor", function(object, ...) standardGeneric("extract.nearest.neighbor"))
setMethod("extract.nearest.neighbor",
          signature = "Iris",
          definition = function(object, min_num_cells=10){
    all_levels <- object@markers
    object@nearest_neighbors <- lapply(object@samples,
                                       nearest.neighbor.sample,
                                       all_levels,
                                       min_num_cells)
    return(object)
})


setGeneric("nearest.neighbor.sample", function(object, ...) standardGeneric("nearest.neighbor.sample"))
setMethod("nearest.neighbor.sample",
          signature = "Sample",
          definition = function(object, all_levels, min_num_cells){
    res <- lapply(object@coordinates,nearest.neighbor.coord.raw,all_levels,min_num_cells)

    means <- lapply(res,function(x)x$means)
    vars <- lapply(res,function(x)x$vars)
    nums <- lapply(res,function(x)x$nums)
    
    #collapse the different coordinates
    means <- collapseMatrices(means,rowMeans)
    vars <- collapseMatrices(vars,rowMeans)
    nums <- collapseMatrices(nums,rowSums)
    
    #there is a special case where there is only 1 cell of a type
    #this leads to an NA variance. In this case the means are set to NA as well
    means[is.na(vars)] <- NA
    
    #calculate the standard error on the combined coordinates
    ses <- sqrt(vars)/sqrt(nums)
    return(list(means=means,SE=ses))
})

setGeneric("nearest.neighbor.coord.raw", function(object, ...) standardGeneric("nearest.neighbor.coord.raw"))
setMethod("nearest.neighbor.coord.raw",
          signature = "Coordinate",
          definition = function(object, all_levels, min_num_cells){
              ppp <- object@ppp
              res <- lapply(all_levels,getToNeighbors,all_levels,ppp,min_num_cells)
              means <- t(sapply(res,function(x)x['means',]))
              vars <- t(sapply(res,function(x)x['vars',]))
              nums <- t(sapply(res,function(x)x['nums',]))
              
              rownames(nums) <- rownames(means) <- rownames(vars) <- colnames(means)
              
              #in some cases there are not all samples represented
              means[!is.finite(means)] <- NA
              vars[!is.finite(vars)] <- NA
              nums[!is.finite(nums)] <- NA
              
              return(list(means=means,vars=vars,nums=nums))
})

getNearestNeighbor <- function(from,to,ppp,min_num_cells){
    if (sum(ppp$marks==from) < min_num_cells | sum(ppp$marks==to) < min_num_cells){
        dis <- NA
    }else{
        dis <- nncross(ppp[ppp$marks==from,],ppp[ppp$marks==to,])[,1]
    }
    means <- mean(dis)
    vars <- var(dis)
    nums <- sum(ppp$marks==from)
    return(c(means=means,vars=vars,nums=nums))
}

getToNeighbors <- function(to,classes,ppp,min_num_cells){
    res <- sapply(classes,getNearestNeighbor,to,ppp,min_num_cells)
    return(res)
}

################################################################
##### Nearest Neighbors getters


#' Get the nearest neighbor for each cell-type
#' @export
#' 
#' 
setGeneric("get.all.nearest.neighbors", function(object, ...) standardGeneric("get.all.nearest.neighbors"))
setMethod("get.all.nearest.neighbors",
          signature = "Iris",
          definition = function(object){
              return(object@nearest_neighbors)
          })

#' Get the nearest neighbor for a specified cell-type
#' @export
#' 
#' 
setGeneric("get.nearest.neighbors", function(object, ...) standardGeneric("get.nearest.neighbors"))
setMethod("get.nearest.neighbors",
          signature = "Iris",
          definition = function(object,marker,normalize=T){
              if (!marker %in% object@markers){
                  stop(paste('There is no celltype: ',marker))
              }
              
              nn <- sapply(object@nearest_neighbors,function(x,y)x$means[,y],marker)
              se <- sapply(object@nearest_neighbors,function(x,y)x$SE[,y],marker)
              return(list(mean=nn,SE=se))
          })


#' Plot nearest neighbor barplots for two cell types
#' @export
#' 
#' 
setGeneric("plot.nearest.neighbor", function(object, ...) standardGeneric("plot.nearest.neighbor"))
setMethod("plot.nearest.neighbor",
          signature = "Iris",
          definition = function(object, from, to, ttest=TRUE, transposed=FALSE){
    marker_names <- object@markers
    
    #grab the relevant markers
    comp <- grep(to,marker_names,fixed = T)
    
    if (length(comp)==0 || !from%in%marker_names){
        stop('One or both selected markers are not included in the dataset')    
    }
    
    #drop distances to self
    comp <- comp[!marker_names[comp]%in%from]
    
    x.mean <- extractNNVals(object@nearest_neighbors,'means',from,comp[1],transposed)
    x.se <- extractNNVals(object@nearest_neighbors,'SE',from,comp[1],transposed)
    
    if (length(comp)>1){
        y.mean <- extractNNVals(object@nearest_neighbors,'means',from,comp[2],transposed)
        y.se <- extractNNVals(object@nearest_neighbors,'SE',from,comp[2],transposed)
        current.mean <- (t(cbind(x.mean,y.mean)))
        current.se <- (t(cbind(x.se,y.se)))
    }else{
        current.mean <- t(x.mean)
        current.se <- t(x.se)
    }
    
    #remove NA values
    current.se[is.na(current.mean)] <- 0
    current.mean[is.na(current.mean)] <- 0
    
    max_idx <- which.max(rowSums(current.mean))
    ord <- order(current.mean[max_idx,],decreasing = T)
    current.mean <- current.mean[,ord]
    current.se <- current.se[,ord]
    
    #Fixing the legend for the plot
    if (length(comp)==1){
        ext <- ''
        leg <- marker_names[comp]
        COLS <- c("lightgrey")
    }else{
        ext <-'+/-'
        leg <- c(marker_names[comp[1]],
                 marker_names[comp[2]])
        COLS <- c("lightgrey","black")
    }
    
    bp <- barplot(current.mean, 
                  main=buildLabel(from,to,ext,transposed),
                  xlab="Samples", 
                  ylab="Avg. distance to NN", 
                  col=COLS,
                  legend = leg, 
                  ylim=c(0,max(current.mean)+max(current.se)),
                  las=2,
                  beside=TRUE)
    
    if (length(comp)>1){
        plotSE(bp,current.mean,current.se,1)
        plotSE(bp,current.mean,current.se,2)
    }else{
        bp <- t(bp)
        current.mean <- t(current.mean)
        current.se <- t(current.se)
        plotSE(bp,current.mean,current.se,1)
    }
    
    #paired t test to test for significance
    if (ttest & length(comp)>1){
        pval <- t.test(current.mean[1,],current.mean[2,],paired = T)$p.value
        mtext(paste('Paired t-test:',format(pval,digits=4)),3)
    }else{
        pval <- NA
    }
    return(list(means=current.mean,
                ses=current.se,
                pval=pval))
})

#plot the arrows for the standard error
plotSE <- function(bp,current.mean,current.se,idx){
    arrows(bp[idx,],
           current.mean[idx,]+current.se[idx,],
           bp[idx,], 
           current.mean[idx,], 
           angle=90, 
           code=3, 
           length=0.02,
           col='black')
}

#extract the marker pos/neg values across the samples
extractNNVals <- function(mat,val,from,to,transposed){
    if (transposed){
        dest <- from
        src <- to
    }else{
        dest <- to
        src <- from
    }
    return(sapply(mat,function(x,y,z)x[[val]][z,y],src,dest))
}

#generate the label of the plot
buildLabel <- function(from,to,ext,transposed){
    if (transposed){
        label <- paste('Distance from',to,ext,'to',from)
    }else{
        label <- paste('Distance from',from,'to',to,ext)
    }   
    return(label)
}

#####################################
#ray plots 


#' Plot nearest neighbor ray plots for each samples
#' @export
#' 
#' 
setGeneric("neighbor.ray.plot", function(object, ...) standardGeneric("neighbor.ray.plot"))
setMethod("neighbor.ray.plot",
          signature = "Iris",
          definition = function(object,
                                from_type,
                                to_type,
                                from_col='#EE7600',
                                to_col='#028482',
                                format='.pdf',
                                plot_dir='./',
                                lineColor='#666666',
                                height=7,
                                width=10){
              #generate the mapping directory
              out_dir <- file.path(getwd(),plot_dir)
              if (!file.exists(out_dir)){
                  dir.create(out_dir, showWarnings = FALSE)
              }
              
              #generate ray plots for each sample
              lapply(object@samples, neighbor.ray.plot.sample, from_type, to_type, 
                     from_col, to_col, out_dir, format, lineColor, height, width)
          })

setGeneric("neighbor.ray.plot.sample", function(object, ...) standardGeneric("neighbor.ray.plot.sample"))
setMethod("neighbor.ray.plot.sample",
          signature = "Sample",
          definition = function(object, from_type, to_type, from_col, to_col, out_dir, format, lineColor, height, width){
              lapply(object@coordinates, 
                     neighbor.ray.plot.coord, 
                     object@sample_name,
                     from_type, 
                     to_type, 
                     from_col, 
                     to_col, 
                     out_dir, 
                     format,
                     lineColor,
                     height,
                     width)
          })    

setGeneric("neighbor.ray.plot.coord", function(object, ...) standardGeneric("neighbor.ray.plot.coord"))
setMethod("neighbor.ray.plot.coord",
          signature = "Coordinate",
          definition = function(object, samp_name, from_type, to_type, from_col, to_col, out_dir, format, lineColor, height, width){
    
    #extract the relevant cells 
    if (sum(object@ppp$marks == from_type)>0 && 
        sum(object@ppp$marks == to_type)>0){
        from <- object@ppp[object@ppp$marks == from_type,]
        to <- object@ppp[object@ppp$marks == to_type,]
        #get limits
        overlap <- superimpose(from,to)
        xlim <- max(data.frame(overlap)$x)
        ylim <- max(data.frame(overlap)$y)
        
        #get distances
        dist <- nncross(from, to)
        nearest <- data.frame(
                         cbind(data.frame(from), 
                         dist$dist, 
                         data.frame(to[dist$which,])))
        names(nearest) <- c('from_x','from_y','from_type','dist','to_x','to_y','to_type')
        
        #get the colors right          
        df <- data.frame(overlap)
        df$cols <- as.character(df$marks)
        df$cols[df$cols==from_type] <- from_col    
        df$cols[df$cols==to_type] <- to_col    
                  
                  
        file_stub <- paste0(samp_name,'_',object@coordinate_name)
        if (format == '.pdf'){
            pdf(file = file.path(out_dir,paste0(file_stub,'.pdf')),width=width,height=height)
        }else if(format == '.png'){
            png(file.path(outdir,paste0(file_stub,'.png')),width = 800, height=600)
        }
    
        par(mar=c(4,4,4,1))
        plot(df$x,
             df$y,
             col=as.character(df$cols),
             pch=18,
             ylab='y (pixels)',
             xlab='x (pixels)',
             main=paste(samp_name,'-',object@coordinate_name))      
        segments(nearest$from_x,
                 nearest$from_y,
                 nearest$to_x,
                 nearest$to_y,
                 col=lineColor)    
        legend('bottomleft',
               col = c(from_col,to_col),
               legend = c(from_type,to_type),
               pch=18,cex = 0.8)
        dev.off()
    }
})


