



setGeneric("extract.nearest.neighbor", function(object, ...) standardGeneric("extract.nearest.neighbor"))
setMethod("extract.nearest.neighbor",
          signature = "Iris",
          definition = function(object){
    all_levels <- object@markers
    object@nearest_neighbors <- lapply(object@samples,nearest.neighbor.sample,all_levels)
    return(object)
})


setGeneric("nearest.neighbor.sample", function(object, ...) standardGeneric("nearest.neighbor.sample"))
setMethod("nearest.neighbor.sample",
          signature = "Sample",
          definition = function(object, all_levels){
    res <- lapply(object@coordinates,nearest.neighbor.coord,all_levels)
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

setGeneric("nearest.neighbor.coord", function(object, ...) standardGeneric("nearest.neighbor.coord"))
setMethod("nearest.neighbor.coord",
          signature = "Coordinate",
          definition = function(object, all_levels){
              ppp <- object@ppp
              res <- lapply(all_levels,getToNeighbors,all_levels,ppp)
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

getNearestNeighbor <- function(from,to,ppp){
    dis <- nncross(ppp[ppp$marks==from,],ppp[ppp$marks==to,])[,1]
    means <- mean(dis)
    vars <- var(dis)
    nums <- length(dis)
    return(c(means=means,vars=vars,nums=nums))
}

getToNeighbors <- function(to,classes,ppp){
    res <- sapply(classes,getNearestNeighbor,to,ppp)
    return(res)
}

################################################################
##### Nearest Neighbors getters

setGeneric("get.all.nearest.neighbors", function(object, ...) standardGeneric("get.all.nearest.neighbors"))
setMethod("get.all.nearest.neighbors",
          signature = "Iris",
          definition = function(object){
              return(object@nearest_neighbors)
          })

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


