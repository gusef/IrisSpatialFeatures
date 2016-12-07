



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

#####################################
#ray plots 

setGeneric("neighbor.ray.plot", function(object, ...) standardGeneric("neighbor.ray.plot"))
setMethod("neighbor.ray.plot",
          signature = "Iris",
          definition = function(object,
                                from_type,
                                to_type,
                                from_col='#EE7600',
                                to_col='#028482',
                                format='.pdf',
                                plot_dir='plots'){
              #generate the mapping directory
              out_dir <- file.path(getwd(),plot_dir)
              if (!file.exists(out_dir)){
                  dir.create(out_dir, showWarnings = FALSE)
              }
              
              #generate ray plots for each sample
              lapply(object@samples, neighbor.ray.plot.sample, from_type, to_type, from_col, to_col, out_dir, format)
              return('Done!')
              
          })

setGeneric("neighbor.ray.plot.sample", function(object, ...) standardGeneric("neighbor.ray.plot.sample"))
setMethod("neighbor.ray.plot.sample",
          signature = "Sample",
          definition = function(object, from_type, to_type, from_col, to_col, out_dir, format){
              lapply(object@coordinates, neighbor.ray.plot.coord, from_type, to_type, from_col, to_col, out_dir, format)
          })    

setGeneric("neighbor.ray.plot.coord", function(object, ...) standardGeneric("neighbor.ray.plot.coord"))
setMethod("neighbor.ray.plot.coord",
          signature = "Coordinate",
          definition = function(object, from_type, to_type, from_col, to_col, out_dir, format){
              
              nams <- paste(samp_name,object@coordinate_name,sep='_')
              #extract all data    
              int <- interactions$ints[[object@coordinate_name]]
              ppp <- interactions$ppp[[object@coordinate_name]]                  
              
              #get the marker prefix
              marker_prefix <- paste0(c(int_markers,silent_markers), collapse = '__')
              
              #extract membrane map and set membranes to -1
              if (is.null(object@raw@mem_seg_map)){
                  stop('The interaction maps can only be created on datasets that include the membrane maps.')
              }
              mem <- t(object@raw@mem_seg_map)
              mem[mem>0] <- -1
              
              if (use_dapi){
                  #extract membrane map and set membranes to -1
                  if (is.null(object@raw@dapi_map)){
                      stop('No DAPI map available, please set dapi_map flag to FALSE')
                  } 
                  dapi_map <- t(as.matrix(object@raw@dapi_map))
                  dapi_map <- dapi_map/max(dapi_map)
              }else{
                  dapi_map <- mem
                  dapi_map[dapi_map!=0] <- 0
              }
              
              if (all(int_markers %in% ppp$marks)){
                  
                  #generate the masks
                  #first the interaction ones
                  int_marker_masks <- lapply(int_markers,generate_mask,mem,ppp)
                  names(int_marker_masks) <- int_markers
                  
                  #then the ones we don't fill out
                  sil_marker_masks <- lapply(silent_markers,generate_mask,mem,ppp)
                  names(sil_marker_masks) <- silent_markers
                  
                  #fill in the marker masks that are relevant for interactions
                  int_marker_masks <- lapply(int_markers,fill_in_maps,int_markers,int_marker_masks,mem,int,ppp)
                  
                  #make the outlines transparent to increase the contrast
                  int_marker_cols2 <- unlist(lapply(int_marker_cols,function(x)c(rgb2hex(col2rgb(x)[,1]*outline_transparency),x)))
                  silent_col <- unlist(lapply(silent_col,function(x)rgb2hex(col2rgb(x)[,1]*outline_transparency)))
                  
                  #color panel
                  cols <- c(colorpanel(10,'black','blue'),silent_col,int_marker_cols2)
                  breaks <- c(seq(0,1,0.1),(1:(length(silent_col)+length(int_marker_cols2)))+1.5)
                  
                  #add up all the outlines
                  if (length(sil_marker_masks)>0){
                      for (i in 1:length(sil_marker_masks)){
                          dapi_map[sil_marker_masks[[i]]==1] <- i+1;
                      }
                  }
                  
                  col_count <- length(sil_marker_masks)
                  #add up all the masks that have interactions
                  for (i in 1:length(int_marker_masks)){
                      dapi_map[int_marker_masks[[i]]==1] <- col_count + i+1;
                      dapi_map[int_marker_masks[[i]]==2] <- col_count + i+2;
                      col_count <- col_count+1
                  }
                  if (format=='.png'){
                      png(file.path(map_dir,paste0(nams,'_',marker_prefix,'.png')),width=nrow(dapi_map),height=ncol(dapi_map))
                      image(dapi_map[,ncol(dapi_map):1],col = cols,breaks=breaks,yaxt='n',xaxt='n')    
                      legend('bottomleft',c(int_markers,silent_markers),col=c(int_marker_cols,silent_col),cex=1.5,pch=18)
                      dev.off()
                  }else if (format=='.tiff'){
                      dapi_map[dapi_map>1] <- 1
                      writeTIFF(t(dapi_map),file.path(map_dir,paste0(nams,'_',marker_prefix,'.tiff')))
                  }
              }
          })



make_NN_plot <- function(pheno1,pheno2,labels,main='',palette=c('#EE7600','#3182bd')){
    #get distances
    distance <- nncross(pheno1,pheno2)
    nn <- cbind(data.frame(pheno1), 
                distance$dist, 
                data.frame(pheno2[distance$which,]))
    colnames(nn) <- c('p1x','p1y','marks_p1','dist','p2x','p2y','marks_p2')
    nn <- data.frame(nn)
    
    #get limits
    superset <- superimpose(pheno1,pheno2)
    xlim <- max(data.frame(superset)$x)
    ylim <- max(data.frame(superset)$y)
    
    #generate the plot
    plot <- nnPlot(nn,
                   pheno2,
                   pnames=labels,
                   pcols=c(palette[1],palette[2]),
                   xlim,
                   ylim,
                   main=main)
}


nnPlot <- function (nn, pheno2, pnames, pcols, xlim, ylim,main,lineColor='gray40') {
    title = paste(main,'NN from',pnames[1],'to',pnames[2])
    p = ggplot(data=data.frame(x=0, y=0), aes(x=x, y=y)) # Fake d.f needed to get background to draw...
    p = p + labs(x='Cell X Position', y='Cell Y Position', title=title)
    addScalesAndBackground(p,xlim,ylim)
    p = p + geom_segment(data = data.frame(nn),
                         aes(x=`p1x`, y=`p1y`, xend=`p2x`, yend=`p2y`), color=lineColor)
    p = p + geom_point(data=data.frame(nn), 
                       aes(x=`p1x`, y=`p1y`),
                       color=pcols[1])
    p = p + geom_point(data=data.frame(pheno2), 
                       aes(`x`, `y`), 
                       color=pcols[2])
    p = p + scale_y_reverse()
    p
}


