
setGeneric("extract.interactions", function(object, ...) standardGeneric("extract.interactions"))
setMethod("extract.interactions",
          signature = "Iris",
          definition = function(object){
              object@interactions <- lapply(object@samples,interactions.per.sample,object@markers)
              names(object@interactions) <- names(object@samples)
              return(object)
})

setGeneric("interactions.per.sample", function(object, ...) standardGeneric("interactions.per.sample"))
setMethod("interactions.per.sample",
          signature = "Sample",
          definition = function(object, all_levels){
              message(paste(object@sample_name,' ... processing...'))              
              interactions <- lapply(object@coordinates,interaction.events,all_levels)
            
              areas_with_counts <- sapply(interactions,length) > 0
              if (sum(areas_with_counts)==0){
                  stop('There has to be at least one coordinate that has cells in it to calculate interactions')
              }
              
              interactions <- interactions[areas_with_counts]
              
              ppps <- lapply(interactions,function(x)x$ppp)
              ints <- lapply(interactions,function(x)x$ints)

              #reshape so we get lists of matrices 
              per_means <- lapply(interactions,function(x)x$stats$per$mean)
              per_vars <- lapply(interactions,function(x)x$stats$per$var)
              avg_means <- lapply(interactions,function(x)x$stats$avg$mean)
              avg_vars <- lapply(interactions,function(x)x$stats$avg$var)
              nums <- rowSums(sapply(interactions,function(x)x$stats$nums))
            
              #get the total number of interactions
              total_ints <- lapply(1:length(per_means),function(x,int)sweep(int[[x]]$stats$avg$mean,2,int[[x]]$stats$nums,'*'),interactions)
              
              #collapse the different coordinates
              per_means <- collapseMatrices(per_means,rowMeans)
              per_vars <- collapseMatrices(per_vars,rowMeans)
              avg_means <- collapseMatrices(avg_means,rowMeans)
              avg_vars <- collapseMatrices(avg_vars,rowMeans)
              total_ints <- collapseMatrices(total_ints,rowSums)
            
              #calculate the standard error on the combined coordinates
              per_SE <- sweep(sqrt(per_vars),2,sqrt(nums),'/')
              avg_SE <- sweep(sqrt(avg_vars),2,sqrt(nums),'/')   
            
              return(list(per=list(mean=per_means,
                                   SE=per_SE),
                          avg=list(mean=avg_means,
                                   SE=avg_SE),
                          total=total_ints,
                          ppp=ppps,
                          ints=ints,
                          nums=nums))
})

setGeneric("interaction.events", function(object, ...) standardGeneric("interaction.events"))
setMethod("interaction.events",
          signature = "Coordinate",
          definition = function(object, all_levels){
              #extract membrane map and set membranes to -1
              if (is.null(object@raw@mem_seg_map)){
                  stop('The interaction analysis can only be run on datasets that include the membrane maps. Try the proximity analysis instead.')
              }
              
              if (length(object@ppp$x)==0){
                  return(list())
              }

              #fill in all of the cells in the membrane map using the cell ID
              ret <- watershed(object)
            
              #update the values
              filled_map <- ret$map
              object <- ret$object
            
              #extract the interactions
              interactions <- getNeighbors(filled_map) 
            
              #extract the means and variances
              inter_stats <- extract.interaction.stats(object,interactions,all_levels)
            
              return(list(stats=inter_stats,
                          ppp=object@ppp,
                          ints=interactions)) 
})


setGeneric("watershed", function(object, ...) standardGeneric("watershed"))
setMethod("watershed",
          signature = "Coordinate",
          definition = function(object){
              mem_map <- t(object@raw@mem_seg_map)
              mem_map[mem_map>0] <- -1
              #watershed filling in all cells with their cell ID
              padded_map <- rbind(-1,cbind(-1,mem_map,-1),-1)
    
              #need to offset the coordinates because of the padding
              cell_coords <- cbind(1:length(object@ppp$x),
                                   object@ppp$x,
                                   object@ppp$y)
    
              #run the watershed algorithm and fill up all cells
              ret <- watershedC(padded_map, cell_coords)
                
              #filled in cells
              padded_map <- ret[[1]]
               
              #updated coordinates 
              object@ppp$x <- ret[[2]][,2]
              object@ppp$y <- ret[[2]][,3]
                
              #remove the padding
              padded_map <- padded_map[-c(1,nrow(padded_map)),-c(1,ncol(padded_map))]
                
              return(list(map=padded_map,object=object))
})


get_single_int <- function(lvl, int, labels, all_levels){
    #generate a per cell summary
    ints <- lapply(int[labels==lvl],function(i,lev)factor(i,levels=lev),all_levels)
    
    if (length(ints)==0){
        num_cells <- 0
        per_means <- rep(NA,length(all_levels))
        names(per_means) <- all_levels
        avg_means <- per_vars <- avg_vars <- per_means
    }else{
        per_cell_summary <- t(sapply(ints,table))
        
        num_cells <- nrow(per_cell_summary)
        
        #calculate the average interaction measurements
        avg_means <- colMeans(per_cell_summary,na.rm = T)
        avg_vars <- apply(per_cell_summary,2,var,na.rm=T)
        
        
        #calculate the percent interaction measurements
        per_cell_summary <- per_cell_summary>0
        per_means <- colMeans(per_cell_summary,na.rm = T)
        per_vars <- apply(per_cell_summary,2,var,na.rm=T)
    }
    
    
    return(list(per=list(mean=per_means,var=per_vars),
                avg=list(mean=avg_means,var=avg_vars),
                nums=num_cells))
}

setGeneric("extract.interaction.stats", function(object, ...) standardGeneric("extract.interaction.stats"))
setMethod("extract.interaction.stats",
          signature = "Coordinate",
          definition = function(object,interactions,all_levels){
        
          labels <- as.character(object@ppp$marks)
            
          #translate the coordinates in the lists to labels
          int <- lapply(interactions,function(x,lab)lab[x],labels)
            
          #get all the stats
          stats <- lapply(all_levels,get_single_int,int,labels,all_levels)
            
          #reshape the data so we get 4 matrices and 1 vector
          per_means <- sapply(stats,function(x)x$per$mean)
          per_vars <- sapply(stats,function(x)x$per$var)
          avg_means <- sapply(stats,function(x)x$avg$mean)
          avg_vars <- sapply(stats,function(x)x$avg$var)
          nums <- sapply(stats,function(x)x$nums)
            
          #make sure we get the right column names
          colnames(per_means) <- colnames(per_vars) <-  all_levels
          colnames(avg_means) <- colnames(avg_vars) <-  all_levels
          names(nums) <- all_levels
            
          return(list(per=list(mean=per_means,var=per_vars),
                      avg=list(mean=avg_means,var=avg_vars),
                      nums=nums))
})


getNeighbors <- function(filled_map){
    interactions <- getInteractionsC(filled_map)[[1]]
    interactions <- interactions[-(1:2)]
    #transform the list so the indices correspond to the names
    interactions <- lapply(as.character(1:max(filled_map)),function(x,int)int[[x]],interactions)
    return(interactions)
}

collapseMatrices <- function(mat,fun){
    #  Make a 3D array from list of matrices
    arr <- array( unlist(mat) , c(nrow(mat[[1]]),nrow(mat[[1]]),length(mat)))
    #  Get mean of third dimension
    collapsed <- fun( arr , dims = 2 ,na.rm = T)
    colnames(collapsed) <- rownames(collapsed) <- colnames(mat[[1]])
    return(collapsed)
}

################################################################
##### Interaction getters

setGeneric("get.all.interactions", function(object, ...) standardGeneric("get.all.interactions"))
setMethod("get.all.interactions",
          signature = "Iris",
          definition = function(object){
              return(object@interactions)
})


setGeneric("get.interactions", function(object, ...) standardGeneric("get.interactions"))
setMethod("get.interactions",
          signature = "Iris",
          definition = function(object,marker,normalize=T){
              if (!marker %in% object@markers){
                  stop(paste('There is no celltype: ',marker))
              }

              int <- lapply(object@interactions,function(x)x$avg$mean)
              marker_int <- sapply(int,function(x)x[,marker])

              if (normalize){
                  marker_int <- sweep(marker_int,2,colSums(marker_int),'/')
              }
              return(marker_int)
})

################################################################
##### Interaction summary plotting functions

setGeneric("plot.interactions", function(object, ...) standardGeneric("plot.interactions"))
setMethod("plot.interactions",
          signature = "Iris",
          definition = function(object, label, ordering=NULL, normalize=T, palette=NULL,
                                celltype_order=NULL, xlim_fix=13, topbar_cols='darkgrey'){
          if (length(object@interactions)==0){
              stop(paste('Please run extract.interactions before plotting the interactions.'))
          }
              
          int <- lapply(object@interactions,function(x)x$avg$mean)
          dat <- sapply(int,function(x)x[,label])
          count <- get.counts.collapsed(object)[label,]
          labels <- rownames(dat)
          
          if (normalize){
              dat <- sweep(dat,2,colSums(dat),'/')
              ylab <- 'Proportions of interactions'
                
          }else{
              ylab <- 'Average interactions'
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
                       main=paste('Interactions with',label)) 
          axis(side = 2,tick = T, labels = T,line = -1,las=1,cex.axis=0.5)
          par(op)
          return(dat)
})

################################################################
##### Interaction maps

setGeneric("interaction.maps", function(object, ...) standardGeneric("interaction.maps"))
setMethod("interaction.maps",
          signature = "Iris",
          definition = function(object,
                                int_markers,
                                int_marker_cols,
                                silent_markers=c(),
                                silent_col=c(),
                                outline_transparency=0.9,
                                use_dapi=F,
                                outdir='interaction_maps',
                                format='.png'){
    #generate the mapping directory
    map_dir <- file.path(getwd(),outdir)
    if (!file.exists(map_dir)){
        dir.create(map_dir, showWarnings = FALSE)
    }

    #generate a map for each sample
    lapply(object@samples, interaction.map.sample, object@interactions, int_markers, int_marker_cols,
           silent_markers, silent_col, map_dir, outline_transparency, use_dapi, format)
    return('Done!')
    
})

setGeneric("interaction.map.sample", function(object, ...) standardGeneric("interaction.map.sample"))
setMethod("interaction.map.sample",
          signature = "Sample",
          definition = function(object, interactions, int_markers, int_marker_cols, silent_markers,
                                silent_col, map_dir, outline_transparency, use_dapi, format){
              
              message("Working on sample: ",object@sample_name)
              lapply(object@coordinates, generate.interaction.map, object@sample_name, 
                     interactions[[object@sample_name]], int_markers, int_marker_cols, 
                     silent_markers, silent_col, map_dir, outline_transparency, use_dapi, format)
})    
    
setGeneric("generate.interaction.map", function(object, ...) standardGeneric("generate.interaction.map"))
setMethod("generate.interaction.map",
              signature = "Coordinate",
              definition = function(object, samp_name, interactions, int_markers, int_marker_cols, silent_markers,
                                    silent_col, map_dir, outline_transparency, use_dapi, format){

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

fill_in_maps <- function(marker,markers,marker_masks,mem,int,ppp){
    mask <- marker_masks[[marker]]
    others <- markers[markers != marker]
    labels <- as.character(ppp$marks)
    
    #translate the interactions into actual labels
    inter <- lapply(int,function(x,lab)lab[x],labels)
    
    #look only at the cells from the current marker
    indicator <- labels == marker
    
    #focus only on the cells with the current marker
    inter <- inter[indicator]
    coords <- ppp[indicator,]
    
    #coordinates to fill
    inter <- sapply(inter,function(x,others)sum(x%in%others)>0,others)
    coords <- data.frame(coords)[inter,1:2]
    coords <- as.matrix(coords)
    
    if (nrow(coords)>0){
        #current mask
        mask <- rbind(-1,cbind(-1,mask,-1),-1)
        padded <- rbind(-1,cbind(-1,mem,-1),-1)
        #c code
        mask <- fillMaskC(mask, padded, coords)
        #remove the padding
        mask <- mask[-c(1,nrow(mask)),-c(1,ncol(mask))]
    }
    return(mask)
}    

#generates a mask for a single marker
generate_mask <- function(lvl,mem,ppp){
    
    #add padding to the membrane map so the watershed cannot go outside of the boundaries
    padded_map <- rbind(-1,cbind(-1,mem,-1),-1)
    
    #extract the coordinates for the current marker
    marker_coords <- ppp[ppp$marks==lvl,]
    #adust the cell coordinates to take the padding into account
    cell_coords <- cbind(marker_coords$x,marker_coords$y)
    
    #make a blank mask
    marker_map <- padded_map
    marker_map[marker_map!=0] <- 0
    marker_map <- generate_maskC(marker_map, padded_map, cell_coords)
    
    #remove the padding
    marker_map <- marker_map[-c(1,nrow(marker_map)),-c(1,ncol(marker_map))]
    
    return(marker_map)
    
}

