# Install and load required packages
required = c('ggplot2',
             'gplots',
             'grid',
             'knitr',
             'raster',
             'RColorBrewer',
             'rmarkdown',
             'spatstat',
             'Rcpp',
             'SpatialTools',
             'tiff')
for (lib in required)
{
    if (!require(lib, character.only=TRUE))
    {
        install.packages(lib, repos="http://cran.rstudio.com")
        suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
    }
}

sourceCpp('C:/work/projects/imaging/R_functions/interaction.cpp')

########################################################################################
################### Interaction algorithm ##############################################
########################################################################################

extract_interactions <- function(dataset,mem_map){
    #some datasets don't have all cell_types, which makes it tricky to collapse them afterwards
    #to avoid this issue we 
    counts <- extractCounts(dataset)
    all_levels <- colnames(counts)
    interactions <- lapply(1:length(dataset),interactions_per_sample,all_levels,dataset,mem_map)
    names(interactions) <- names(dataset)
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

interactions_per_sample <- function(id,all_levels,dataset,mem_map){
    samp <- dataset[[id]]
    mem <- mem_map[[id]]
    interactions <- lapply(1:length(samp),interaction_events,all_levels,samp,mem)

    ppps <- lapply(interactions,function(x)x$ppp)
    ints <- lapply(interactions,function(x)x$ints)
    names(ints) <- names(ppps) <- names(samp)    
        
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
}

interaction_events <- function(id,all_levels,samp,mem){
    cat('New coordinate... processing...\n')
    
    #extract membrane map and set membranes to -1
    coord <- samp[[id]]
    mem_map <- t(mem[[id]])
    mem_map[mem_map>0] <- -1

    #fill in all of the cells in the membrane map using the cell ID
    ret <- watershed(mem_map,coord)
    
    #update the values
    filled_map <- ret$map
    coord <- ret$coord
    
    #extract the interactions
    interactions <- getNeighbors(filled_map) 
    
    #extract the means and variances
    inter_stats <- extract_interaction_stats(interactions,coord,all_levels)
    
    return(list(stats=inter_stats,
                ppp=coord$ppp,
                ints=interactions)) 
}

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
        
        #caluclate the average interaction measurements
        avg_means <- colMeans(per_cell_summary,na.rm = T)
        avg_vars <- apply(per_cell_summary,2,var,na.rm=T)
    
        
        #caluclate the percent interaction measurements
        per_cell_summary <- per_cell_summary>0
        per_means <- colMeans(per_cell_summary,na.rm = T)
        per_vars <- apply(per_cell_summary,2,var,na.rm=T)
    }
    
    
    return(list(per=list(mean=per_means,var=per_vars),
                avg=list(mean=avg_means,var=avg_vars),
                nums=num_cells))
}

extract_interaction_stats <- function(interactions,coord,all_levels){

    labels <- as.character(coord$ppp$marks)
    
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
}

#generates a mask for each membrane marker
generate_marker_masks <- function(mem_map,ppp){
    lvls <- levels(ppp$marks)
    marker_masks <- lapply(lvls,generate_mask,mem_map,ppp)
    names(marker_masks) <- lvls
    return(marker_masks)
    
}

#watershed algorithm that fills in all cells with 
watershed <- function(mem_map,coord){
    
    #watershed filling in all cells with their cell ID
    padded_map <- rbind(-1,cbind(-1,mem_map,-1),-1)
    
    #need to offset the coordinates because of the padding
    cell_coords <- cbind(1:length(coord$ppp$x),
                         coord$ppp$x,
                         coord$ppp$y)
    
    #run the watershed algorithm and fill up all cells
    ret <- watershedC(padded_map, cell_coords)
    
    #filled in cells
    padded_map <- ret[[1]]
    
    #updated coordinates 
    coord$ppp$x <- ret[[2]][,2]
    coord$ppp$y <- ret[[2]][,3]

    #remove the padding
    padded_map <- padded_map[-c(1,nrow(padded_map)),-c(1,ncol(padded_map))]

    return(list(map=padded_map,coord=coord))
}


getNeighbors <- function(filled_map){
    interactions <- getInteractionsC(filled_map)[[1]]
    interactions <- interactions[-(1:2)]
    #transform the list so the indices correspond to the names
    interactions <- lapply(as.character(1:max(filled_map)),function(x,int)int[[x]],interactions)
    return(interactions)
}


########################################################################################
###################### Old touching algorithm ##########################################
########################################################################################

extract_proximity <- function(d, fr, tr, radii, uncertainty_margin,only_closest){
    if (class(radii) =='character'){
        
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

touching_events <- function(coord,
                            all_levels,
                            radii,
                            uncertainty_margin,
                            only_closest){
    gc()
    cat('New coordinate... processing...\n')
    
    total <- matrix(NA,nrow=length(all_levels),ncol=length(all_levels))
    colnames(total) <- rownames(total) <- all_levels
    
    for (from in all_levels){
        
        #extract the 'from' cells
        f <- coord$ppp[coord$ppp$marks==from]
        fr <- coord$raw$data[coord$ppp$marks==from,]
        
        #count only if there actually cells of the 'from' type present. Otherwise just use the NA's 
        if (f$n > 0){
            for (to in all_levels){
                #extract the 'to' cells
                t <- coord$ppp[coord$ppp$marks==to]
                tr <- coord$raw$data[coord$ppp$marks==to,]
                
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
    counts <- table(coord$ppp$marks)[all_levels]
    avg_touching_events <- sweep(total,2,counts,'/')

    return(list(avg_touching=avg_touching_events,
                total_touching=total))
}

touches_per_sample <- function(samp,
                               radii,
                               uncertainty_margin,
                               all_levels,
                               only_closest){
    touches <- lapply(samp,touching_events,all_levels,radii,uncertainty_margin,only_closest)
    
    avg_touching <- collapseMatrices(lapply(touches,function(x)x$avg_touching),rowMeans)
    total <- collapseMatrices(lapply(touches,function(x)x$total_touching),rowSums)
    
    return(list(avg_touching=avg_touching,
                total=total))
    
}

extract_touches <- function(dataset,
                            radii=c('Entire.Cell.Major.Axis',
                                    'Entire.Cell.Minor.Axis'),
                            uncertainty_margin=1,
                            raw=NA,
                            only_closest=T){
    if (any(!is.na(raw))){
        #combine raw and ppp data
        for (i in 1:length(dataset)){
            for (j in 1: length(dataset[[i]])){
                dataset[[i]][[j]]$raw <- raw[[i]][[j]]
            }
        }
    }    
    
    #some datasets don't have all cell_types, which makes it tricky to collapse them afterwards
    #to avoid this issue we
    counts <- extractCounts(dataset)
    all_levels <- sort(colnames(counts))
    touches <- lapply(dataset,touches_per_sample,radii,uncertainty_margin,all_levels,only_closest)
    return(touches)
}



########################################################################################
############################ Interaction maps ##########################################
########################################################################################


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


interaction_maps <- function(interactions,
                             mem_map,
                             dapi,
                             int_markers,
                             int_marker_cols,
                             silent_markers=c(),
                             silent_col=c(),
                             outline_transparency=0.9,
                             use_dapi=T,
                             outdir='html_reports/interaction_maps',
                             format='.png'){
    #generate the mapping directory
    map_dir <- file.path(getwd(),outdir)

    #get the marker prefix
    marker_prefix <- paste0(c(int_markers,silent_markers), collapse = '__')
    
    #extract all data    
    mem_maps <- unlist(mem_map,recursive = F)    
    dapi_maps <- unlist(dapi,recursive = F)    
    ints <- lapply(interactions,function(x)x$ints)
    ints <- unlist(ints,recursive = F)    
    ppps <- lapply(interactions,function(x)x$ppp)
    ppps <- unlist(ppps,recursive = F) 
    
    #generate a map for each sample
    lapply(names(mem_maps),
           generate_interaction_map,
           mem_maps,
           dapi_maps,
           ints,
           ppps,
           int_markers,
           int_marker_cols,
           silent_markers,
           silent_col,
           outdir,
           marker_prefix,
           map_dir,
           outline_transparency,
           use_dapi,
           format)
    
}


generate_interaction_map <- function(nams,
                                     mem_maps,
                                     dapi_maps,
                                     ints,
                                     ppps,
                                     int_markers,
                                     int_marker_cols,
                                     silent_markers,
                                     silent_col,
                                     outdir,
                                     marker_prefix,
                                     map_dir,
                                     outline_transparency,
                                     use_dapi,
                                     format){
    cat("Working on sample: ",nams,'\n')
    
    #generate the image directory
    samp_dir <- file.path(map_dir)
    #dir.create(samp_dir, showWarnings = FALSE)
    
    all_markers <- c(int_markers,silent_markers)
    
    mem <- t(mem_maps[[nams]])
    mem[mem>0] <- -1
    
    dapi_map <- t(as.matrix(dapi_maps[[nams]]))
    if (use_dapi){
        dapi_map <- dapi_map/max(dapi_map)
    }else{
        dapi_map[dapi_map>0] <- 0
    }
    int <- ints[[nams]]
    ppp <- ppps[[nams]]
    
    if (all(all_markers %in% ppp$marks)){
    
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
            png(file.path(samp_dir,paste0(nams,'_',marker_prefix,'.png')),width=nrow(dapi_map),height=ncol(dapi_map))
            image(dapi_map[,ncol(dapi_map):1],col = cols,breaks=breaks,yaxt='n',xaxt='n')    
            legend('bottomleft',all_markers,col=c(int_marker_cols,silent_col),cex=1.5,pch=18)
            dev.off()
        }else if (format=='.tiff'){
            dapi_map[dapi_map>1] <- 1
            writeTIFF(t(dapi_map),file.path(samp_dir,paste0(nams,'_',marker_prefix,'.tiff')))
        }
    }
}

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
##########################################################################################
############################## Plotting function #########################################
##########################################################################################

#extract the marker pos/neg values across the samples
extractTTVals <- function(mat,val,from,to,transposed,type){
    if (transposed){
        dest <- from
        src <- to
    }else{
        dest <- to
        src <- from
    }
    return(sapply(mat,function(x,y,z)x[[type]][[val]][z,y],src,dest))
}

#generate the label of the plot
buildLabel <- function(from,to,ext,transposed){
    if (transposed){
        label <- paste('Avg. # of ',to,ext,' interacting with',from)
    }else{
        label <- paste('Avg. # of ',from,'interacting with',to,ext)
    }   
    return(label)
}

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

plot_touching <- function(mat,from,to,ttest=TRUE,transposed=F,type='avg'){
    marker_names <- rownames(mat[[1]]$avg$mean)
    
    #grab the 2 rows that are relevant
    comp <- grep(to,marker_names,fixed = T)
    
    x.mean <- extractTTVals(mat,'mean',from,comp[1],transposed,type)
    x.se <- extractTTVals(mat,'SE',from,comp[1],transposed,type)
    
    if (length(comp)>1){
        y.mean <- extractTTVals(mat,'mean',from,comp[2],transposed,type)
        y.se <- extractTTVals(mat,'SE',from,comp[2],transposed,type)
        current.mean <- (t(cbind(x.mean,y.mean)))
        current.se <- (t(cbind(x.se,y.se)))
    }else{
        current.mean <- t(x.mean)
        current.se <- t(x.se)
    }
    
    current.mean[is.na(current.mean)] <- 0
    current.se[is.na(current.se)] <- 0
        
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
                  ylab="Avg. # of interactions", 
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
    return(pval)
}


plot_interactions <- function(label,int,counts,ordering=NULL,normalize=T,palette=NULL,celltype_order=NULL,xlim_fix=13,topbar_cols='darkgrey'){
    dat <- sapply(int,function(x)x[,label])
    count <- counts[label,]

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
    
    
    labels <- rownames(dat)
    
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
}


###########################################################################################################
########### Additional helper functions that make use of the watershed algorithm ##########################
###########################################################################################################

extract_masks_based_on_single_cell <- function(ppp,mem_maps){
    res <- ppp
    #for each samples
    for (i in 1:length(res)){
        #for each coordinate
        for (j in 1:length(res[[i]])){
            if (length(res[[i]][[j]]$ppp$marks)>0){
                coord <- res[[i]][[j]]
                mem_coord <- t(mem_map[[i]][[j]])
                mem_coord[mem_coord>0] <- -1
                
                #run watershed to fix the coordinate issues and get a filled map
                ret <- watershed(mem_coord,coord)
                res[[i]][[j]]$ppp <- ret$coord$ppp
                res[[i]][[j]]$ppp <- res[[i]][[j]]$ppp[res[[i]][[j]]$ppp$x != -1,]
                map <- ret$map
                map[map<0] <- 0
                map[map>0] <- 1
                
            }else{
                map <- mem_map[[i]][[j]]
                map[map!=0] <- 0 
                
            }
            res[[i]][[j]]$map <- map
        }
    }
    return(res)
}



