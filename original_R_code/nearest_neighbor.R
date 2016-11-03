# Install and load required packages
required = c('ggplot2',
             'grid',
             'knitr',
             'raster',
             'RColorBrewer',
             'rmarkdown',
             'spatstat')
for (lib in required)
{
    if (!require(lib, character.only=TRUE))
    {
        install.packages(lib, repos="http://cran.rstudio.com")
        suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
    }
}

MicronsPerPixel <- 0.496

# Columns that must be in every data set
fixedColumns = c('Cell ID', 'Cell X Position', 'Cell Y Position', 'phenotype')



########################################################################################
################## Nearest neighbor algorithm ##########################################
########################################################################################


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

nearestNeighborCoord <- function(coord,classes){
    ppp <- coord$ppp
    res <- lapply(classes,getToNeighbors,classes,ppp)
    means <- t(sapply(res,function(x)x['means',]))
    vars <- t(sapply(res,function(x)x['vars',]))
    nums <- t(sapply(res,function(x)x['nums',]))
    
    rownames(nums) <- rownames(means) <- rownames(vars) <- colnames(means)
    
    #in some cases there are not all samples represented
    means[!is.finite(means)] <- NA
    vars[!is.finite(vars)] <- NA
    nums[!is.finite(nums)] <- NA
    
    return(list(means=means,vars=vars,nums=nums))
}

collapseMatrices <- function(mat,fun){
    #  Make a 3D array from list of matrices
    arr <- array( unlist(mat) , c(nrow(mat[[1]]),nrow(mat[[1]]),length(mat)))
    #  Get mean of third dimension
    collapsed <- fun( arr , dims = 2 ,na.rm = T)
    colnames(collapsed) <- rownames(collapsed) <- colnames(mat[[1]])
    return(collapsed)
}

nearestNeighborSample <- function(sample,classes){
    res <- lapply(sample,nearestNeighborCoord,classes)
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
}


generate_NN <- function(dataset){
    classes <- sort(unique(unlist(lapply(unlist(dataset,recursive = F),function(x)levels(x$ppp$marks)))))
    nn <- lapply(dataset,nearestNeighborSample,classes)
    return(nn)
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

plot_NN <- function(mat,from,to,ttest=TRUE,transposed=F){
    marker_names <- sort(unique(unlist(lapply(mat,function(x)rownames(x$means)))))
    
    #grab the 2 rows that are relevant
    comp <- grep(to,marker_names,fixed = T)
    
    #drop distances to self
    comp <- comp[!marker_names[comp]%in%from]
    x.mean <- extractNNVals(mat,'means',from,comp[1],transposed)
    x.se <- extractNNVals(mat,'SE',from,comp[1],transposed)
    
    if (length(comp)>1){
        y.mean <- extractNNVals(mat,'means',from,comp[2],transposed)
        y.se <- extractNNVals(mat,'SE',from,comp[2],transposed)
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
}


