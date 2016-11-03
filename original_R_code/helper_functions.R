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
##################### Get all count data################################################
########################################################################################
extractCounts <- function(dataset){
    images <- unlist(dataset,recursive = F)
    counts <- extractCountsSample(images)
    counts <- counts[,sort(colnames(counts))]
    return(counts)
}

extractCountsF <- function(x,counter){
    counter[match(names(x),names(counter))] <- x
    return(counter)
}

extractCountsSample <- function(sample){
    counts <- lapply(sample,function(x)table(x$ppp$marks))
    nams <- unique(unlist(lapply(counts,names)))
    counter <- rep(0,length(nams))
    names(counter) <- nams
    counts <- t(sapply(counts,extractCountsF,counter))
    return(counts)
}


extractCountsNonCollapsed <- function(dataset){
    counts <- lapply(dataset,extractCountsSample)
    nams <- unique(unlist(lapply(counts,colnames)))
    
    standardize <- function(x,nams){
        y <- matrix(0,nrow=nrow(x),ncol=length(nams))
        colnames(y) <- nams
        rownames(y) <- rownames(x)
        y[,colnames(x)] <- x
        return(y)
    }    
    counts <- lapply(counts,standardize,nams)
    return(counts)
}


extractCountsCollapsed <- function(dataset){
    counts <- lapply(dataset,extractCountsSample)
    nams <- unique(unlist(lapply(counts,colnames)))
    for (i in 1:length(counts)){
        counts[[i]] <- counts[[i]][,match(nams,colnames(counts[[i]]))]
    }
    combined <- sapply(counts,colSums,na.rm=T)
    
    if (class(combined) != 'matrix'){
        nams <- sort(unique(unlist(lapply(combined,names))))
        counter <- rep(0,length(nams))
        names(counter) <- nams
        counts <- sapply(combined,extractCountsF,counter)
    }else{
        counts <- combined
        rownames(counts) <- nams
    }
    counts <- counts[order(rownames(counts)),]
    return(counts)
}


########################################################################################
######## Report and plotting functions #################################################
########################################################################################

#generates a summary report of all provided files in the imgdata ppp object list
runReport <- function(images,base_dir='html_reports',sample_report='sample_report.Rmd',output_dir=getwd()){

    
    #add a slash to the base directory in case there isn't already one
    if (substr(base_dir,nchar(base_dir),nchar(base_dir)) != '/'){
        base_dir <- paste0(base_dir,'/')
    }
    
    base_dir <- file.path(output_dir,base_dir)
    
    #if base directory doesn't exist create it
    if (!file.exists(sub('/$','',base_dir))){
        dir.create(base_dir)
    }
    
    #if htmls directory doesn't exist create it
    if (!file.exists(file.path(base_dir,'htmls'))){
        dir.create(file.path(base_dir,'htmls'))
    }
    
    #make a html report
    overview_file <- file.path(base_dir,'index.html')
    unlink(overview_file)
    file_con <- file(overview_file,'w')
    write('<h1>Dataset overview</h1>',append=T, file_con)
    
    for (idx in 1:length(images)){
        #get data, title and the filename for the subreport
        imgdata <- images[[idx]]$ppp
        main <- sub('\\.',' - Coordinate: ',names(images)[idx])
        fname <- sub(',','_',names(images)[idx])
        html_link <- paste0('htmls/report_',fname,'.html')
        img_file <- paste0('htmls/',fname,'.png')
        
        #write into report
        write(paste0('<h2>Sample ',main,'</h2>'),append=T, file_con)
        write(paste0('<img src="',img_file,'">'),append=T, file_con)
        write(paste0('<a href="',html_link,'">Full report</a>'),append=T, file_con)
        
        #plot overview figure
        png(file.path(base_dir,img_file),width = 800, height=600)
        p <- plotOverview(imgdata)
        print(p)
        dev.off()
        
        render(sample_report, 
               'html_document', 
               file.path(base_dir,html_link))
    }
    close(file_con)
}


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


plotOverview <- function(imgdata){
    xlim <- max(data.frame(imgdata)$x)
    ylim <- max(data.frame(imgdata)$y)
    
    p <- ggplot(data.frame(imgdata), aes(x=`x`, y=`y`, color=marks))
    p <- p + scale_color_brewer(palette='Spectral')
    p <- addScalesAndBackground(p,xlim,ylim)
    p <- p + geom_point()
    p + labs(title='Locations of all cells')
    p
}


make_highlight_plot <- function(lab1,lab2,ppp,main='',palette=c('#3157bd','#EE7600')){
    #get limits
    xlim <- max(data.frame(ppp)$x)
    ylim <- max(data.frame(ppp)$y)
    
    title <- paste('Overview',lab1,'/',lab2)
    
    p = ggplot(data=data.frame(x=0, y=0), aes(x=x, y=y))
    p = p + labs(x='Cell X Position', y='Cell Y Position', title=title)
    
    p = p + geom_point(data=data.frame(ppp), 
                       aes(`x`, `y`),
                       color='lightgrey')
    
    p = p + geom_point(data=data.frame(ppp)[ppp$marks==lab1,], 
                       aes(`x`, `y`),
                       color=palette[1])
    
    p = p + geom_point(data=data.frame(ppp)[ppp$marks==lab2,], 
                       aes(`x`, `y`), 
                       color=palette[2])
    p = p + theme_bw()
    p = p + scale_y_reverse()
    p
    
}


# Add scales, scale line and background image to a ggplot object
# background, xlim, ylim are taken from the environment <cough> <cough>
addScalesAndBackground = function(p,xlim,ylim){
    
    # Add scales at the image limits. Reverse the y scale to match the image
    p = p + scale_x_continuous(limits=c(0, xlim)) + scale_y_reverse(limits=c(ylim, 0))
    
    # Force square aspect ratio
    p = p + coord_fixed()
    
    # Add background image if we have one
    if (exists('background'))
    {
        p = p + annotation_raster(background, xmin=0, xmax=xlim, ymin=0, ymax=-ylim)
    }
    
    # Add a 200-micron line segment for scale reference
    p = p + geom_segment(aes(x=xlim-50-200, xend=xlim-50, y=ylim-100, yend=ylim-100), color='black', size=1)
    p = p + geom_text(aes(x=xlim-50-200/2, y=ylim-90, label=paste(200, '~mu*m')), 
                      size=3, hjust=0.5, vjust=1, color='black', parse=TRUE)
    p
}

theme_reports = function(){
    theme_update(
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.minor = element_line(colour = "grey99", size = 0.25),
        strip.background = element_rect(fill = "grey90", color="grey50"))
}

# Multiple plot function, to print a grid of ggplot objects
# From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(plotlist, cols=1, layout=NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- plots
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


################################################################################################
########################### Extract data within a mask #########################################
################################################################################################

filter_coord <- function(coord,mask){
    dat <- coord$ppp
    msk <- coord$masks[[mask]]
    filter <- sapply(1:length(dat$x),function(i,dat,msk)msk[dat$x[i],dat$y[i]]==1,dat,msk)
    coord$ppp <- coord$ppp[filter,]
    coord$raw$data <- coord$raw$data[filter,]
    return(coord)    
}

filter_sample <- function(samp,mask){
    samp <- lapply(samp,filter_coord,mask)
    return(samp)
}

extractMaskedData <- function(dataset,mask='filled_margin'){
    masked_set <- lapply(dataset,filter_sample,mask)
    return(masked_set)
}

################################################################################################
############# barplotter that plots a matrix with SE, each row a bar ###########################
################################################################################################


plot_bar_se <- function(mat,cols,main,ylim=NULL){
    means <- rowMeans(mat,na.rm = T)
    ses <- apply(mat,1,function(x) sd(x,na.rm = T)/sqrt(length(x)))

    if (is.null(main)){
        ylim=c(0,max(means)+max(ses))
    }
    bp <- barplot(means, 
                  main=main,
                  xlab='', 
                  ylab='Average normalized interactions', 
                  col=cols,
                  ylim=ylim,
                  las=2,
                  beside=TRUE)
    arrows(bp,
           means+ses,
           bp, 
           means, 
           angle=90, 
           code=3, 
           length=0.02,
           col='black')
}


