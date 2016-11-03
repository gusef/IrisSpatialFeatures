# Install and load required packages
required = c('ggplot2',
             'grid',
             'raster',
             'RColorBrewer',
             'spatstat',
             'tiff')
for (lib in required)
{
    if (!require(lib, character.only=TRUE))
    {
        install.packages(lib, repos="http://cran.rstudio.com")
        suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
    }
}

#########################################################################################
######################## File reading function ##########################################
#########################################################################################

#This function assumes each sample has its own directory 
#where it has one or more images. The image names include the coordinate of the image
#in a format [xxxxx,yyyyy], which is the standard output of inForm.
#For each coordinate there are 3 files:
#   ***_cell_seg_data.txt
#   ***_cell_seg_data_summary.txt
#   ***_score_data.txt
#Input: raw directory
#       list of label fixes, where each list element has 2 strings to change [from,to]
readRawData <- function(raw_dir_name='raw_data',label_fix,format='Vectra',dir_filter=''){
    #parse the sample names
    raw_directories <- dir(raw_dir_name)
    raw_data <- lapply(raw_directories,extract_sample,raw_dir_name,label_fix,format,dir_filter)
    names(raw_data) <- toupper(raw_directories)
    return(raw_data)
}



#function to extract a single sample, figures out the different coordinate based on []
extract_sample <- function(sample,raw_dir_name,label_fix,format,dir_filter){
    cat('Sample:', sample,'\n')
    #get sample directory
    
    sample_dir <- file.path(raw_dir_name,sample)
    image_names <- dir(sample_dir,recursive = T)
    
    #directory filter in case there are different projects for each sample
    if (dir_filter!=''){
        image_names <- image_names[grep(dir_filter,image_names)]
    }

    #figure out the different coordinates for each sample    
    if (format == 'Vectra'){
        image_names <- image_names[grep('\\[.*\\]',image_names)]
        coordinates <- unique(sub('\\].+$','',sub('^.+\\[','',image_names)))
    }else if (format == 'Mantra'){
        coords <- image_names[grep("_cell_seg_data.txt",image_names)]
        coords <- sub("_cell_seg_data.txt",'',coords)
        if (length(grep('MULTI',coords))>0){
            coordinates <- sub('^.+MULTI_','',coords)
        }else{
            coordinates <- sub('^[^_]+_','',coords)
        }
    }

    #extract the data for each coordinate
    coord_data <- lapply(coordinates,extract_coordinate,sample_dir,image_names,label_fix,format)
    #label each coordinate appropriately and return
    names(coord_data) <- coordinates
    return(coord_data)
}

extract_mask <- function(filename){
    mask <- readTIFF(filename)
    mask <- as.matrix((mask[,,1]+mask[,,2]+mask[,,3])>0)
    mask <- t(mask)
    return(mask)    
}

#this is just a function stub, which can be overwritten to read masks for specific use-cases.
extract_mask_data <- function(img_names, sample_dir){
    return(NULL)
}


#extracts all the raw data (output from inform) and generates a spatstat ppp class object
extract_coordinate <- function(coord,sample_dir,image_names,label_fix,format){
    #focus only on the files for the current coordinate
    if (format == 'Vectra'){
        img_names <- image_names[grep(coord,image_names)]
    }else if (format == 'Mantra'){
        #file_parts <- sapply(strsplit(sapply(strsplit(image_names,'_cell'),function(x)x[1]),'_score'),function(y)y[1])
        file_parts <- image_names
        if (length(grep('MULTI',coord))>0){
            file_parts <- sub('^.+MULTI_','',file_parts)
        }else{
            file_parts <- sub('^[^_]+_','',file_parts)
        }
        img_names <- image_names[grep(paste0(coord,'_'),file_parts)]
    }   
    
    #grab all of the data files and put them into a list
    image_data <- list()
    image_data$data <- read.csv(file.path(sample_dir,
                                          img_names[grep('cell_seg_data.txt$',img_names)]),
                                sep='\t',
                                as.is=T)
    
    if (length(grep('_cell_seg_data_summary.txt$',img_names))>0){
        image_data$summary <- t(read.csv(file.path(sample_dir,
                                                   img_names[grep('_cell_seg_data_summary.txt$',img_names)]),
                                         sep='\t',
                                         as.is=T))
    }

    if (length(grep('_memb_seg_map.tif',img_names))>0){
        image_data$memb_seg_map <- readTIFF(file.path(sample_dir,
                                                      img_names[grep('_memb_seg_map.tif',img_names)]))
    }
    
    if (length(grep('_nuc_seg_map.tif',img_names))>0){
        image_data$nuc_seg_map <- readTIFF(file.path(sample_dir,
                                                     img_names[grep('_nuc_seg_map.tif',img_names)]))
    }
    
    if (length(grep('_component_data.tif',img_names))>0){
        dapi <- readTIFF(file.path(sample_dir,
                           img_names[grep('_component_data.tif',img_names)]),all = T)[[1]]
        if (class(dapi) == 'array'){
            dapi <- dapi[,,4]
        }
        image_data$dapi_map <- dapi
    }
    
    image_data$score <- t(read.csv(file.path(sample_dir,
                                             img_names[grep('_score_data.txt$',img_names)]),
                                   sep='\t',
                                   as.is=T))
    
    #fix the labels if necessary
    if(length(label_fix)>0){
        #for each label fix
        for (fix in label_fix){
            image_data$data$Phenotype[grep(fix[1],image_data$data$Phenotype,fixed = T)] <- fix[2]
        }
    }
    
    if (length(image_data$memb_seg_map)>0){
        x_max <- ncol(image_data$memb_seg_map)
        y_max <- nrow(image_data$memb_seg_map)
    }else{
        x_max <- max(image_data$data$Cell.X.Position)
        y_max <- max(image_data$data$Cell.Y.Position)
    }
    
    #estimate the window size
    window = owin(xrange = c(0, x_max), 
                  yrange = c(0, y_max))
    
    #sqeeze the data into the spatstats package format
    ppp  <- with(image_data$data, ppp(`Cell.X.Position`, 
                                      `Cell.Y.Position`, 
                                      window = window, 
                                      marks=factor(image_data$data$Phenotype)))
    
    #extract mask data
    masks <- extract_mask_data(img_names, sample_dir)
    
    
    return(list(ppp=ppp,raw=image_data,masks=masks))
}

########################################################################################
############### Thresholding functions #################################################
########################################################################################

#Apply the thresholds on the entire dataset
#raw data: is the object generated by readRawData
#marker: is the name of the marker used in the score file
#marker_name: is the corresponding name for the marker that should be used in the Phenotype column
#base: is the cells for which a certain marker should be called
#pheno_name: is the name of the phenotype column to be used
#remove_blanks: indicates whether or not not called cells are to be removed
threshold_dataset <- function (raw_data,marker,marker_name,base=NULL,pheno_name='Phenotype',remove_blanks=T){
    #for each sample
    thresholded <- lapply(raw_data, threshold_samples,marker,marker_name,base,pheno_name,remove_blanks)
    names(thresholded) <- names(raw_data)
    return(thresholded)
}

#apply threshold to all coordinates of a single sample
threshold_samples <- function (sample,marker,marker_name,base,pheno_name,remove_blanks){
    #for each coordinate
    thresholded <- lapply(sample,threshold_coords,marker,marker_name,base,pheno_name,remove_blanks)
    names(thresholded) <- names(sample)
    return(thresholded)
}

#apply threshold to a single coordinate
threshold_coords <- function (coord,marker,marker_name,base,pheno_name,remove_blanks){
    #remove the cells that are not called by inForm (usually not very many!)
    if (remove_blanks){
        coord$raw$data <- coord$raw$data[coord$raw$data[[pheno_name]]!='',] 
        coord$ppp <- coord$ppp[coord$ppp$marks!='',]
        coord$ppp$marks <- droplevels(coord$ppp$marks)
    }
    
    #if no cell types were specified 
    if (is.null(base)){
        base <- levels(coord$ppp$marks)
    }
    
    #make a combined phenotype column in the rawdata
    if (!'Phenotype.combined'%in%colnames(coord$raw$data)){
        coord$raw$data$Phenotype.combined <- coord$raw$data[[pheno_name]]
    }
    
    #get the thresholds for the marker we want to score
    scoring <- getScoring(coord$raw$score)
    #r automatically replaces special characters with '.' when they are used for names so I'm fixing this
    marker <- gsub('[ \\(\\)]','.',marker)
    marker <- gsub('-','.',marker)
    
    scoring <- scoring[scoring$Component==marker,]
    
    #extract the current marker expression
    expression <- coord$raw$data[[paste(scoring$Compartment,
                                        scoring$Component,
                                        'Mean..Normalized.Counts..Total.Weighting.',
                                        sep='.')]]
    #deterine the positive cells
    positive_cells <- expression > scoring$Threshold
    
    #use only the cells that are within base
    current <- coord$raw$data[[pheno_name]] %in% base
    current_base <- coord$raw$data$Phenotype.combined[current]
    
    #fetch the cells we are currently working on
    coord$raw$data$Phenotype.combined[current & !positive_cells] <- paste0(current_base[!positive_cells[current]],' ',marker_name,'-')
    coord$raw$data$Phenotype.combined[current & positive_cells] <- paste0(current_base[positive_cells[current]],' ',marker_name,'+')
    coord$ppp$marks <- as.factor(coord$raw$data$Phenotype.combined)
    
    return(coord)
}

#extracts the scoring information automatically, for all scored markers
#works well with all of Chris' data, might cause some problems if the score file is not consistent
getScoring <- function(scoring){
    
    scores <- matrix(nrow=0,ncol=3)
    colnames(scores) <- c('Compartment','Component','Threshold')
    #if there were more than one additional markers scores:    
    if (length(grep('First',rownames(scoring))>0)){
        tab <- c('First','Second','Third')
        for (i in seq(length(tab))){
            #only if indicator actually exists (as of now I'm not sure if inForm allows for more than 2 markers)
            if (length(grep(tab[i],rownames(scoring))>0)){
                compartment <- scoring[paste0(tab[i],'.Cell.Compartment'),1]
                component <- scoring[paste0(tab[i],'.Stain.Component'),1]
                #r automatically replaces special characters with '.' when they are used for names so I'm fixing this
                component <- gsub('[ \\(\\)]','.',component)
                component <- gsub('-','.',component)
                threshold <- scoring[paste0(component,'.Threshold'),1]
                scores <- rbind(scores,c(compartment,component,threshold))
            }
        }
    #if there was only one marker scored
    }else{
        compartment <- scoring['Cell.Compartment',1]
        component <- scoring['Stain.Component',1]
        #r automatically replaces special characters with '.' when they are used for names so I'm fixing this
        component <- gsub('[ \\(\\)]','.',component)
        component <- gsub('-','.',component)
        threshold <- scoring['Positivity.Threshold',1]
        scores <- rbind(scores,c(compartment,component,threshold))
    }
    scores <- data.frame(scores,stringsAsFactors = F)
    scores$Threshold <- as.numeric(scores$Threshold)
    return(scores)
}

########################################################################################
##################### Post processing functions#########################################
########################################################################################


extract_raw <- function(coord){
    raw <- list(data= coord$raw$data,
                score= coord$raw$score)
    return(raw)
}

extract_ppp <- function(coord){
    ppp <- list(ppp=coord$ppp)
    return(ppp)
}


exclude_low_freq <- function(sample,no_cell_exclusion){
    for (idx in seq(length(sample))){
        tab <- table(sample[[idx]]$ppp$marks)
        exclude <- names(tab)[tab < no_cell_exclusion]
        if (length(exclude)>0){
            #modify ppp object
            exclude <- sample[[idx]]$ppp$marks %in% exclude
            sample[[idx]]$ppp <- sample[[idx]]$ppp[!exclude,]
            sample[[idx]]$ppp$marks <- droplevels(sample[[idx]]$ppp$marks)
            sample[[idx]]$raw$data <- sample[[idx]]$raw$data[!exclude,]
        }
    }
    return(sample)
}

splitDataset <- function(dataset,no_cell_exclusion=15,prefix=''){
    dir.create(file.path(getwd(),'datasets'), showWarnings = FALSE)
    nuc_map <- lapply(dataset,function(x)lapply(x,function(y)y$raw$nuc_seg_map))
    saveRDS(nuc_map,file=paste0('datasets/',prefix,'nuc_map.RDS'))
    mem_map <- lapply(dataset,function(x)lapply(x,function(y)y$raw$memb_seg_map))
    saveRDS(mem_map,file=paste0('datasets/',prefix,'mem_map.RDS'))
    if (length(dataset[[1]][[1]]$raw$dapi_map)>0){
        dapi_map <- lapply(dataset,function(x)lapply(x,function(y)y$raw$dapi_map))
        saveRDS(dapi_map,file=paste0('datasets/',prefix,'dapi_images.RDS'))
    }
    #exclude all the cell types that have less than 15 of a cell-count
    dataset <- lapply(dataset, exclude_low_freq, no_cell_exclusion)
    
    raw_data <- lapply(dataset,function(x)lapply(x,extract_raw))
    saveRDS(raw_data,file=paste0('datasets/',prefix,'raw_data.RDS'))
    ppp <- lapply(dataset,function(x)lapply(x,extract_ppp))
    saveRDS(ppp,file=paste0('datasets/',prefix,'ppp.RDS'))
    
    #extract the masks
    if (length(dataset[[1]][[1]]$masks)>0){
        masks <- lapply(dataset,function(x)lapply(x,function(y)y$masks))
        saveRDS(masks,file=paste0('datasets/',prefix,'masks.RDS'))
    }
}



########################################################################################
##################### Writing function #################################################
########################################################################################

#write all thresholded sets
output_thresholded <- function(dataset,outdir){
    if (!file.exists(outdir)){
        dir.create(outdir)
    }
    
    #for each sample
    for (samp in names(dataset)){
        #for each coordinate
        for (coord in names(dataset[[samp]])){
            data <- dataset[[samp]][[coord]]$raw$data
            #fix the column names
            colnames(data) <- sub('^Phenotype$','Phenotype.inForm',colnames(data))
            colnames(data) <- sub('^Phenotype.combined$','Phenotype',colnames(data))
            #output the file
            write.table(data,
                        file = file.path(outdir,paste0(samp,'_[',coord,']_cell_seq_data.txt')),
                        quote = F, 
                        row.names = F,
                        sep = '\t')
        }
    }
}


