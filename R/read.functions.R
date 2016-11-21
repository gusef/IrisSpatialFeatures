###########################################
######### Thresholding functions


#This function assumes each sample has its own directory 
#where it has one or more images. The image names include the coordinate of the image
#in a format [xxxxx,yyyyy], which is the standard output of inForm.
#For each coordinate there are 3 files:
#   ***_cell_seg_data.txt
#   ***_cell_seg_data_summary.txt
#   ***_score_data.txt

setGeneric("read.raw", function(object, ...) standardGeneric("read.raw"))
setMethod("read.raw",
          signature = "Iris",
          definition = function(object,
                                raw_dir_name,
                                label_fix=list(),
                                format='Vectra',
                                dir_filter='',
                                read_nuc_seg_map=F,
                                read_dapi_map=F,
                                MicronsPerPixel=0.496){
              object@microns_per_pixel=MicronsPerPixel
              raw_directories <- dir(raw_dir_name)
              object@samples <- lapply(raw_directories,function(x)Sample(sample_name=x))
              object@samples <- lapply(object@samples,read.raw.sample,raw_dir_name,label_fix,
                                       format,dir_filter,read_nuc_seg_map,read_dapi_map)
              names(object@samples) <- toupper(raw_directories)
              
              if(length(object@samples)==0){
                  stop('No images files found in :',raw_dir_name)
              }
              #automatically extract the counts
              object <- extract.counts(object)
              return(object)
})



setGeneric("read.raw.sample", function(object, ...) standardGeneric("read.raw.sample"))
setMethod("read.raw.sample",
          signature = "Sample",
          definition = function(object,
                                raw_dir_name,
                                label_fix,
                                format,
                                dir_filter,
                                read_nuc_seg_map,
                                read_dapi_map){
              print(paste('Sample:', object@sample_name))
        
              #get sample directory
              sample_dir <- file.path(raw_dir_name,object@sample_name)
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
              object@coordinates <- lapply(coordinates,function(x)Coordinate(coordinate_name=x))
              object@coordinates <-lapply(object@coordinates, read.raw.coordinate, sample_dir, image_names,
                                          label_fix, format, read_nuc_seg_map, read_dapi_map)
                                  
              names(object@coordinates) <- coordinates
              return(object)
})

setGeneric("read.raw.coordinate", function(object,...) standardGeneric("read.raw.coordinate"))
setMethod("read.raw.coordinate",
          signature = "Coordinate",
          definition = function(object,
                                sample_dir,
                                image_names,
                                label_fix,
                                format,
                                read_nuc_seg_map,
                                read_dapi_map){
              if (format == 'Vectra'){
                  img_names <- image_names[grep(object@coordinate_name,image_names)]
              }else if (format == 'Mantra'){
                  file_parts <- image_names
                  if (length(grep('MULTI',object@coordinate_name))>0){
                      file_parts <- sub('^.+MULTI_','',file_parts)
                  }else{
                      file_parts <- sub('^[^_]+_','',file_parts)
                  }
                  img_names <- image_names[grep(paste0(object@coordinate_name,'_'),file_parts)]
              }   
                  
              #grab all of the data files and put them into a list
              object@raw@data <- read.csv(file.path(sample_dir,
                                                    img_names[grep('cell_seg_data.txt$',img_names)]),
                                          sep='\t',
                                          as.is=T)
                  
              if (length(grep('_cell_seg_data_summary.txt$',img_names))>0){
                  object@raw@summary <- t(read.csv(file.path(sample_dir,
                                                             img_names[grep('_cell_seg_data_summary.txt$',img_names)]),
                                                   sep='\t',
                                                   as.is=T))
              }
                  
              if (length(grep('_memb_seg_map.tif',img_names))>0){
                  object@raw@mem_seg_map <- readTIFF(file.path(sample_dir,
                                                                img_names[grep('_memb_seg_map.tif',img_names)]))
              }
                  
              if (read_nuc_seg_map && length(grep('_nuc_seg_map.tif',img_names))>0){
                  object@raw@nuc_seg_map <- readTIFF(file.path(sample_dir,
                                                               img_names[grep('_nuc_seg_map.tif',img_names)]))
              }
                  
              if (read_dapi_map && length(grep('_component_data.tif',img_names))>0){
                  dapi <- readTIFF(file.path(sample_dir,
                                             img_names[grep('_component_data.tif',img_names)]),all = T)[[1]]
                  if (class(dapi) == 'array'){
                      dapi <- dapi[,,4]
                  }
                  object@raw@dapi_map <- dapi
              }
                  
              object@raw@score <- t(read.csv(file.path(sample_dir,
                                                       img_names[grep('_score_data.txt$',img_names)]),
                                             sep='\t',
                                             as.is=T))
                  
              #fix the labels if necessary
              if(length(label_fix)>0){
                  #for each label fix
                  for (fix in label_fix){
                      object@raw@data$Phenotype[grep(fix[1],object@raw@data$Phenotype,fixed = T)] <- fix[2]
                  }
              }
              
              if (length(object@raw@data$memb_seg_map)>0){
                  x_max <- ncol(object@raw@memb_seg_map)
                  y_max <- nrow(object@raw@memb_seg_map)
              }else{
                  x_max <- max(object@raw@data$Cell.X.Position)
                  y_max <- max(object@raw@data$Cell.Y.Position)
              }
                  
              #estimate the window size
              window = owin(xrange = c(0, x_max), 
                            yrange = c(0, y_max))
              
              #sqeeze the data into the spatstats package format
              object@ppp  <- with(object@raw@data, ppp(`Cell.X.Position`, 
                                                       `Cell.Y.Position`, 
                                                   window = window, 
                                                   marks=factor(object@raw@data$Phenotype)))
              
              #extract mask data
              object <- extract_mask_data(object, img_names, sample_dir)
              return(object)
              
})


extract_mask <- function(filename){
    mask <- readTIFF(filename)
    mask <- as.matrix((mask[,,1]+mask[,,2]+mask[,,3])>0)
    mask <- t(mask)
    return(mask)    
}

setGeneric("extract_mask_data", function(object,...) standardGeneric("extract_mask_data"))
setMethod("extract_mask_data",
          signature = "Coordinate",
          definition = function(object, img_names, sample_dir){
             
             inv_mar <- img_names[grep('_Invasive Margin.tif',img_names)]
             if (length(inv_mar)>0){
                 object@mask$invasive_margin <- extract_mask(file.path(sample_dir, inv_mar))
                 object@mask$filled_margin <- growMarginC(object@mask$invasive_margin, 80)
             }
             
             tumor_tif <- img_names[grep('_Tumor.tif', img_names)]
             if (length(tumor_tif) > 0){
                 object@mask$tumor <- extract_mask(file.path(sample_dir,tumor_tif))
             }
             
             if (length(tumor_tif) > 0 && length(inv_mar) > 0){
                 object@mask$tumor[object@mask$filled_margin > 0] <- 0
             }
             return(object)
    
})

###########################################
######### Thresholding functions

#Apply the thresholds on the entire dataset
#raw data: is the object generated by readRawData
#marker: is the name of the marker used in the score file
#marker_name: is the corresponding name for the marker that should be used in the Phenotype column
#base: is the cells for which a certain marker should be called
#pheno_name: is the name of the phenotype column to be used
#remove_blanks: indicates whether or not not called cells are to be removed

setGeneric("threshold.dataset", function(object,...) standardGeneric("threshold.dataset"))
setMethod("threshold.dataset",
          signature = "Iris",
          definition = function (object, 
                                 marker, 
                                 marker_name, 
                                 base=NULL, 
                                 pheno_name='Phenotype', 
                                 remove_blanks=T){
              #for each sample
              object@samples <- lapply(object@samples, threshold_samples, marker, marker_name, base, pheno_name, remove_blanks)
              names(object@samples) <- sapply(object@samples,function(x)x@sample_name)
              object <- extract.counts(object)
              return(object)
})


setGeneric("threshold_samples", function(object,...) standardGeneric("threshold_samples"))
setMethod("threshold_samples",
          signature = "Sample",
          definition = function (object, marker, marker_name, base, pheno_name, remove_blanks){
              #for each coordinate
              object@coordinates <- lapply(object@coordinates,threshold_coords,marker,marker_name,base,pheno_name,remove_blanks)
              names(object@coordinates) <- sapply(object@coordinates,function(x)x@coordinate_name)
              return(object)
})

setGeneric("threshold_coords", function(object,...) standardGeneric("threshold_coords"))
setMethod("threshold_coords",
          signature = "Coordinate",
          definition = function (object, marker, marker_name, base, pheno_name, remove_blanks){
              #remove the cells that are not called by inForm (usually not very many!)
              if (remove_blanks){
                  object@raw@data <- object@raw@data[object@raw@data[[pheno_name]]!='',] 
                  object@ppp <- object@ppp[object@ppp$marks!='',]
                  object@ppp$marks <- droplevels(object@ppp$marks)
              }
              
              #if no cell types were specified 
              if (is.null(base)){
                  base <- levels(object@ppp$marks)
              }
            
              #make a combined phenotype column in the rawdata
              if (!'Phenotype.combined'%in%colnames(object@raw@data)){
                  object@raw@data$Phenotype.combined <- object@raw@data[[pheno_name]]
              }
              
              #get the thresholds for the marker we want to score
              scoring <- getScoring(object)
              #r automatically replaces special characters with '.' when they are used for names so I'm fixing this
              marker <- gsub('[ \\(\\)]','.',marker)
              marker <- gsub('-','.',marker)
              
              scoring <- scoring[scoring$Component==marker,]
            
              #extract the current marker expression
              expression <- object@raw@data[[paste(scoring$Compartment,
                                                   scoring$Component,
                                                   'Mean..Normalized.Counts..Total.Weighting.',
                                                   sep='.')]]
              #deterine the positive cells
              positive_cells <- expression > scoring$Threshold
            
              #use only the cells that are within base
              current <- object@raw@data[[pheno_name]] %in% base
              current_base <- object@raw@data$Phenotype.combined[current]
             
              #fetch the cells we are currently working on
              object@raw@data$Phenotype.combined[current & !positive_cells] <- paste0(current_base[!positive_cells[current]],' ',marker_name,'-')
              object@raw@data$Phenotype.combined[current & positive_cells] <- paste0(current_base[positive_cells[current]],' ',marker_name,'+')
              object@ppp$marks <- as.factor(object@raw@data$Phenotype.combined)
            
              return(object)
})


setGeneric("getScoring", function(object,...) standardGeneric("getScoring"))
setMethod("getScoring",
          signature = "Coordinate",
          definition = function(object){
              scoring <- object@raw@score
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
})




