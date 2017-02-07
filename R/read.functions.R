###########################################
######### Thresholding functions


#This function assumes each sample has its own directory 
#where it has one or more images. The image names include the coordinate of the image
#in a format [xxxxx,yyyyy], which is the standard output of inForm.
#For each coordinate there are 3 files:
#   ***_cell_seg_data.txt
#   ***_cell_seg_data_summary.txt
#   ***_score_data.txt





#' Read inForm output and store it in an Iris object.
#' 
#' @param object Iris object.
#' @param raw_dir_name Directory that contains the raw files
#' @param label_fix List of length 2 character vector that is used to fix filenames.
#' @param format Output format: Currently only "Vectra" and "Mantra" are supported. 
#' @param dir_filter Filter to select only certain directory names.
#' @param read_nuc_seg_map Flag indicating whether the nuclear map should be read.
#' @param read_dapi_map Flag indicating whether the dapi layer should be read.
#' @param MicronsPerPixel Length of one pixel. Default: 0.496, corresponding to a 20x Mantra/Vectra images
#' 
#' @return Iris object.
#' @examples
#' raw_data <- Iris()
#' raw_data<- read.raw(raw_data,
#'                     raw_dir_name=system.file("extdata", package = "Iris"),
#'                     format='Mantra')
#' @docType methods
#' @export
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
                                MicronsPerPixel=0.496,
                                invasive_margin_in_px=80,
                                readMasks=T){
              object@microns_per_pixel=MicronsPerPixel
              raw_directories <- dir(raw_dir_name)
              object@samples <- lapply(raw_directories,function(x)Sample(sample_name=x))
              object@samples <- lapply(object@samples,read.raw.sample,raw_dir_name,label_fix,
                                       format,dir_filter,read_nuc_seg_map,read_dapi_map,invasive_margin_in_px,readMasks)
              names(object@samples) <- toupper(raw_directories)
              
              if(length(object@samples)==0){
                  stop('No images files found in :',raw_dir_name)
              }
              #automatically extract the counts
              object <- extract.counts(object)
              return(object)
})


#' Read inForm output and store it in an Iris object.
#' 
#' @param object Iris object.
#' @param raw_dir_name Directory that contains the raw files
#' @param label_fix List of length 2 character vector that is used to fix filenames.
#' @param format Output format: Currently only "Vectra" and "Mantra" are supported. 
#' @param dir_filter Filter to select only certain directory names.
#' @param read_nuc_seg_map Flag indicating whether the nuclear map should be read.
#' @param read_dapi_map Flag indicating whether the dapi layer should be read.
#' @param MicronsPerPixel Length of one pixel. Default: 0.496, corresponding to a 20x Mantra/Vectra images
#' 
#' @return Iris object.
#' @examples
#' raw_data <- Iris()
#' raw_data<- read.raw(raw_data,
#'                     raw_dir_name=system.file("extdata", package = "Iris"),
#'                     format='Mantra')
setGeneric("read.raw.sample", function(object, ...) standardGeneric("read.raw.sample"))
setMethod("read.raw.sample",
          signature = "Sample",
          definition = function(object,
                                raw_dir_name,
                                label_fix,
                                format,
                                dir_filter,
                                read_nuc_seg_map,
                                read_dapi_map,
                                invasive_margin_in_px,
                                readMasks){
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
              }else{
                  stop('Unknown image format')
              }
              object@coordinates <- lapply(coordinates,function(x)Coordinate(coordinate_name=x))
              object@coordinates <-lapply(object@coordinates, read.raw.coordinate, sample_dir, image_names,
                                          label_fix, format, read_nuc_seg_map, read_dapi_map,invasive_margin_in_px,readMasks)
                                  
              names(object@coordinates) <- coordinates
              return(object)
})

#' Read inForm output from a single coordinate
#' 
#' @importFrom tiff readTIFF
#' @importFrom spatstat owin
setGeneric("read.raw.coordinate", function(object,...) standardGeneric("read.raw.coordinate"))
setMethod("read.raw.coordinate",
          signature = "Coordinate",
          definition = function(object,
                                sample_dir,
                                image_names,
                                label_fix,
                                format,
                                read_nuc_seg_map,
                                read_dapi_map,
                                invasive_margin_in_px,
                                readMasks){
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
              object@raw@data <- object@raw@data[object@raw@data$Phenotype!='',]
                  
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
              
              if (readMasks){
                  #extract mask data
                  object <- extract_mask_data(object, img_names, sample_dir, object@coordinate_name, invasive_margin_in_px)
              }
              
              if (length(object@mask)>0){
                  di <- dim(object@mask[[1]])
                  object@size_in_px <- di[1]*di[2]
              }else{
                  object@size_in_px <- c(max(object@ppp$x),max(object@ppp$y))
              }
              return(object)
              
})


extract_mask <- function(filename){
    mask <- readTIFF(filename)
    mask <- as.matrix((mask[,,1]+mask[,,2]+mask[,,3])>0)
    mask <- t(mask)
    return(mask)    
}

#' Read and process all masks
#' 
#' @useDynLib Iris
#' @importFrom Rcpp sourceCpp
setGeneric("extract_mask_data", function(object,...) standardGeneric("extract_mask_data"))
setMethod("extract_mask_data",
          signature = "Coordinate",
          definition = function(object, img_names, sample_dir, coordinate_name, invasive_margin_in_px){
              
             #invasive margin
             inv_mar <- img_names[grep('_Invasive_Margin.tif',img_names)]
             if (length(inv_mar)==0){
                 stop(paste('_Invasive_Margin.tif for',coordinate_name, 'in',
                            sample_dir, 'does not exists!'))
             }
             inv_mar <- file.path(sample_dir, inv_mar)
             object@mask$invasive_margin <- extract_mask(inv_mar)
             object@mask$invasive_margin <- growMarginC(object@mask$invasive_margin, invasive_margin_in_px)
             
             #tumor mask
             tumor_tif <- img_names[grep('_Tumor.tif', img_names)]
             if (length(tumor_tif)==0){
                 stop(paste('_Tumor.tif for',coordinate_name, 'in',
                            sample_dir, 'does not exists!'))
             }
             tumor_tif <- file.path(sample_dir,tumor_tif)
             object@mask$tumor <- extract_mask(tumor_tif)
             object@mask$tumor[object@mask$invasive_margin > 0] <- 0
                 
             #also add the stroma
             object@mask$stroma <- object@mask$tumor
             object@mask$stroma[object@mask$stroma==0] <- 1
             object@mask$stroma[object@mask$tumor>0] <- 0
             object@mask$stroma[object@mask$invasive_margin>0] <- 0

             return(object)
})

###########################################
######### Thresholding functions

#' This function reads the manually determined thresholds of certain markers (e.g. PD1, PD-L1) and splits selected celltypes into marker+ and marker- celltypes.
#' 
#' @param object Iris object.
#' @param marker Name of the marker used in the score file.
#' @param marker_name corresponding name, which should be appended at the selected cell types.
#' @param base Vector of cell types for which the marker should be used.
#' @param pheno_name Name of the phenotype column to be used. Default from inForm is "Phenotype".
#' @param remove_blanks Flag that indicates whether or not not called cells are to be removed.
#' 
#' @return Iris object.
#' @examples
#' dataset <- threshold.dataset(raw_data,
#'                              marker='PD-Ligand-1 (Opal 690)',
#'                              marker_name='PDL1',
#'                              base=c('SOX10+'))
#' dataset <- threshold.dataset(dataset,
#'                              marker='PD-1 (Opal 540)',
#'                              marker_name='PD1',
#'                              base=c('CD8+','OTHER'))
#' @export
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
              object@coordinates <- lapply(object@coordinates,threshold_coords,marker,marker_name,
                                           base,pheno_name,object@sample_name,remove_blanks)
              names(object@coordinates) <- sapply(object@coordinates,function(x)x@coordinate_name)
              return(object)
})

setGeneric("threshold_coords", function(object,...) standardGeneric("threshold_coords"))
setMethod("threshold_coords",
          signature = "Coordinate",
          definition = function (object, marker, marker_name, base, pheno_name, sample_name, remove_blanks){
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
              if (!paste0(pheno_name,'.combined')%in%colnames(object@raw@data)){
                  object@raw@data[[paste0(pheno_name,'.combined')]] <- object@raw@data[[pheno_name]]
              }
              
              #get the thresholds for the marker we want to score
              scoring <- getScoring(object)
              #r automatically replaces special characters with '.' when they are used for names so I'm fixing this
              mark <- gsub('[ \\(\\)]','.',marker)
              mark <- gsub('-','.',mark)
              
              if(sum(scoring$Component==mark) == 0){
                 stop('Could not find Score for: ',marker,' for Sample: ',sample_name,' Coordinate: ',object@coordinate_name)
              }
              scoring <- scoring[scoring$Component==mark,]
            
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

########################################################################################################################
###ROI functions


#' Extract ROI
#' @export
#' 
#' 
setGeneric("extract.ROI", function(object, ...) standardGeneric("extract.ROI"))
setMethod("extract.ROI",
          signature = "Iris",
          definition = function(object, ROI='invasive_margin'){
              all_levels <- sort(unique(unlist(lapply(object@counts,colnames))))
              object@samples <- lapply(object@samples, extract.ROI.sample, ROI, all_levels)
              
              #update the counts
              object <- extract.counts(object)
              
              #reset all spatial stats
              object@nearest_neighbors <- list()
              object@interactions <- list()
              object@proximity <- list()
              
              return(object)
          })


setGeneric("extract.ROI.sample", function(object, ...) standardGeneric("extract.ROI.sample"))
setMethod("extract.ROI.sample",
          signature = "Sample",
          definition = function(object, ROI, all_levels){
              object@coordinates <- lapply(object@coordinates, extract.ROI.Coordinate, ROI, all_levels)
              return(object)
          })


setGeneric("extract.ROI.Coordinate", function(object,...) standardGeneric("extract.ROI.Coordinate"))
setMethod("extract.ROI.Coordinate",
          signature = "Coordinate",
          definition = function(object, ROI, all_levels){
              #reduce to the filter
              mask <- object@mask[[ROI]]
              filter <- sapply(1:length(object@ppp$x),function(i,dat,mask)mask[dat$x[i],dat$y[i]]==1,object@ppp,mask)
              object@ppp <- object@ppp[filter,]
              object@raw@data <- object@raw@data[filter,]
              object@size_in_px <- sum(mask>0) 
              
              return(object)
          })
