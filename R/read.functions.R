
#This function assumes each sample has its own directory 
#where it has one or more images. The image names include the coordinate of the image
#in a format [xxxxx,yyyyy], which is the standard output of inForm.
#For each coordinate there are 3 files:
#   ***_cell_seg_data.txt
#   ***_score_data.txt
#   ***_cell_seg_data_summary.txt (which is optional)


#' Read inForm output and store it in an Iris object.
#' 
#' @param x Iris object.
#' @param raw_dir_name Directory that contains the raw files
#' @param label_fix List of length 2 character vector that is used to fix filenames.
#' @param format Output format: Currently only "Vectra" and "Mantra" are supported. 
#' @param dir_filter Filter to select only certain directory names.
#' @param read_nuc_seg_map Flag indicating whether the nuclear map should be read.
#' @param read_dapi_map Flag indicating whether the dapi layer should be read.
#' @param MicronsPerPixel Length of one pixel. Default: 0.496, corresponding to a 20x Mantra/Vectra images
#' @param invasive_margin_in_px The width of the invasive margin in pixels
#' @param readMasks Flag indicating whether the "_Tumor.tif" and "_Invasive_Margin.tif" should be read (default: True)
#' @param ignore_scoring Flag indicating whether the scoring file should be ignored (default: False)
#' @param ... Additional arguments  
#' 
#' @return Iris object.
#' @examples
#' raw_data <- Iris()
#' raw_data <- read_raw(raw_data,
#'                      raw_dir_name=system.file("extdata", package = "Iris"),
#'                      format='Mantra')
#' @docType methods
#' @export
#' @importFrom methods new 
#' @rdname read_raw
setGeneric("read_raw", function(x, ...) standardGeneric("read_raw"))

#' @rdname read_raw
#' @aliases read_raw,ANY,ANY-method
setMethod("read_raw",
          signature = "Iris",
          definition = function(x,
                                raw_dir_name,
                                label_fix=list(),
                                format='Vectra',
                                dir_filter='',
                                read_nuc_seg_map=F,
                                read_dapi_map=F,
                                MicronsPerPixel=0.496,
                                invasive_margin_in_px=80,
                                readMasks=T,
                                ignore_scoring=F){
              x@microns_per_pixel=MicronsPerPixel
              raw_directories <- dir(raw_dir_name)
              x@samples <- lapply(raw_directories,function(x)Sample(sample_name=x))
              x@samples <- lapply(x@samples,read_raw_sample,raw_dir_name,label_fix,
                                       format,dir_filter,read_nuc_seg_map,read_dapi_map,
                                       invasive_margin_in_px, readMasks, ignore_scoring)
              names(x@samples) <- toupper(raw_directories)
              
              if(length(x@samples)==0){
                  stop('No images files found in :',raw_dir_name)
              }
              #automatically extract the counts
              x <- extract_counts(x)
              return(x)
})

setGeneric("read_raw_sample", function(x, ...) standardGeneric("read_raw_sample"))
setMethod("read_raw_sample",
          signature = "Sample",
          definition = function(x,
                                raw_dir_name,
                                label_fix,
                                format,
                                dir_filter,
                                read_nuc_seg_map,
                                read_dapi_map,
                                invasive_margin_in_px,
                                readMasks,
                                ignore_scoring){
              print(paste('Sample:', x@sample_name))
        
              #get sample directory
              sample_dir <- file.path(raw_dir_name,x@sample_name)
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
              x@coordinates <- lapply(coordinates,function(x)Coordinate(coordinate_name=x))
              x@coordinates <-lapply(x@coordinates, read_raw_coordinate, sample_dir, image_names,
                                          label_fix, format, read_nuc_seg_map, read_dapi_map,
                                          invasive_margin_in_px, readMasks, ignore_scoring)
                                  
              names(x@coordinates) <- coordinates
              return(x)
})

#' @importFrom tiff readTIFF
#' @importFrom spatstat owin
#' @importFrom utils read.csv
setGeneric("read_raw_coordinate", function(x, ...) standardGeneric("read_raw_coordinate"))
setMethod("read_raw_coordinate",
          signature = "Coordinate",
          definition = function(x,
                                sample_dir,
                                image_names,
                                label_fix,
                                format,
                                read_nuc_seg_map,
                                read_dapi_map,
                                invasive_margin_in_px,
                                readMasks,
                                ignore_scoring){
              if (format == 'Vectra'){
                  img_names <- image_names[grep(x@coordinate_name,image_names)]
              }else if (format == 'Mantra'){
                  file_parts <- image_names
                  if (length(grep('MULTI',x@coordinate_name))>0){
                      file_parts <- sub('^.+MULTI_','',file_parts)
                  }else{
                      file_parts <- sub('^[^_]+_','',file_parts)
                  }
                  img_names <- image_names[grep(paste0(x@coordinate_name,'_'),file_parts)]
              }   
                  
              seg_data <- img_names[grep('cell_seg_data.txt$',img_names)]
              if (length(seg_data)!=1){
                  stop('Could not find a single *_cell_seg_data.txt for ',x@coordinate_name,' in ',sample_dir)
              }
              #grab all of the data files and put them into a list
              x@raw@data <- read.csv(file.path(sample_dir,seg_data),
                                          sep='\t',
                                          as.is=T)

              x@raw@data <- x@raw@data[x@raw@data$Phenotype!='',]
              
              if (length(grep('_cell_seg_data_summary.txt$',img_names))>0){
                  x@raw@summary <- t(read.csv(file.path(sample_dir,
                                                             img_names[grep('_cell_seg_data_summary.txt$',img_names)]),
                                                   sep='\t',
                                                   as.is=T))
              }
                  
              if (length(grep('_memb_seg_map.tif',img_names))>0){
                  x@raw@mem_seg_map <- readTIFF(file.path(sample_dir,
                                                                img_names[grep('_memb_seg_map.tif',img_names)]))
              }
                  
              if (read_nuc_seg_map && length(grep('_nuc_seg_map.tif',img_names))>0){
                  x@raw@nuc_seg_map <- readTIFF(file.path(sample_dir,
                                                               img_names[grep('_nuc_seg_map.tif',img_names)]))
              }
                  
              if (read_dapi_map && length(grep('_component_data.tif',img_names))>0){
                  dapi <- readTIFF(file.path(sample_dir,
                                             img_names[grep('_component_data.tif',img_names)]),all = T)[[1]]
                  if (class(dapi) == 'array'){
                      dapi <- dapi[,,4]
                  }
                  x@raw@dapi_map <- dapi
              }
                  
              if (!ignore_scoring){
                  score_data <- img_names[grep('_score_data.txt$',img_names)]
                  if (length(score_data) != 1){
                      stop('Could not find a single *_score_data.txt for ',x@coordinate_name,' in ',sample_dir, 
                           'if no markers were scored please set ignore_scoring flag to TRUE.')
                  }    
                  x@raw@score <- t(read.csv(file.path(sample_dir,score_data),
                                                 sep='\t',
                                                 as.is=T))
              }
                  
              #fix the labels if necessary
              if(length(label_fix)>0){
                  #for each label fix
                  for (fix in label_fix){
                      x@raw@data$Phenotype[grep(fix[1],x@raw@data$Phenotype,fixed = T)] <- fix[2]
                  }
              }
              
              if (length(x@raw@data$memb_seg_map)>0){
                  x_max <- ncol(x@raw@memb_seg_map)
                  y_max <- nrow(x@raw@memb_seg_map)
              }else{
                  x_max <- max(x@raw@data$Cell.X.Position)
                  y_max <- max(x@raw@data$Cell.Y.Position)
              }
                  
              #estimate the window size
              window = owin(xrange = c(0, x_max), 
                            yrange = c(0, y_max))
              
              #sqeeze the data into the spatstats package format
              x@ppp  <- with(x@raw@data, ppp(`Cell.X.Position`, 
                                             `Cell.Y.Position`, 
                                             window = window, 
                                             marks=factor(x@raw@data$Phenotype)))
              
              if (readMasks){
                  #extract mask data
                  x <- extract_mask_data(x, img_names, sample_dir, x@coordinate_name, invasive_margin_in_px)
              }
              
              if (length(x@mask)>0){
                  di <- dim(x@mask[[1]])
                  x@size_in_px <- di[1]*di[2]
              }else{
                  x@size_in_px <- x_max * y_max
              }
              return(x)
              
})


#' Read inForm output from a single coordinate
#' 
#' @param filename Name of the .tif file that contains the mask.
#' 
#' @return Mask matrix
#' @export
extract_mask <- function(filename){
    mask <- readTIFF(filename)
    mask <- as.matrix((mask[,,1]+mask[,,2]+mask[,,3])>0)
    mask <- t(mask)
    return(mask)    
}


#' @useDynLib Iris
#' @importFrom Rcpp sourceCpp
setGeneric("extract_mask_data", function(x, ...) standardGeneric("extract_mask_data"))
setMethod("extract_mask_data",
          signature = "Coordinate",
          definition = function(x, img_names, sample_dir, coordinate_name, invasive_margin_in_px){
              
             #invasive margin
             inv_mar <- img_names[grep('_Invasive.Margin.tif',img_names)]
             if (length(inv_mar)==0){
                 stop(paste('_Invasive_Margin.tif for',coordinate_name, 'in',
                            sample_dir, 'does not exists!'))
             }
             inv_mar <- file.path(sample_dir, inv_mar)
             x@mask$margin_line <- extract_mask(inv_mar)
             if (sum(x@mask$margin_line)==0){
                 stop(paste('_Invasive.Margin.tif for',coordinate_name, 'in',
                            sample_dir, 'has no pixels, please check the mask file!'))
             }
             x@mask$invasive_margin <- growMarginC(x@mask$margin_line, invasive_margin_in_px)
             
             #tumor mask
             tumor_tif <- img_names[grep('_Tumor.tif', img_names)]
             if (length(tumor_tif)==0){
                 stop(paste('_Tumor.tif for',coordinate_name, 'in',
                            sample_dir, 'does not exists!'))
             }
             tumor_tif <- file.path(sample_dir,tumor_tif)
             x@mask$tumor <- extract_mask(tumor_tif)
             if (sum(x@mask$margin_line)==0){
                 stop(paste('_Tumor.tif for',coordinate_name, 'in',
                            sample_dir, 'has no pixels, please check the mask file!'))
             }
             if (!all(dim(x@mask$invasive_margin)==dim(x@mask$tumor))){
                 stop(paste('_Tumor.tif for',coordinate_name, 'in',
                            sample_dir, 'does not have the same dimension as the invasive margin mask!'))
             }
             x@mask$tumor[x@mask$invasive_margin > 0] <- 0
                 
             #also add the stroma
             x@mask$stroma <- x@mask$tumor
             x@mask$stroma[x@mask$stroma==0] <- 1
             x@mask$stroma[x@mask$tumor>0] <- 0
             x@mask$stroma[x@mask$invasive_margin>0] <- 0

             return(x)
})

###########################################
######### Thresholding functions

#' This function reads the manually determined thresholds of certain markers (e.g. PD1, PD-L1) and splits selected celltypes into marker+ and marker- celltypes.
#' 
#' @param x Iris object.
#' @param marker Name of the marker used in the score file.
#' @param marker_name corresponding name, which should be appended at the selected cell types.
#' @param base Vector of cell types for which the marker should be used.
#' @param pheno_name Name of the phenotype column to be used. (Default from inForm is "Phenotype")
#' @param remove_blanks Flag that indicates whether or not not called cells are to be removed. (Default: TRUE)
#' @param ... Additional arguments  
#' 
#' @docType methods
#' @return Iris object.
#' @examples
#' raw_data <- Iris()
#' raw_data <- read_raw(raw_data,
#'                      raw_dir_name=system.file("extdata", package = "Iris"),
#'                      format='Mantra')
#' dataset <- threshold_dataset(raw_data,
#'                              marker='PD-Ligand-1 (Opal 690)',
#'                              marker_name='PDL1',
#'                              base=c('SOX10+'))
#' dataset <- threshold_dataset(dataset,
#'                              marker='PD-1 (Opal 540)',
#'                              marker_name='PD1',
#'                              base=c('CD8+','OTHER'))
#' @export
#' @rdname threshold_dataset
setGeneric("threshold_dataset", function(x, ...) standardGeneric("threshold_dataset"))

#' @rdname threshold_dataset
#' @aliases threshold_dataset,ANY,ANY-method
setMethod("threshold_dataset",
          signature = "Iris",
          definition = function (x, 
                                 marker, 
                                 marker_name, 
                                 base=NULL, 
                                 pheno_name='Phenotype', 
                                 remove_blanks=T){
              #for each sample
              x@samples <- lapply(x@samples, threshold_samples, marker, marker_name, base, pheno_name, remove_blanks)
              names(x@samples) <- sapply(x@samples,function(x)x@sample_name)
              x <- extract_counts(x)
              return(x)
})

setGeneric("threshold_samples", function(x, ...) standardGeneric("threshold_samples"))
setMethod("threshold_samples",
          signature = "Sample",
          definition = function (x, marker, marker_name, base, pheno_name, remove_blanks){
              #for each coordinate
              x@coordinates <- lapply(x@coordinates,threshold_coords,marker,marker_name,
                                           base,pheno_name,x@sample_name,remove_blanks)
              names(x@coordinates) <- sapply(x@coordinates,function(x)x@coordinate_name)
              return(x)
})


setGeneric("threshold_coords", function(x, ...) standardGeneric("threshold_coords"))
setMethod("threshold_coords",
          signature = "Coordinate",
          definition = function (x, marker, marker_name, base, pheno_name, sample_name, remove_blanks){
              #remove the cells that are not called by inForm (usually not very many!)
              if (remove_blanks){
                  x@raw@data <- x@raw@data[x@raw@data[[pheno_name]]!='',] 
                  x@ppp <- x@ppp[x@ppp$marks!='',]
                  x@ppp$marks <- droplevels(x@ppp$marks)
              }
              
              #if no cell types were specified 
              if (is.null(base)){
                  base <- levels(x@ppp$marks)
              }
            
              #make a combined phenotype column in the rawdata
              if (!paste0(pheno_name,'.combined')%in%colnames(x@raw@data)){
                  x@raw@data[[paste0(pheno_name,'.combined')]] <- x@raw@data[[pheno_name]]
              }
              
              #get the thresholds for the marker we want to score
              scoring <- getScoring(x)
              #r automatically replaces special characters with '.' when they are used for names so I'm fixing this
              mark <- gsub('[ \\(\\)]','.',marker)
              mark <- gsub('-','.',mark)
              
              if(sum(scoring$Component==mark) == 0){
                 stop('Could not find Score for: ',marker,' for Sample: ',sample_name,' Coordinate: ',x@coordinate_name)
              }
              scoring <- scoring[scoring$Component==mark,]
            
              #extract the current marker expression
              expression <- x@raw@data[[paste(scoring$Compartment,
                                                   scoring$Component,
                                                   'Mean..Normalized.Counts..Total.Weighting.',
                                                   sep='.')]]
              #deterine the positive cells
              positive_cells <- expression > scoring$Threshold
            
              #use only the cells that are within base
              current <- x@raw@data[[pheno_name]] %in% base
              current_base <- x@raw@data$Phenotype.combined[current]
             
              #fetch the cells we are currently working on
              x@raw@data$Phenotype.combined[current & !positive_cells] <- paste0(current_base[!positive_cells[current]],' ',marker_name,'-')
              x@raw@data$Phenotype.combined[current & positive_cells] <- paste0(current_base[positive_cells[current]],' ',marker_name,'+')
              x@ppp$marks <- as.factor(x@raw@data$Phenotype.combined)
            
              return(x)
})

setGeneric("getScoring", function(x, ...) standardGeneric("getScoring"))
setMethod("getScoring",
          signature = "Coordinate",
          definition = function(x){
              scoring <- x@raw@score
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

#' Method that reduces the current dataset to a specific region of interest, discarding all cell coordinates outside of that region
#'  
#' @param x Iris object 
#' @param ROI Region of interest (default: 'invasive_margin')
#' @param ... Additional arguments  
#' 
#' @return Iris object
#' @examples
#' raw_data <- Iris()
#' raw_data <- read_raw(raw_data,
#'                      raw_dir_name=system.file("extdata", package = "Iris"),
#'                      format='Mantra')
#' dataset <- threshold_dataset(raw_data,
#'                              marker='PD-Ligand-1 (Opal 690)',
#'                              marker_name='PDL1',
#'                              base=c('SOX10+'))
#' dataset <- threshold_dataset(dataset,
#'                              marker='PD-1 (Opal 540)',
#'                              marker_name='PD1',
#'                              base=c('CD8+','OTHER'))                     
#' @docType methods
#' @export
#' @rdname extract_ROI
setGeneric("extract_ROI", function(x, ...) standardGeneric("extract_ROI"))

#' @rdname extract_ROI
#' @aliases extract_ROI,ANY,ANY-method
setMethod("extract_ROI",
          signature = "Iris",
          definition = function(x, ROI='invasive_margin'){
              if (length(x@samples[[1]]@coordinates[[1]]@mask[[ROI]])==0){
                  stop('There is no mask for "',ROI,'"')
              }
             
              x@samples <- lapply(x@samples, extract_ROI_sample, ROI)
              
              #update the counts
              x <- extract_counts(x)
              
              #reset all spatial stats
              x@nearest_neighbors <- list()
              x@interactions <- list()
              x@proximity <- list()
              
              return(x)
          })

setGeneric("extract_ROI_sample", function(x, ...) standardGeneric("extract_ROI_sample"))
setMethod("extract_ROI_sample",
          signature = "Sample",
          definition = function(x, ROI){
              x@coordinates <- lapply(x@coordinates, extract_ROI_Coordinate, ROI)
              return(x)
          })


setGeneric("extract_ROI_Coordinate", function(x, ...) standardGeneric("extract_ROI_Coordinate"))
setMethod("extract_ROI_Coordinate",
          signature = "Coordinate",
          definition = function(x, ROI){
              #reduce to the filter
              mask <- x@mask[[ROI]]
              filter <- sapply(1:length(x@ppp$x),function(i,dat,mask)mask[dat$x[i],dat$y[i]]==1,x@ppp,mask)
              x@ppp <- x@ppp[filter,]
              x@raw@data <- x@raw@data[filter,]
              x@size_in_px <- sum(mask>0) 
              
              return(x)
          })
