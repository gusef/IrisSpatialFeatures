#This function assumes each sample has its own directory
#where it has one or more images. The image names include the coordinate of the image
#in a format [xxxxx,yyyyy], which is the standard output of inForm.
#For each coordinate there are 3 files:
#   ***_cell_seg_data.txt
#   ***_score_data.txt
#   ***_cell_seg_data_summary.txt (which is optional)


#' Read inForm output and store it in an IrisSpatialFeatures ImageSet object.
#'
#' @param path Directory that contains the raw files
#' @param label_fix List of length 2 character vector that is used to fix filenames.
#' @param format Output format: Currently only "Vectra" and "Mantra" are supported.
#' @param dir_filter Filter to select only certain directory names.
#' @param read_nuc_seg_map Flag indicating whether the nuclear map should be read.
#' @param read_dapi_map Flag indicating whether the dapi layer should be read.
#' @param use_binary_seg_maps Flag indicating whether the dapi layer should be read.
#' @param MicronsPerPixel Length of one pixel. Default: 0.496, corresponding to a 20x Mantra/Vectra images
#' @param invasive_margin_in_px The width of the invasive margin in pixels
#' @param readMasks Flag indicating whether the "_Tumor.tif" and "_Invasive_Margin.tif" should be read (default: True)
#' @param ROI Flag indicating whether ROI was specified so the image is only analyzed within that region a '_ROI.tif' needs to be present. (default: False)
#' @param ignore_scoring Flag indicating whether the scoring file should be ignored (default: False)
#' @param read_only_relevant_markers Flag that indicates whether all read inform output should be kept or only the relevant markers
#'
#' @return IrisSpatialFeatures ImageSet object.
#' @examples
#'  raw_data <- new("ImageSet")
#'  raw_data <- read_raw(path=system.file("extdata", package = "IrisSpatialFeatures"),
#'                       format='Mantra')
#' @docType methods
#' @export
#' @importFrom methods new
#' @rdname read_raw
setGeneric("read_raw",
           function(path,
                    label_fix = list(),
                    format = 'Vectra',
                    dir_filter = '',
                    read_nuc_seg_map = FALSE,
                    read_dapi_map = FALSE,
                    use_binary_seg_maps = FALSE,
                    MicronsPerPixel = 0.496,
                    invasive_margin_in_px = 100,
                    readMasks = TRUE,
                    ROI = FALSE,
                    ignore_scoring = FALSE,
                    read_only_relevant_markers = TRUE,
                    ...) standardGeneric("read_raw"),
           valueClass = "ImageSet")

#' @rdname read_raw
#' @aliases read_raw,ANY,ANY-method
setMethod(
    "read_raw",
    signature = c(path="character"),
    definition = function(path,
                          label_fix,
                          format,
                          dir_filter,
                          read_nuc_seg_map,
                          read_dapi_map,
                          use_binary_seg_maps,
                          MicronsPerPixel,
                          invasive_margin_in_px,
                          readMasks,
                          ROI,
                          ignore_scoring,
                          read_only_relevant_markers) {
        x <- new("ImageSet")
        x@microns_per_pixel = MicronsPerPixel
        raw_directories <- dir(path)
        x@samples <-
            lapply(raw_directories, function(x)
                Sample(sample_name = x))
        x@samples <-
            lapply(
                x@samples,
                read_raw_sample,
                path,
                label_fix,
                format,
                dir_filter,
                read_nuc_seg_map,
                read_dapi_map,
                use_binary_seg_maps,
                invasive_margin_in_px,
                readMasks,
                ROI,
                ignore_scoring,
                read_only_relevant_markers
            )
        names(x@samples) <- toupper(raw_directories)

        if (length(x@samples) == 0) {
            stop('No images files found in :', path)
        }
        #automatically extract the counts
        x <- extract_counts(x)
        return(x)
    }
)

setGeneric("read_raw_sample", function(x, ...)
    standardGeneric("read_raw_sample"))
setMethod(
    "read_raw_sample",
    signature = "Sample",
    definition = function(x,
                          raw_dir_name,
                          label_fix,
                          format,
                          dir_filter,
                          read_nuc_seg_map,
                          read_dapi_map,
                          use_binary_seg_maps,
                          invasive_margin_in_px,
                          readMasks,
                          ROI,
                          ignore_scoring,
                          read_only_relevant_markers) {
        print(paste('Sample:', x@sample_name))

        #get sample directory
        sample_dir <- file.path(raw_dir_name, x@sample_name)
        image_names <- dir(sample_dir, recursive = TRUE)

        #directory filter in case there are different projects for each sample
        if (dir_filter != '') {
            image_names <- image_names[grep(dir_filter, image_names)]
        }

        #figure out the different coordinates for each sample
        if (format == 'Vectra') {
            image_names <- image_names[grep('\\[.*\\]', image_names)]
            coordinates <-
                unique(sub('\\].+$', '', sub('^.+\\[', '', image_names)))
        } else if (format == 'Mantra') {
            coords <- image_names[grep("_cell_seg_data.txt", image_names)]
            coords <- sub("_cell_seg_data.txt", '', coords)
            if (length(grep('MULTI', coords)) > 0) {
                coordinates <- sub('^.+MULTI_', '', coords)
            } else{
                coordinates <- sub('^[^_]+_', '', coords)
            }
        } else{
            stop('Unknown image format')
        }
        x@coordinates <-
            lapply(coordinates, function(x)
                Coordinate(coordinate_name = x))
        x@coordinates <-
            lapply(
                x@coordinates,
                read_raw_coordinate,
                sample_dir,
                image_names,
                label_fix,
                format,
                read_nuc_seg_map,
                read_dapi_map,
                use_binary_seg_maps,
                invasive_margin_in_px,
                readMasks,
                ROI,
                ignore_scoring,
                read_only_relevant_markers
            )

        names(x@coordinates) <- coordinates
        return(x)
    }
)

#' @importFrom tiff readTIFF
#' @importFrom spatstat owin
#' @importFrom utils read.csv
setGeneric("read_raw_coordinate", function(x, ...)
    standardGeneric("read_raw_coordinate"))
setMethod(
    "read_raw_coordinate",
    signature = "Coordinate",
    definition = function(x,
                          sample_dir,
                          image_names,
                          label_fix,
                          format,
                          read_nuc_seg_map,
                          read_dapi_map,
                          use_binary_seg_maps,
                          invasive_margin_in_px,
                          readMasks,
                          ROI,
                          ignore_scoring,
                          read_only_relevant_markers) {
        if (format == 'Vectra') {
            img_names <- image_names[grep(x@coordinate_name, image_names)]
        } else if (format == 'Mantra') {
            file_parts <- image_names
            if (length(grep('MULTI', x@coordinate_name)) > 0) {
                file_parts <- sub('^.+MULTI_', '', file_parts)
            } else{
                file_parts <- sub('^[^_]+_', '', file_parts)
            }
            img_names <-
                image_names[grep(paste0(x@coordinate_name, '_'), file_parts)]
        }

        seg_data <-
            img_names[grep('cell_seg_data.txt$', img_names)]
        if (length(seg_data) != 1) {
            stop(
                'Could not find a single *_cell_seg_data.txt for ',
                x@coordinate_name,
                ' in ',
                sample_dir
            )
        }

        #grab all of the data files and put them into a list
        dat <- read.csv(file.path(sample_dir, seg_data),
                               sep = '\t',
                               as.is = TRUE)

        #use only entries that have a phenotype assigned
        x@raw@data <- dat[dat$Phenotype != '', ]

        if (length(grep('_cell_seg_data_summary.txt$', img_names)) >
            0) {
            x@raw@summary <- t(read.csv(
                file.path(sample_dir,
                          img_names[grep('_cell_seg_data_summary.txt$', img_names)]),
                sep = '\t',
                as.is = TRUE
            ))
        }


        if (use_binary_seg_maps) {
            maps <- readTIFF(file.path(sample_dir,
                                       img_names[grep('_binary_seg_maps.tif', img_names)]), all =
                                 TRUE)
            x@raw@mem_seg_map <- maps[[2]]
            if (read_nuc_seg_map) {
                x@raw@nuc_seg_map <- maps[[1]]
            }


        } else{
            if (length(grep('_memb_seg_map.tif', img_names)) > 0) {
                x@raw@mem_seg_map <- readTIFF(file.path(sample_dir,
                                                        img_names[grep('_memb_seg_map.tif', img_names)]))
            }

            if (read_nuc_seg_map &&
                length(grep('_nuc_seg_map.tif', img_names)) > 0) {
                x@raw@nuc_seg_map <- readTIFF(file.path(sample_dir,
                                                        img_names[grep('_nuc_seg_map.tif', img_names)]))
            }
        }
        if (read_dapi_map &&
            length(grep('_component_data.tif', img_names)) > 0) {
            dapi <- readTIFF(file.path(sample_dir,
                                       img_names[grep('_component_data.tif', img_names)]), all = TRUE)[[1]]
            if (class(dapi) == 'array') {
                dapi <- dapi[, , 4]
            }
            x@raw@dapi_map <- dapi
        }

        if (!ignore_scoring) {
            score_data <- img_names[grep('_score_data.txt$', img_names)]
            if (length(score_data) != 1) {
                stop(
                    'Could not find a single *_score_data.txt for ',
                    x@coordinate_name,
                    ' in ',
                    sample_dir,
                    'if no markers were scored please set ignore_scoring flag to TRUE.'
                )
            }
            x@raw@score <-
                t(read.csv(
                    file.path(sample_dir, score_data),
                    sep = '\t',
                    as.is = TRUE
                ))
        }

        #in most cases we need only a fraction of columns in the segmentation file
        if (read_only_relevant_markers){
            markers <- c('Phenotype',
                         'Cell.ID',
                         'Cell.X.Position',
                         'Cell.Y.Position',
                         'Nucleus.Area..pixels.',
                         'Nucleus.Compactness',
                         'Nucleus.Minor.Axis',
                         'Nucleus.Major.Axis',
                         "Membrane.Area..pixels.",
                         "Membrane.Compactness",
                         "Membrane.Minor.Axis",
                         "Membrane.Major.Axis",
                         "Entire.Cell.Area..pixels.",
                         "Entire.Cell.Minor.Axis",
                         "Entire.Cell.Major.Axis",
                         "Confidence")

            if (!is.null(x@raw@score)){
                #extract all relevant markers
                scores <- x@raw@score
                scores <- scores[grep('Component',rownames(scores)),]
                scores <- gsub('[ ()-]','.',scores)
                scores <- colnames(x@raw@data)[unlist(lapply(scores,
                                                             function(x,y)
                                                                 grep(x,y),colnames(x@raw@data)))]
                scores <- scores[grep('Mean',scores)]
                markers <- c(markers, scores)
            }
            x@raw@data <- x@raw@data[,markers]
        }

        #fix the labels if necessary
        if (length(label_fix) > 0) {
            #for each label fix
            for (fix in label_fix) {
                x@raw@data$Phenotype[x@raw@data$Phenotype==fix[1]] <- fix[2]
            }
        }

        if (length(x@raw@data$memb_seg_map) > 0) {
            x_max <- ncol(x@raw@memb_seg_map)
            y_max <- nrow(x@raw@memb_seg_map)
        } else{
            x_max <- max(x@raw@data$Cell.X.Position)
            y_max <- max(x@raw@data$Cell.Y.Position)
        }

        #estimate the window size
        window = owin(xrange = c(0, x_max),
                      yrange = c(0, y_max))

        #sqeeze the data into the spatstats package format
        x@ppp  <- with(
            x@raw@data,
            ppp(
                `Cell.X.Position`,
                `Cell.Y.Position`,
                window = window,
                marks = factor(x@raw@data$Phenotype)
            )
        )

        if (readMasks) {
            #extract mask data
            x <-
                extract_mask_data(x,
                                  img_names,
                                  sample_dir,
                                  x@coordinate_name,
                                  invasive_margin_in_px)
        }

        if (length(x@mask) > 0) {
            di <- dim(x@mask[[1]])
            x@size_in_px <- di[1] * di[2]
        } else{
            x@size_in_px <- x_max * y_max
        }

        if (ROI) {
            #extract ROI mask
            roi <- img_names[grep('_ROI.tif', img_names)]
            if (length(roi) == 0) {
                stop(
                    paste(
                        '_ROI.tif for',
                        x@coordinate_name,
                        'in',
                        sample_dir,
                        'does not exists!'
                    )
                )
            }
            roi <- file.path(sample_dir, roi)
            roi <- extract_mask(roi)

            #drop all coordinates outside of the image
            filter <-
                sapply(1:length(x@ppp$x), function(i, dat, mask)
                    mask[dat$x[i], dat$y[i]] == 1, x@ppp, roi)
            x@ppp <- x@ppp[filter, ]
            x@raw@data <- x@raw@data[filter, ]

            #change the size of the image
            x@size_in_px <- sum(roi > 0)

            #if there were other masks read we set all these masks to 0
            if (readMasks) {
                #reduce the other masks
                for (i in 1:length(x@masks)) {
                    x@masks[[i]] <- x@masks[[i]][roi == 0] <- 0
                }

            }
            x@mask$ROI <- roi

        }
        return(x)
    }
)

#' Read inForm output from a single coordinate
#'
#' @param filename Name of the .tif file that contains the mask.
#'
#' @return Mask matrix
#' @export
#' @examples
#' extract_mask(system.file("extdata",
#'                          "MEL2","MEL2_080416_2_Invasive_Margin.tif",
#'                          package = "IrisSpatialFeatures"))
extract_mask <- function(filename) {
    mask <- readTIFF(filename)
    mask <- as.matrix((mask[, , 1] + mask[, , 2] + mask[, , 3]) > 0)
    mask <- t(mask)
    return(mask)
}


#' @useDynLib IrisSpatialFeatures
#' @importFrom Rcpp sourceCpp
setGeneric("extract_mask_data", function(x, ...)
    standardGeneric("extract_mask_data"))
setMethod(
    "extract_mask_data",
    signature = "Coordinate",
    definition = function(x,
                          img_names,
                          sample_dir,
                          coordinate_name,
                          invasive_margin_in_px) {
        #invasive margin
        inv_mar <-
            img_names[grep('_Invasive.Margin.tif', img_names)]
        if (length(inv_mar) == 0) {
            stop(
                paste(
                    '_Invasive_Margin.tif for',
                    coordinate_name,
                    'in',
                    sample_dir,
                    'does not exists!'
                )
            )
        }
        inv_mar <- file.path(sample_dir, inv_mar)
        x@mask$margin_line <- extract_mask(inv_mar)
        if (sum(x@mask$margin_line) == 0) {
            stop(
                paste(
                    '_Invasive.Margin.tif for',
                    coordinate_name,
                    'in',
                    sample_dir,
                    'has no pixels, please check the mask file!'
                )
            )
        }
        x@mask$invasive_margin <-
            growMarginC(x@mask$margin_line, invasive_margin_in_px)

        #tumor mask
        tumor_tif <- img_names[grep('_Tumor.tif', img_names)]
        if (length(tumor_tif) == 0) {
            stop(
                paste(
                    '_Tumor.tif for',
                    coordinate_name,
                    'in',
                    sample_dir,
                    'does not exists!'
                )
            )
        }
        tumor_tif <- file.path(sample_dir, tumor_tif)
        x@mask$tumor <- extract_mask(tumor_tif)
        if (sum(x@mask$margin_line) == 0) {
            stop(
                paste(
                    '_Tumor.tif for',
                    coordinate_name,
                    'in',
                    sample_dir,
                    'has no pixels, please check the mask file!'
                )
            )
        }
        if (!all(dim(x@mask$invasive_margin) == dim(x@mask$tumor))) {
            stop(
                paste(
                    '_Tumor.tif for',
                    coordinate_name,
                    'in',
                    sample_dir,
                    'does not have the same dimension as the invasive margin mask!'
                )
            )
        }
        x@mask$tumor[x@mask$invasive_margin > 0] <- 0

        #also add the stroma
        x@mask$stroma <- x@mask$tumor
        x@mask$stroma[x@mask$stroma == 0] <- 1
        x@mask$stroma[x@mask$tumor > 0] <- 0
        x@mask$stroma[x@mask$invasive_margin > 0] <- 0

        return(x)
    }
)

