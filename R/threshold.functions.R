###########################################
######### Thresholding functions

#' This function reads the manually determined thresholds of certain markers (e.g. PD1, PD-L1) and splits selected celltypes into marker+ and marker- celltypes.
#'
#' @param image_set IrisSpatialFeatures ImageSet object.
#' @param marker Name of the marker used in the score file.
#' @param marker_name corresponding name, which should be appended at the selected cell types.
#' @param base Vector of cell types for which the marker should be used.
#' @param pheno_name Name of the phenotype column to be used. (Default from inForm is "Phenotype")
#' @param remove_blanks Flag that indicates whether or not not called cells are to be removed. (Default: TRUE)
#'
#' @docType methods
#' @return IrisSpatialFeatures ImageSet object.
#' @examples
#' dataset <- read_raw(path=system.file("extdata", package = "IrisSpatialFeatures"),
#'                      format='Mantra')
#' dataset <- threshold_dataset(dataset,
#'                              marker='PD-Ligand-1 (Opal 690)',
#'                              marker_name='PDL1',
#'                              base=c('SOX10+'))
#' dataset <- threshold_dataset(dataset,
#'                              marker='PD-1 (Opal 540)',
#'                              marker_name='PD1',
#'                              base=c('CD8+','OTHER'))
#' @export
#' @rdname threshold_dataset
setGeneric("threshold_dataset",
           function(image_set, marker, marker_name, base=NULL, pheno_name='Phenotype', remove_blanks=TRUE) {
               standardGeneric("threshold_dataset")
           },
    valueClass = "ImageSet")

#' @rdname threshold_dataset
#' @aliases threshold_dataset,ANY,ANY-method
setMethod(
    "threshold_dataset",
    signature = "ImageSet",
    definition = function (image_set,
                           marker,
                           marker_name,
                           base = NULL,
                           pheno_name = 'Phenotype',
                           remove_blanks = TRUE) {
        x <- image_set
        #for each sample
        x@samples <-
            lapply(
                x@samples,
                threshold_samples,
                marker,
                marker_name,
                base,
                pheno_name,
                remove_blanks
            )
        names(x@samples) <-
            sapply(x@samples, function(x)
                x@sample_name)
        x <- extract_counts(x)
        x@nearest_neighbors <- list()
        x@interactions <- list()

        return(x)
    }
)

setGeneric("threshold_samples", function(x, ...)
    standardGeneric("threshold_samples"))
setMethod(
    "threshold_samples",
    signature = "Sample",
    definition = function (x,
                           marker,
                           marker_name,
                           base,
                           pheno_name,
                           remove_blanks) {
        #for each coordinate
        x@coordinates <-
            lapply(
                x@coordinates,
                threshold_coords,
                marker,
                marker_name,
                base,
                pheno_name,
                x@sample_name,
                remove_blanks
            )
        names(x@coordinates) <-
            sapply(x@coordinates, function(x)
                x@coordinate_name)
        return(x)
    }
)


setGeneric("threshold_coords", function(x, ...)
    standardGeneric("threshold_coords"))
setMethod(
    "threshold_coords",
    signature = "Coordinate",
    definition = function (x,
                           marker,
                           marker_name,
                           base,
                           pheno_name,
                           sample_name,
                           remove_blanks) {

        # Create a list of phenotypes we will be using after this
        newlist = list()
        i <- 0
        for (phenotype in levels(x@ppp$marks)) {
            i <- i + 1
            if (phenotype %in% base) {
                newlist[[i]] <- paste0(phenotype,' ',marker_name,'+')
                i <- i+1
                newlist[[i]] <- paste0(phenotype,' ',marker_name,'-')
            } else {
                newlist[[i]] <- phenotype
            }
        }
        new_phenotypes <- unlist(newlist)

        #remove the cells that are not called by inForm (usually not very many!)
        if (remove_blanks) {
            x@raw@data <- x@raw@data[x@raw@data[[pheno_name]] != '', ]
            x@ppp <- x@ppp[x@ppp$marks != '', ]
            x@ppp$marks <- droplevels(x@ppp$marks)
        }

        #if no cell types were specified
        if (is.null(base)) {
            base <- levels(x@ppp$marks)
        }

        #make a combined phenotype column in the rawdata
        if (!paste0(pheno_name, '.combined') %in% colnames(x@raw@data)) {
            x@raw@data[[paste0(pheno_name, '.combined')]] <-
                x@raw@data[[pheno_name]]
        }

        #get the thresholds for the marker we want to score
        scoring <- getScoring(x)
        #r automatically replaces special characters with '.'
        #when they are used for names so I'm fixing this
        mark <- gsub('[ \\(\\)]', '.', marker)
        mark <- gsub('-', '.', mark)

        if (sum(scoring$Component == mark) == 0) {
            stop(
                'Could not find Score for: ',
                marker,
                ' for Sample: ',
                sample_name,
                ' Coordinate: ',
                x@coordinate_name
            )
        }
        scoring <- scoring[scoring$Component == mark, ]

        #extract the current marker expression
        expression <- x@raw@data[[paste(
            scoring$Compartment,
            scoring$Component,
            'Mean..Normalized.Counts..Total.Weighting.',
            sep = '.'
        )]]


        #sometimes there are #N/A from Inform
        if (class(expression) == 'character') {
            keep <- expression != '#N/A'
            x@raw@data <- x@raw@data[keep, ]
            x@ppp <- x@ppp[keep, ]
            expression <- as.numeric(expression[keep])
        }

        expression <- as.numeric(expression)

        #deterine the positive cells
        positive_cells <- expression > scoring$Threshold

        #use only the cells that are within base
        current <- x@raw@data[[pheno_name]] %in% base
        current_base <- x@raw@data$Phenotype.combined[current]

        #fetch the cells we are currently working on
        x@raw@data$Phenotype.combined[current &
                                          !positive_cells] <-
            paste0(current_base[!positive_cells[current]], ' ', marker_name, '-')
        x@raw@data$Phenotype.combined[current &
                                          positive_cells] <-
            paste0(current_base[positive_cells[current]], ' ', marker_name, '+')
        x@ppp$marks <-
            as.factor(x@raw@data$Phenotype.combined)

        # fix our levels if they are explicitly defined
        x@ppp$marks <- factor(x@ppp$marks,levels=new_phenotypes)

        return(x)
    }
)

setGeneric("getScoring", function(x, ...)
    standardGeneric("getScoring"))
setMethod(
    "getScoring",
    signature = "Coordinate",
    definition = function(x) {
        scoring <- x@raw@score
        scores <- matrix(nrow = 0, ncol = 3)
        colnames(scores) <-
            c('Compartment', 'Component', 'Threshold')
        #if there were more than one additional markers scores:
        if (length(grep('First', rownames(scoring)) > 0)) {
            tab <- c('First', 'Second', 'Third')
            for (i in seq(length(tab))) {
                #only if indicator actually exists (as of now I'm not sure if inForm allows for more than 2 markers)
                if (length(grep(tab[i], rownames(scoring)) > 0)) {
                    compartment <- scoring[paste0(tab[i], '.Cell.Compartment'), 1]
                    component <-
                        scoring[paste0(tab[i], '.Stain.Component'), 1]
                    #r automatically replaces special characters with '.' when they are used for names so I'm fixing this
                    component <-
                        gsub('[ \\(\\)]', '.', component)
                    component <- gsub('-', '.', component)
                    threshold <-
                        scoring[paste0(component, '.Threshold'), 1]
                    scores <-
                        rbind(scores, c(compartment, component, threshold))
                }
            }
            #if there was only one marker scored
        } else{
            compartment <- scoring['Cell.Compartment', 1]
            component <- scoring['Stain.Component', 1]
            #r automatically replaces special characters with '.' when they are used for names so I'm fixing this
            component <- gsub('[ \\(\\)]', '.', component)
            component <- gsub('-', '.', component)
            threshold <- scoring['Positivity.Threshold', 1]
            scores <-
                rbind(scores, c(compartment, component, threshold))
        }
        scores <- data.frame(scores, stringsAsFactors = FALSE)
        scores$Threshold <- as.numeric(scores$Threshold)
        return(scores)
    }
)

###########################################
######### Collapse functions

#' This function collapses two markers into one, and reruns the counting of cells.
#' Mostly a convenience function for the Shiny interface so we start with a completely split set
#' successively adding more markers
#'
#' @param image_set IrisSpatialFeatures ImageSet object.
#' @param marker1 Name of the first marker that should be collapsed.
#' @param marker2 Name of the second marker that should be collapsed.
#' @param combined Name of the combined marker.
#'
#' @docType methods
#' @return IrisSpatialFeatures ImageSet object.
#' @examples
#' dataset <- IrisSpatialFeatures_data
#' ds <- collapse_markers(dataset,marker1="SOX10+ PDL1+",marker2="SOX10+ PDL1-", combined="SOX10+")
#'
#' @export
#' @importFrom methods .valueClassTest
#' @rdname collapse_markers
setGeneric("collapse_markers",
           function(image_set, marker1, marker2, combined) {
               standardGeneric("collapse_markers")
           },
           valueClass = "ImageSet")

#' @rdname collapse_markers
#' @aliases collapse_markers,ANY,ANY-method
setMethod(
    "collapse_markers",
    signature = "ImageSet",
    definition = function (image_set,
                           marker1,
                           marker2,
                           combined) {
        x <- image_set

        #check if the markers are actually in the Iris object
        if (!all(c(marker1, marker2) %in% x@markers)){
            stop('At least one of the specified markers is not included in the dataset')
        }

        #collapse the marker
        for (idx in 1:length(x@samples)){
            for (jdx in 1:length(x@samples[[idx]]@coordinates)){

                #fix the labeling
                pheno <- x@samples[[idx]]@coordinates[[jdx]]@raw@data$Phenotype
                pheno[pheno %in% c(marker1, marker2)] <- combined

                #add it to the ppp and raw data object
                x@samples[[idx]]@coordinates[[jdx]]@raw@data$Phenotype <- pheno
                x@samples[[idx]]@coordinates[[jdx]]@ppp$marks <- as.factor(pheno)
            }
        }

        #fix the markers in the object
        x@markers <- sort(c(x@markers[!x@markers %in% c(marker1, marker2)], combined))
        x <- extract_counts(x)
        x@nearest_neighbors <- list()
        x@interactions <- list()

        return(x)
    }
)


##############################################################################
# ROI Functions

#' Method that reduces the current dataset to a specific region of interest, discarding all cell coordinates outside of that region
#'
#' @param x IrisSpatialFeatures ImageSet object
#' @param ROI Region of interest (default: 'invasive_margin')
#' @param verbose Print some details (default: False)
#' @param ... Additional arguments
#'
#' @return IrisSpatialFeatures ImageSet object
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' im_area <- extract_ROI(dataset,ROI='invasive_margin')
#'
#' @docType methods
#' @export
#' @rdname extract_ROI
setGeneric("extract_ROI", function(x, ...)
    standardGeneric("extract_ROI"))

#' @rdname extract_ROI
#' @aliases extract_ROI,ANY,ANY-method
setMethod(
    "extract_ROI",
    signature = "ImageSet",
    definition = function(x, ROI = 'invasive_margin', verbose = FALSE) {
        if (length(x@samples[[1]]@coordinates[[1]]@mask[[ROI]]) == 0) {
            stop('There is no mask for "', ROI, '"')
        }

        x@samples <- lapply(x@samples, extract_ROI_sample, ROI, verbose)

        #update the counts
        x <- extract_counts(x)

        #reset all spatial stats
        x@nearest_neighbors <- list()
        x@interactions <- list()
        x@proximity <- list()

        return(x)
    }
)

setGeneric("extract_ROI_sample", function(x, ...)
    standardGeneric("extract_ROI_sample"))
setMethod(
    "extract_ROI_sample",
    signature = "Sample",
    definition = function(x, ROI, verbose) {
        x@coordinates <- lapply(x@coordinates, extract_ROI_Coordinate, ROI, verbose)
        return(x)
    }
)


setGeneric("extract_ROI_Coordinate", function(x, ...)
    standardGeneric("extract_ROI_Coordinate"))
setMethod(
    "extract_ROI_Coordinate",
    signature = "Coordinate",
    definition = function(x, ROI, verbose) {
        #reduce to the filter
        mask <- x@mask[[ROI]]
        filter <-
            sapply(1:length(x@ppp$x), function(i, dat, mask)
                mask[dat$x[i], dat$y[i]] == 1, x@ppp, mask)
        x@ppp <- x@ppp[filter, ]
        x@raw@data <- x@raw@data[filter, ]
        if (verbose) { print(paste0("Old pixel size ",x@size_in_px)) }
        x@size_in_px <- sum(mask > 0)
        if (verbose)  { print(paste0("New pixel size ",x@size_in_px)) }
        return(x)
    }
)
