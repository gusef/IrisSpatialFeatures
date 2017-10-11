##############################################################################

#' Function to extract all numeric features
#' @param dat A data matrix with features as rows and samples as columns
#' @param lab Label annotation that contains 2 classes,
#'            which corresponds to the samples in the column
#'
#' @return t-test and wilcox test btween the 2 classes
#'
#'
#' @examples
#' dat <- cbind(matrix(runif(400),ncol=10),matrix(runif(400)+0.2,ncol=10))
#' lab <- c(rep('classA',10),rep('classB',10))
#' rownames(dat) <- paste0('F',1:nrow(dat))
#' feature_selection(dat,lab)
#'
#' @export
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' @importFrom stats p.adjust
feature_selection <- function(dat, lab) {
    if (length(unique(lab)) != 2) {
        stop('the variable lab must have 2 classes')
    }
    res <- sapply(1:nrow(dat),
                  function(x, dat, lab)
                      t.test(dat[x, lab == unique(lab)[1]],
                             dat[x, lab == unique(lab)[2]])$p.value,
                  dat,
                  lab)
    ttest_res <- data.frame(
        Feature = rownames(dat),
        p.value = res,
        adj.p.Val = p.adjust(res,
                             method = 'BH',
                             n = length(res))
    )
    ttest_res <- ttest_res[order(ttest_res$p.value), ]

    #using a wilcoxon test
    res <- sapply(1:nrow(dat),
                  function(x, dat, lab)
                      wilcox.test(dat[x, lab == unique(lab)[1]],
                                  dat[x, lab == unique(lab)[2]])$p.value,
                  dat,
                  lab)
    wilcox_res <- data.frame(
        Feature = rownames(dat),
        p.value = res,
        adj.p.Val = p.adjust(res,
                             method = 'BH',
                             n = length(res))
    )
    wilcox_res <- wilcox_res[order(wilcox_res$p.value), ]
    rownames(wilcox_res) <- NULL
    return(list(t_test = ttest_res,
                wilcox = wilcox_res))
}

#' Extract all spatial features
#' @param x IrisSpatialFeatures ImageSet object
#' @param name Prefix for all features, e.g. 'invasive_margin' (Default: '')
#' @param rm.na Should features with NA values be removed (Default: FALSE)
#' @param ... Additional arguments
#' @return dataframe of features
#' @examples
#'
#' #' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' dataset <- extract_nearest_neighbor(dataset,min_num_cells=2)
#' dataset <- extract_proximity(dataset,only_closest=TRUE,radii=25)
#' dataset <- extract_interactions(dataset)
#' extract_features(dataset)
#' @docType methods
#' @export
#' @rdname extract_features
setGeneric("extract_features", function(x, ...)
    standardGeneric("extract_features"))

#' @rdname extract_features
#' @aliases extract_features,ANY,ANY-method
setMethod(
    "extract_features",
    signature = "ImageSet",
    definition = function(x, name = '', rm.na = FALSE) {
        counts <- get_counts_per_mm2_noncollapsed(x)
        counts <- sapply(counts, colMeans)

        rownames(counts) <-
            paste('Counts -', rownames(counts))
        count_ratios <- extractRatios(counts, 'Counts')
        dat <- rbind(counts,
                     count_ratios)

        if (length(x@interactions) > 0) {
            f_inter <-
                lapply(x@interactions,
                       extract_interaction_features,
                       'Interaction')
            f_inter <- do.call(rbind, f_inter)
            dat <- rbind(dat, t(f_nn))
        } else{
            message(paste('Skipping interactions .. ',
                   'please run extract_interactions to include them.',
                   sep=""))
        }

        if (length(x@nearest_neighbors) > 0) {
            f_nn <-
                lapply(x@nearest_neighbors,
                       extract_interaction_features,
                       'NN')
            f_nn <- do.call(rbind, f_nn)
            dat <- rbind(dat, t(f_nn))
        } else{
            message(paste(
            'Skipping nearest neighbors .. ',
            'please run extract_nearest_neighbor to include them.',
            sep=""))
        }
        dat <- dat[!duplicated(rownames(dat)), ]
        #some of the ratios cause infinte values
        dat[is.infinite(dat)] <- NA
        if (rm.na) {
            dat <- dat[rowSums(is.na(dat)) == 0, ]
        }
        return(dat)
    }
)


#extracts the values for NN and interaction analysis
extractSimpleValues <- function(mat, remove_self = TRUE) {
    big_mat <- array(unlist(mat), dim = c(dim(mat[[1]]), length(mat)))
    phenos <- colnames(mat[[1]])
    combinations <-
        expand.grid(seq(length(phenos)), seq(length(phenos)))
    combinations <- as.matrix(combinations)
    combinations <-
        rbind(combinations, cbind(combinations[, 2], combinations[, 1]))
    if (remove_self) {
        #remove combinations where both values are the sames
        #(nearest neighbor with itself doesn't make sense)
        combinations <-
            combinations[combinations[, 1] != combinations[, 2], ]
    }
    collapsed <-
        t(apply(combinations, 1, function(x, bm)
            bm[x[2], x[1], ], big_mat))
    combinations <- apply(combinations, 2, function(x, p)
        p[x], phenos)
    rownames(collapsed) <-
        paste(combinations[, 1], ' -> ', combinations[, 2])
    colnames(collapsed) <- names(mat)
    collapsed[is.nan(collapsed)] <- 0

    #remove all the duplicates
    collapsed <- collapsed[!duplicated(rownames(collapsed)),]

    return(collapsed)
}

getPaired <- function(nams) {
    tab <- table(nams)
    tab <- names(tab)[tab == 2]
    return(nams %in% tab)
}

extractRatios <- function(mat, nam) {
    #if we look at interactions or nearest neighbors
    if (length(grep(' -> ', rownames(mat), fixed = TRUE)) == nrow(mat)) {
        nams <-
            t(sapply(strsplit(
                sapply(strsplit(rownames(mat), ' - '), function(x)
                    x[2]), ' -> '
            ), function(x)
                x))
        nams[, 1] <- sub('[+-] $', '', nams[, 1])
        paired <- getPaired(paste(nams[, 1], nams[, 2]))
        nams <- nams[paired, ]
        COUNTS <- FALSE
    } else{
        COUNTS <- TRUE
        nams <-
            sub('[+-]$', '', sapply(strsplit(rownames(mat), ' - '),
                function(x)
                x[2]))
        paired <- getPaired(nams)
        nams <- nams[paired]
    }
    if (sum(paired) == 0) {
        ratios <- matrix(nrow = 0, ncol = ncol(mat))
        colnames(ratios) <- colnames(mat)
    } else{
        mat <- mat[paired, ]

        #get the ratios
        num_pairs <- seq(nrow(mat) / 2)
        indices <- sort(rep(num_pairs, 2))
        ratios <- t(sapply(num_pairs,
                function(x, indices, mat)
                    mat[grep(x, indices)[1], ] / mat[grep(x, indices)[2], ],
                        indices,
                        mat))
        #log2 to get a nicer behavior
        ratios <- log2(ratios)
        if (COUNTS) {
            rownames(ratios) <- paste('ratio -', unique(paste0(nams, '+/-')))
        } else{
            rownames(ratios) <-
                paste('ratio -', unique(paste0(nams[, 1], '+/- -> ',
                    nams[, 2])))
        }
        rownames(ratios) <- paste(nam, '-', rownames(ratios))
    }
    return(ratios)
}

extract_interaction_features <- function(interactions, nam) {
    f_interactions <- extractSimpleValues(interactions)
    rownames(f_interactions) <-
        paste(nam, '-', rownames(f_interactions))
    f_int_ratios <- extractRatios(mat = f_interactions, nam)
    dat <- rbind(f_interactions,
                 f_int_ratios)
    return(dat[,'means'])
}

extract_count_features <- function(f_counts, nam) {
    rownames(f_counts) <- paste(nam, '-', rownames(f_counts))
    f_count_ratios <- extractRatios(mat = f_counts, nam)
    dat <- rbind(f_counts, f_count_ratios)
    return(dat)
}

