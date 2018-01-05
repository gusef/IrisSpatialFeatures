#' Plot all coordinates in a given dataset
#'
#' @param x Sample ImageSet object of the IrisSpatialFeatures package.
#' @param outdir Output directory (default: './')
#' @param palette Color palette used for the different cell-types. (default: NULL)
#' @param type File format for the plots. Can bei either 'pdf' or 'png'. (default: 'pdf')
#' @param width Width of the plot in inches for pdf and pixels for png. (default: 10)
#' @param height Heigth of the plot in inches for pdf and pixels for png. (default: 7)
#' @param ... Additional arguments
#' @return plot of all coordinates
#'
#' @docType methods
#' @export
#'
#' @examples
#'
#' #loading pre-read dataset
#' dataset <- IrisSpatialFeatures_data
#' overview_plot(dataset)
#'
#' @rdname overview_plot
setGeneric("overview_plot", function(x, ...)
    standardGeneric("overview_plot"))

#' @rdname overview_plot
#' @aliases overview_plot,ANY,ANY-method
setMethod(
    "overview_plot",
    signature = "ImageSet",
    definition = function(x,
                          outdir = './',
                          palette = NULL,
                          type = 'pdf',
                          width = 10,
                          height = 7) {
        lapply(
            x@samples,
            overview_plot_sample,
            all_levels = x@markers,
            outdir,
            palette,
            type,
            width,
            height
        )
    }
)

setGeneric("overview_plot_sample", function(x, ...)
    standardGeneric("overview_plot_sample"))
setMethod(
    "overview_plot_sample",
    signature = "Sample",
    definition = function(x,
                          all_levels,
                          outdir,
                          palette,
                          type,
                          width,
                          height) {
        lapply(
            x@coordinates,
            overview_plot_coord,
            sample_name = x@sample_name,
            all_levels,
            outdir,
            palette,
            type,
            width,
            height
        )
    }
)

#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics layout
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom RColorBrewer brewer.pal
setGeneric("overview_plot_coord", function(x, ...)
    standardGeneric("overview_plot_coord"))
setMethod(
    "overview_plot_coord",
    signature = "Coordinate",
    definition = function(x,
                          sample_name,
                          all_levels,
                          outdir,
                          palette,
                          type,
                          width,
                          height) {
        if (is.null(palette)) {
            palette <- brewer.pal(length(all_levels), 'Spectral')
        }
        mapping <- data.frame(col = palette, lvl = all_levels)

        #preprocess data
        df <- data.frame(x@ppp)
        df$cols <- mapping$col[match(df$marks, mapping$lvl)]

        #build the filename
        file_stub <- paste0(sample_name, '_', x@coordinate_name)

        if (type == 'pdf') {
            pdf(
                file = file.path(outdir, paste0(file_stub, '.pdf')),
                width = width,
                height = height
            )
        } else if (type == 'png') {
            png(file.path(outdir, paste0(file_stub, '.png')),
                width = 800,
                height = 600)
        }
        par(mar = c(4, 4, 4, 0))
        layout(matrix(c(1, 2), ncol = 2),
               widths = c(0.8, 0.25),
               heights = 1)
        plot(
            df$x,
            df$y,
            col = as.character(df$cols),
            pch = 18,
            ylab = 'y',
            xlab = 'x',
            main = paste(sample_name, '-', x@coordinate_name)
        )
        par(mar = c(1, 1, 1, 1))
        image(matrix(c(1, 1, 1, 1), ncol = 2), col = 'white', axes = FALSE)
        legend(
            'left',
            col = palette,
            legend = all_levels,
            pch = 18,
            cex = 0.8
        )
        dev.off()
    }
)

#' Plot all frames of an IrisSpatialFeatures ImageSet object.
#'
#' @param x Iris ImageSet object
#'
#' @return A plot
#' @examples
#' dataset <- IrisSpatialFeatures_data
#' plot(dataset)
#' @docType methods
#' @export
#' @rdname plot
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @import RColorBrewer
#' @import tibble
#'
setMethod("plot",
    signature = c(x="ImageSet"),
    function(x) {
        tbl <- as_tibble(as.data.frame(x))
        sample <- NA
        frame_index <- NA
        lownum <- NA
        frame_index <- NA
        marks <- NA
        frame <- NA
        rowid <- NA
        numbered <- tbl %>% group_by(sample,frame) %>% slice(1) %>% ungroup %>% arrange(sample,frame) %>% rowid_to_column()
        low_ids <- numbered %>% group_by(sample) %>% filter(rowid==min(rowid)) %>% ungroup %>% rename(lownum=rowid) %>% select(lownum,sample)
        numbered_frames <- numbered %>% left_join(low_ids,by="sample") %>% mutate(frame_index=rowid-lownum) %>% select(frame,sample,frame_index)
        idxtbl <- tbl %>% left_join(numbered_frames,by=c("sample","frame"))
        v <- ggplot(data=idxtbl,mapping=aes_string(x='x',y='y'))+
            geom_point(aes_string(color='marks'),alpha=0.8,stroke=0)+
            facet_grid(sample~frame_index)+
            theme_minimal()+
            scale_color_brewer(palette="Dark2")+
            theme(strip.text.x = element_blank())+
            geom_label(data = numbered_frames,
                       mapping=aes(x=0,y=0,label =frame),
                       inherit.aes=FALSE,
                       hjust=0,
                       vjust=0)+
            coord_fixed(ratio = 1)
        return(v)
})
