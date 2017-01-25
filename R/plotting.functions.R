
setGeneric("overview.plot", function(object, ...) standardGeneric("overview.plot"))
setMethod("overview.plot",
          signature = "Iris",
          definition = function(object,outdir='./',palette=NULL,type='pdf',width=10,height=7){
          lapply(object@samples,
                 overview.plot.sample,
                 all_levels=object@markers,
                 outdir,
                 palette,
                 type,
                 width,
                 height)
})

setGeneric("overview.plot.sample", function(object, ...) standardGeneric("overview.plot.sample"))
setMethod("overview.plot.sample",
          signature = "Sample",
          definition = function(object,all_levels,outdir,palette,type,width,height){
          lapply(object@coordinates,
                 overview.plot.coord,
                 sample_name=object@sample_name,
                 all_levels,
                 outdir,
                 palette,
                 type,
                 width,
                 height)        
})

setGeneric("overview.plot.coord", function(object, ...) standardGeneric("overview.plot.coord"))
setMethod("overview.plot.coord",
          signature = "Coordinate",
          definition = function(object,sample_name,all_levels,outdir,palette,type,width,height){
  
     if (is.null(palette)){
         palette <- brewer.pal(length(all_levels),'Spectral')
     }
     mapping <- data.frame(col=palette,lvl=all_levels)
     
     #preprocess data
     df <- data.frame(object@ppp)
     df$cols <- mapping$col[match(df$marks,mapping$lvl)]
     
     #build the filename
     file_stub <- paste0(sample_name,'_',object@coordinate_name)
            
     if (type == 'pdf'){
         pdf(file = file.path(outdir,paste0(file_stub,'.pdf')),width=width,height=height)
     }else if(type == 'png'){
         png(file.path(outdir,paste0(file_stub,'.png')),width = 800, height=600)
     }
     par(mar=c(4,4,4,0))
     layout(matrix(c(1,2), ncol = 2), widths = c(0.8,0.25), heights = 1)
     plot(df$x,
          df$y,
          col=as.character(df$cols),
          pch=18,
          ylab='y',
          xlab='x',
          main=paste(sample_name,'-',object@coordinate_name))             
     par(mar=c(1,1,1,1))
     image(matrix(c(1,1,1,1),ncol=2),col='white',axes=F)
     legend('left',col = palette,legend = all_levels,pch=18,cex = 0.8)
     dev.off()
})





