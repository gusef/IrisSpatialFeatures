
# tools::vignetteEngine("knitr", weave = vweave, tangle = vtangle,
#                       pattern = "[.]Rmd$", package = "knitr")

setwd('c:/work/projects/imaging/Iris/test/')
rm(list=ls())
gc()
require(Iris)


raw_data <- Iris()
raw_data<- read.raw(raw_data,
                    raw_dir_name='../inst/extdata/',
                    format='Mantra')

#apply all the thresholds PD1 for T and other cells, PD-L1 for macrophages and tumor cells
dataset <- threshold.dataset(raw_data,
                             marker='PD-Ligand-1 (Opal 690)',
                             marker_name='PDL1',
                             base=c('SOX10+'))
dataset <- threshold.dataset(dataset,
                             marker='PD-1 (Opal 540)',
                             marker_name='PD1',
                             base=c('CD8+','OTHER'))

#generate overview plots
p <- overview.plot(dataset,outdir='overview/',palette=NULL,type='pdf')
p <- overview.plot(dataset,outdir='overview/',palette=NULL,type='png')

#get the counts
get.counts.per.mm2(dataset)
get.counts.per.mm2.noncollapsed(dataset)
get.count.ratios(dataset,'SOX10+ PDL1-','SOX10+ PDL1+')


#run the nearest neighbor analysis
dataset <- extract.nearest.neighbor(dataset,min_num_cells=2)
get.nearest.neighbors(dataset,"SOX10+ PDL1+")
p <- plot.nearest.neighbor(dataset,'CD8+ PD1+','SOX10+ PDL1')
p <- plot.nearest.neighbor(dataset,'SOX10+ PDL1+','SOX10+ PDL1-')

dataset <- extract.nearest.neighbor(dataset,min_num_cells=2)
get.nearest.neighbors(dataset,"SOX10+ PDL1+")
p <- plot.nearest.neighbor(dataset,'CD8+ PD1+','SOX10+ PDL1')
p <- plot.nearest.neighbor(dataset,'SOX10+ PDL1+','SOX10+ PDL1-')

#ray plots for 
neighbor.ray.plot(dataset,from_type='SOX10+ PDL1-',to_type='OTHER PD1-',plot_dir='plots/')

                  
#run the interaction analysis
dataset <- extract.interactions(dataset)
get.interactions(dataset,'CD8+ PD1+')

#plotting interaction summaries
p <- plot.interactions(dataset,'CD8+ PD1+',xlim_fix=4)

#plotting interaction maps
int_markers <- c('CD8+ PD1+','SOX10+ PDL1+')
int_marker_cols <- c('#dd1c77','#99d8c9')
silent_markers <- c('CD8+ PD1-')
silent_col=c('yellow')
p <- interaction.maps(dataset,int_markers,int_marker_cols,silent_markers,silent_col)

#running the proximity / touching analysis
dataset <- extract.proximity(dataset)
p <- plot.proximities(dataset,"SOX10+ PDL1-",xlim_fix=4)

#get interactions again
tumor_area <- extract.ROI(dataset,ROI='tumor')
get.counts.per.mm2(tumor_area)
tumor_area <- extract.interactions(tumor_area)

p <- plot.interactions(tumor_area,"SOX10+ PDL1-",xlim_fix=4)

#get counts in all vs different areas
get.counts.per.mm2(dataset)
get.counts.per.mm2.noncollapsed(dataset)
get.counts.per.mm2(extract.ROI(dataset,ROI='tumor'))
get.counts.per.mm2.noncollapsed(extract.ROI(dataset,ROI='tumor'))
get.counts.per.mm2(extract.ROI(dataset,ROI='invasive_margin'))
get.counts.per.mm2.noncollapsed(extract.ROI(dataset,ROI='invasive_margin'))
get.counts.per.mm2(extract.ROI(dataset,ROI='stroma'))
get.counts.per.mm2.noncollapsed(extract.ROI(dataset,ROI='stroma'))


############################################################################
### Extract spatial features for downstream processing
spat_features <- extract.features(dataset)

