setwd("C:/work/projects/imaging/Iris_repo/Iris/test")
rm(list=ls())
gc()
require(Iris)
require(RColorBrewer)

raw_data <- Iris()
raw_data<- read.raw(raw_data,
                    raw_dir_name='../../../R_package/test/',
                    format='Mantra')


#apply all the thresholds PD1 for T and other cells, PD-L1 for macrophages and tumor cells
dataset <- threshold.dataset(raw_data,
                             marker='PD-Ligand-1 (Opal 690)',
                             marker_name='PDL1',
                             base=c('SOX10+','OTHER'))
dataset <- threshold.dataset(dataset,
                             marker='PD-1 (Opal 540)',
                             marker_name='PD1',
                             base=c('CD8+','OTHER'))

#get the counts
get.counts(dataset)

#run the interaction analysis
dataset <- extract.interactions(dataset)
get.interactions(dataset,'CD8+ PD1+')

#plotting interaction summaries
plot.interactions(dataset,"SOX10+ PDL1+",xlim_fix=4)

#plotting interaction maps
int_markers <- c('CD8+ PD1+','SOX10+ PDL1+')
int_marker_cols <- c('#dd1c77','#99d8c9')
silent_markers <- c('CD8+ PD1-')
silent_col=c('yellow')
interaction.maps(dataset,int_markers,int_marker_cols,silent_markers,silent_col)

#running the proximity / touching analysis
touching <- extract_touches(dataset,raw=raw)

