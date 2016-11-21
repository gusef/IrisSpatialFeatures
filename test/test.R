setwd("C:/work/projects/imaging/Iris/test")
rm(list=ls())
gc()
require(Iris)


raw_data <- Iris()
raw_data<- read.raw(raw_data,
                    raw_dir_name='data/',
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
get.counts.per.mm2(dataset)

get.counts.collapsed(dataset)
extract.counts.noncollapsed(dataset)



#run the nearest neighbor analysis
dataset <- extract.nearest.neighbor(dataset)
get.nearest.neighbors(dataset,"SOX10+ PDL1+")
plot.nearest.neighbor(dataset,'CD8+ PD1+','SOX10+ PDL1')
plot.nearest.neighbor(dataset,'SOX10+ PDL1+','SOX10+ PDL1-')

#run the interaction analysis
dataset <- extract.interactions(dataset)
get.interactions(dataset,'CD8+ PD1+')

#plotting interaction summaries
plot.interactions(dataset,'CD8+ PD1+',xlim_fix=4)

#plotting interaction maps
int_markers <- c('CD8+ PD1+','SOX10+ PDL1+')
int_marker_cols <- c('#dd1c77','#99d8c9')
silent_markers <- c('CD8+ PD1-')
silent_col=c('yellow')
interaction.maps(dataset,int_markers,int_marker_cols,silent_markers,silent_col)

#running the proximity / touching analysis
dataset <- extract.proximity(dataset)
plot.proximities(dataset,"SOX10+ PDL1-",xlim_fix=4)

#get interactions again
tumor_area <- extract.ROI(dataset,ROI='tumor')
get.counts.per.mm2(tumor_area)
tumor_area <- extract.interactions(tumor_area)
plot.interactions(tumor_area,"SOX10+ PDL1-",xlim_fix=4)



############################################################################
### Extract annotation

#correlation with outcome
annotation <- read.xlsx('Case list with response.xlsx',1)
annotation$sample_name <- paste0('MEL',annotation$Case)
annotation <- annotation[match(colnames(dat),annotation$sample_name),]


res <- feature_selection(dat,annotation$Response)

write.xlsx(res$t_test,file='automatic_feature_selection/Automatic_feature_selection_Mike_bond.xlsx',sheetName = 'Limma')
write.xlsx(res$wilcox,file='automatic_feature_selection/Automatic_feature_selection_Mike_bond.xlsx',sheetName = 'Wilcox',append = T)

topbar_cols <- col <- c('#af8dc3','#7fbf7b')[as.numeric(annotation$Response)]
markers <- as.character(res$wilcox[1:20,1])




