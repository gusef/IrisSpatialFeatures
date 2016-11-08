setwd("C:/work/projects/imaging/Iris_repo/Iris/test")
rm(list=ls())
gc()


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
                             base=c('OTHER'))





