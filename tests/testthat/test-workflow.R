context("workflow")
test_that("steps of a workflow excute", {
    # Test read_raw
    raw_data <- new("ImageSet")
    raw_data<- read_raw(raw_data,
                        raw_dir_name=system.file("extdata", package = "IrisSpatialFeatures"),
                        format='Mantra')
    df = data.frame(observed=raw_data@markers,expected=c('CD8+','OTHER','SOX10+'))
    diff = length(df[df$expected!=df$observed,1])
    expect_that(0,equals(diff))
    # Test threshold_dataset
    dataset <- threshold_dataset(raw_data,
                                 marker='PD-Ligand-1 (Opal 690)',
                                 marker_name='PDL1',
                                 base=c('SOX10+'))
    dataset <- threshold_dataset(dataset,
                                 marker='PD-1 (Opal 540)',
                                 marker_name='PD1',
                                 base=c('CD8+','OTHER'))
    df = data.frame(observed=dataset@markers,expected=c("CD8+ PD1-",
                                                        "CD8+ PD1+",
                                                        "OTHER PD1-",
                                                        "OTHER PD1+",
                                                        "SOX10+ PDL1-",
                                                        "SOX10+ PDL1+"))
    diff = length(df[df$expected!=df$observed,1])
    expect_that(0,equals(diff))
    #
})
