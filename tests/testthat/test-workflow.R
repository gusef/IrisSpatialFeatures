context("workflow")
test_that("steps of a workflow excute", {
    # Test read_raw
    raw_data <- new("ImageSet")
    raw_data<- read_raw(raw_data,
                        raw_dir_name=system.file("extdata", package = "IrisSpatialFeatures"),
                        format='Mantra')
    expect_that(length(raw_data@markers),equals(3))
    # Test threshold_dataset
    dataset <- threshold_dataset(raw_data,
                                 marker='PD-Ligand-1 (Opal 690)',
                                 marker_name='PDL1',
                                 base=c('SOX10+'))
    dataset <- threshold_dataset(dataset,
                                 marker='PD-1 (Opal 540)',
                                 marker_name='PD1',
                                 base=c('CD8+','OTHER'))

    expect_that(length(dataset@markers),equals(6))
    #
})
