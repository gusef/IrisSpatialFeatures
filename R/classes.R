raw_data <- setClass("raw_data",
                     slots = c(data='data.frame',
                               summary='matrix',
                               score='matrix',
                               mem_seg_map='matrix',
                               nuc_seg_map='matrix',
                               dapi_map='matrix'))

setOldClass("ppp")
Coordinate <- setClass("Coordinate",
                       slots = c(ppp = "ppp",
                                 raw = "raw_data",
                                 mask = "list",
                                 coordinate_name="character"),
                       prototype=list(
                           ppp = ppp()))

Sample <- setClass("Sample",
                   slots = c(coordinates = "list",
                             sample_name="character"))

#' @export
Iris <- setClass("Iris",
                 slots = c(samples = "list",
                           pData = "data.frame",
                           counts = "list",
                           nearest_neighbors = "list",
                           interactions = "list",
                           proximity = "list",
                           microns_per_pixel="numeric",
                           markers="character"))

