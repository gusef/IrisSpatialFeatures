

Iris <- setClass("Iris",
                 slots = c(marray = "matrix",
                           fmeta = "data.frame",
                           pmeta = "data.frame"))


# raw_data$MEL12164$`081016_10`$
#     ppp
#     raw - data
#         - summary
#         - mem_seg_map
#         - nuc_seg_map
#         - dapi map
#         - score
#     masks
#        - invasive margin
#        - tumor
#        - filled margin