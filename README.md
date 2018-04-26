# Iris Spatial Features


## How to install in R:

``` r
#install devtools from cran
install.packages('devtools')
 
#load devtools and install the package
library(devtools)
install_github("gusef/IrisSpatialFeatures")

```

## How to use the package:
There is a vignette included in the 'vignette' directory of the package. 

### Common task walkthroughs

#### 1. Read in all samples

Uses the default mask to define the ROI if one is included

```r
data <- read_raw('Test ROI')
```

#### 2. Read in all samples and a custom ROI mask

Assumes every frame of every sample has *_ROI.tif with a defined area.

```r
data <- read_raw('Test ROI',customMask="ROI")
```

#### 3. Extract counts/mm2 from the data for samples and individual frames

```r
sample_count_density <- counts_per_mm2_sample_data_frame(data)
frame_count_density <- counts_per_mm2_data_frame(data)
```

#### 4. Extract raw counts from the data for samples and individual frames

```r
sample_count <- counts_data_sample_frame(data)
frame_count <- counts_data_frame(data)
```

#### 5. Analyze counts in samples with a tumor and margin defined

Assumes each sample and frame has *_Tumor.tif and *_Invasive_Margin.tif with defined areas

```r
data <- read_raw('Test tumor IM mixed case copy',
                 readTumorAndMarginMasks=TRUE)
tumor <- extract_ROI(data,'tumor')
stroma <- extract_ROI(data,'stroma')
invasive_margin <- extract_ROI(data,'invasive_margin')

tumor_counts <- counts_per_mm2_sample_data_frame(tumor)
stroma_counts <- counts_per_mm2_sample_data_frame(stroma)
invasive_margin_counts <- counts_per_mm2_sample_data_frame(invasive_margin)
```

#### 6. Analyze the complete tumor in samples with a tumor and margin defined

Only requires *_Tumor.tif for each frame. We use it as a custom mask and don't concern ourselves with the Invasive Margine files.

```r
tumor <- read_raw('Test tumor IM mixed case copy', 
                 customMask='Tumor')
tumor_counts <- counts_per_mm2_sample_data_frame(tumor)
```
