\name{log_dat}
\alias{log_dat}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("log_dat")}
\format{
  The format is:
Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
  .. .. ..@ name             : chr ""
  .. .. ..@ lab              : chr ""
  .. .. ..@ contact          : chr ""
  .. .. ..@ title            : chr ""
  .. .. ..@ abstract         : chr ""
  .. .. ..@ url              : chr ""
  .. .. ..@ pubMedIds        : chr ""
  .. .. ..@ samples          : list()
  .. .. ..@ hybridizations   : list()
  .. .. ..@ normControls     : list()
  .. .. ..@ preprocessing    : list()
  .. .. ..@ other            : list()
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ assayData        :<environment: 0x000000000f7eb0a8> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	17 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:17] "BAM.File.Name" "Accession" "Acc.Tube" "USUBJID" ...
  .. .. ..@ data             :'data.frame':	36 obs. of  17 variables:
  .. .. .. ..$ BAM.File.Name: Factor w/ 157 levels "x355154.bam",..: 58 105 116 9 16 45 97 92 50 76 ...
  .. .. .. ..$ Accession    : Factor w/ 157 levels "AD505621","AE103210",..: 13 42 46 9 39 45 65 11 50 48 ...
  .. .. .. ..$ Acc.Tube     : Factor w/ 157 levels "AD505621-6","AE103210-5",..: 13 42 46 9 39 45 65 11 50 48 ...
  .. .. .. ..$ USUBJID      : Factor w/ 111 levels "CA209064-1-64001",..: 47 24 57 109 58 59 60 103 100 62 ...
  .. .. .. ..$ SITEID       : int [1:36] 3 2 4 8 4 4 4 7 6 4 ...
  .. .. .. ..$ SAFFL        : Factor w/ 1 level "Y": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ ACTARM       : Factor w/ 2 levels "IPILIMUMAB-NIVOLUMAB",..: 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ ACTARMCD     : Factor w/ 2 levels "A","B": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ Visit        : Factor w/ 6 levels "SCREEN PRE AMD 03",..: 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ Visit.Code   : Factor w/ 4 levels "OTHER","SCREEN",..: 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ RACE         : Factor w/ 3 levels "BLACK OR AFRICAN AMERICAN",..: 3 3 3 3 3 3 3 1 3 3 ...
  .. .. .. ..$ SEX          : Factor w/ 2 levels "F","M": 2 1 2 1 1 2 2 2 2 1 ...
  .. .. .. ..$ AGE          : int [1:36] 62 62 65 27 42 74 58 57 79 64 ...
  .. .. .. ..$ BRAFSTA      : Factor w/ 3 levels "MUTANT","NOT REPORTED",..: 1 3 3 1 3 3 3 1 3 3 ...
  .. .. .. ..$ MELSUBT      : Factor w/ 5 levels "ACRAL","CUTANEOUS",..: 2 3 2 2 2 2 2 2 5 2 ...
  .. .. .. ..$ SUBJID       : chr [1:36] "64008" "64021" "64034" "64036" ...
  .. .. .. ..$ response     : Factor w/ 4 levels "CR","PD","PR",..: 2 2 1 3 3 3 3 2 3 3 ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	10 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:10] NA NA NA NA ...
  .. .. ..@ data             :'data.frame':	55094 obs. of  10 variables:
  .. .. .. ..$ ensembl_gene_id     : chr [1:55094] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457" ...
  .. .. .. ..$ hgnc_symbol         : chr [1:55094] "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
  .. .. .. ..$ wikigene_name       : chr [1:55094] "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
  .. .. .. ..$ wikigene_description: chr [1:55094] "tetraspanin 6" "tenomodulin" "dolichyl-phosphate mannosyltransferase subunit 1, catalytic" "SCY1 like pseudokinase 3" ...
  .. .. .. ..$ chromosome_name     : chr [1:55094] "X" "X" "20" "1" ...
  .. .. .. ..$ strand              : int [1:55094] -1 1 -1 -1 1 -1 1 -1 -1 1 ...
  .. .. .. ..$ band                : chr [1:55094] "q22.1" "q22.1" "q13.13" "q24.2" ...
  .. .. .. ..$ start_position      : int [1:55094] 100627109 100584802 50934867 169849631 169662007 27612064 196651878 143494811 53497341 41072945 ...
  .. .. .. ..$ end_position        : int [1:55094] 100639991 100599885 50958555 169894267 169854080 27635277 196747504 143511690 53616970 41099976 ...
  .. .. .. ..$ transcript_length   : int [1:55094] 2206 1339 1075 6364 3849 2729 4127 2356 3774 3811 ...
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ annotation       : chr "RNASeq raw counts"
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	36 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. ..@ .Data:List of 4
  .. .. .. ..$ : int [1:3] 3 2 2
  .. .. .. ..$ : int [1:3] 2 30 0
  .. .. .. ..$ : int [1:3] 1 3 0
  .. .. .. ..$ : int [1:3] 1 0 0
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(log_dat)
## maybe str(log_dat) ; plot(log_dat) ...
}
\keyword{datasets}
