# Iris


## How to install in R:

This is a private repository so you need to generate a personal access token (PAT) here: https://github.com/settings/tokens which lets you install R packages from private Github repositories.

In R they can then just install it using devtools:


install.packages('devtools')
 
#load devtools and install the package
library(devtools)
install_github("gusef/Iris",auth_token = 'COPY_PASTE_YOUR_TOKEN_HERE')



## How to use the package:
There is a vignette included in the 'vignette' directory of the package. In addition we have 2 example datasets on Dropbox that show the full capability of the package. 