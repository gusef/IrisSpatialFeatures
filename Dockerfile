FROM zamora/r-devtools

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends \
		libtiff5-dev \
		fftw-dev \
	&& rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'source("https://bioconductor.org/biocLite.R");biocLite()'
RUN Rscript -e 'library(devtools);install_github("gusef/Iris")'
RUN echo "library(IrisSpatialFeatures)" > $HOME/.Rprofile
