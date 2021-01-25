# load required packages
required_packages <- c("rgdal","readr",
                       "tidyverse",
                       "spdep",
                       "sf",
                       "surveillance",
                       "hhh4contacts",
                       "socialmixr",
                       "ggplot2",
                       "viridis",
                       "ggrepel")

for (i in required_packages) {
  if (!i %in% rownames(installed.packages())){
    install.packages(i)
  }
}

lapply(required_packages, library, character.only = TRUE)
