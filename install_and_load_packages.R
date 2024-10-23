# List of required CRAN packages
cran_packages <- c(
  "plyr", "dplyr", "readr", "ggplot2", "reshape", "pheatmap",
  "gridExtra", "grid", "cowplot", "ggrepel", "hexbin", "tidyr", "tibble"
)

# List of required Bioconductor packages
bioc_packages <- c("vsn", "CARNIVAL", "progeny", "dorothea", "limma", "viper", "cosmosR")

# Install BiocManager if not already installed
if (!require("BiocManager")) {
  install.packages("BiocManager")
  library(BiocManager)
}

# Function to check if a CRAN package is installed, and install it if it is not
install_if_missing_cran <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# Function to check if a Bioconductor package is installed, and install it if it is not
install_if_missing_bioc <- function(package) {
  if (!require(package, character.only = TRUE)) {
    BiocManager::install(package, ask = FALSE, update = FALSE)
    library(package, character.only = TRUE)
  }
}

# Install/load Bioconductor packages
for (pkg in bioc_packages) {
  install_if_missing_bioc(pkg)
}

# Install OmnipathR from the devel version if not already installed
if (!require("OmnipathR", character.only = TRUE)) {
  BiocManager::install("OmnipathR", force=TRUE)
  library(OmnipathR, character.only = TRUE)
}

# Install/load CRAN packages
for (pkg in cran_packages) {
  install_if_missing_cran(pkg)
}

# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }

# Install OmnipathR version 3.9.6 from GitHub
# library(remotes)
# remotes::install_github("saezlab/OmnipathR")

