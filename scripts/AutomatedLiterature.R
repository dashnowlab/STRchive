#install.packages("jsonlite")
library(jsonlite)
library(dplyr)
library(biomaRt)
library(rentrez)
# solution for potential error based on library versions
#devtools::install_version("dbplyr", version = "2.3.4")
#library(dbplyr)
library(easyPubMed)
library(stringr)
library(biomaRt)
library(rentrez)
# solution for potential error based on library versions
#devtools::install_version("dbplyr", version = "2.3.4")

### Data Setup
# change to STRchive directory
data <- fromJSON("/Users/quinlan/Documents/Git/STRchive-1/STRchive/data/STRchive-database.json")


# Extract STRchive gene names into a list
gene_list <- as.character(unique(data$gene))

# Filter out NA values, if there are any, from the 'gene' column
filtered_data <- data[!is.na(data$gene), ]

# Extract STRchive gene names into a list
gene_list <- as.character(unique(filtered_data$gene))

### Use biomaRt to get gene synonyms
# This is necessary because gene names have evolved over time and can differ
# by publication
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
                host = "https://www.ensembl.org")
