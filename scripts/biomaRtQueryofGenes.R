library(biomaRt)
#install.packages("rentrez")
library(rentrez)
library(dplyr)

# Define the list of genes
gene_list <- c("AFF2", "FMR2", "AFF3", "AR", "ARX", "ARX", "ATN1", "ATXN1", "ATXN10", "ATXN2",
               "ATXN3", "ATXN7", "ATXN8OS", "ATXN8", "BEAN1", "C9orf72", "CACNA1A", "CBL2", "COMP",
               "DAB1", "DIP2B", "DMD", "DMPK", "FMR1", "FOXL2", "FXN", "GIPC1", "GLS", "HOXA13",
               "HOXA13", "HOXA13", "HOXD13", "HTT", "JPH3", "LRP12", "MARCH6", "NIPA1", "NOP56",
               "NOTCH2NLC", "NUTM2B-AS1", "PABPN1", "PHOX2B", "POLG", "PPP2R2B", "PRDM12", "RAPGEF2",
               "RFC1", "RILPL1", "RUNX2", "SAMD12", "SOX3", "STARD7", "TBP", "TBX1", "TCF4", "TNRC6A",
               "XYLT1", "YEATS2", "ZIC2", "ZIC3", "ZNF713", "ZNF9", "CNBP1", "CSTB", "EIF4A3", "PRNP",
               "VWA1", "ABCD3", "FGF14")

# Specify the dataset to use (in this case, HGNC in human, using GRCh38 assembly)
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

# Retrieve the HGNC symbol and Gene Synonym for the genes in the list
gene_info <- getBM(attributes = c("hgnc_symbol", "external_synonym"),
                   filters = "hgnc_symbol",
                   values = gene_list,
                   mart = mart)

#if I want to copy this whole thing
#clipr::write_clip(gene_info)
# Replace missing values in 'external_synonym' with 'hgnc_symbol'
gene_info$external_synonym[is.na(gene_info$external_synonym) | gene_info$external_synonym == ""] <- gene_info$hgnc_symbol[is.na(gene_info$external_synonym) | gene_info$external_synonym == ""]

## making the query to get full article data
perform_pubmed_query <- function(gene_info) {
  #set up publications to be filled out
  publications <- list()

  #I have to do the queries separately or it breaks things
  for (i in 1:nrow(gene_info)) {
    gene_name <- gene_info$hgnc_symbol[i]
    gene_synonyms <- unlist(strsplit(gene_info$external_synonym[i], ", "))

    #basically, I want to group all gene synonyms with the root gene name
    search_terms <- c(gene_name, gene_synonyms)
    search_terms <- unique(unlist(search_terms))


    query <- paste0('("repeat expansion" OR "tandem repeat") AND (', paste(search_terms, collapse = ' OR '), ') AND "english"[Language] AND ("disease" OR "disorder")')
    # Remove unnecessary slashes from search terms
    query <- gsub("\"", "", query, fixed = TRUE)
    search_results <- entrez_search(db = "pubmed", term = query)
    if (search_results$count > 0) {
      article_ids <- entrez_search(db = "pubmed", term = query)$ids
      articles <- entrez_fetch(db = "pubmed", id = article_ids, rettype = "xml")
      publications[[gene_name]] <- articles
    } else {
      publications[[gene_name]] <- NULL
    }
  }

  return(publications)
}

# Assuming 'gene_info_subset' contains your gene names and synonyms
pubmed_results <- perform_pubmed_query(gene_info_subset)

# Printing the PubMed query results for each gene
for (gene_name in names(pubmed_results)) {
  print(paste("Gene:", gene_name))
  print(pubmed_results[[gene_name]])
}


gene_info_subset <- gene_info[1:6, ]

get_pmids <- function(gene_info) {
  pmid_list <- list()

  for (i in 1:nrow(gene_info)) {
    gene_name <- gene_info$hgnc_symbol[i]
    gene_synonyms <- unlist(strsplit(gene_info$external_synonym[i], ", "))

    search_terms <- c(gene_name, gene_synonyms)
    search_terms <- unique(unlist(search_terms))

    query <- paste0('("repeat expansion" OR "tandem repeat") AND (',
                    paste(search_terms, collapse = ' OR '), ') AND "english"[Language
                    ] AND ("disease"[Title/Abstract] OR "disorder"[Title/Abstract] OR "
                    diseases"[Title/Abstract] OR "disorders"[Title/Abstract] OR "
                    syndrome"[Title/Abstract] OR "syndromes"[Title/Abstract]) AND "
                    journal article"[Publication Type])')
    query <- gsub("\"", "", query, fixed = TRUE)

    search_results <- entrez_search(db = "pubmed", term = query)

    if (search_results$count > 0) {
      article_ids <- entrez_search(db = "pubmed", term = query, retmax = 10000)$ids
      pmid_list[[gene_name]] <- article_ids
    } else {
      pmid_list[[gene_name]] <- NULL
    }
  }

  return(pmid_list)
}

# Assuming 'gene_info_subset' contains your gene names and synonyms
pmid_results <- get_pmids(gene_info_subset)

# Printing the retrieved PMIDs for each gene
for (gene_name in names(pmid_results)) {
  print(paste("Gene:", gene_name))
  print(pmid_results[[gene_name]])
}
