library(biomaRt)
#install.packages("rentrez")
library(rentrez)
library(dplyr)
#install.packages("easyPubMed")
library(easyPubMed)
library(stringr)


# Define the list of genes
gene_list <- c("AFF2", "FMR2", "AFF3", "AR", "ARX", "ARX", "ATN1", "ATXN1", "ATXN10", "ATXN2",
               "ATXN3", "ATXN7", "ATXN8OS", "ATXN8", "BEAN1", "C9orf72", "CACNA1A", "CBL2", "COMP",
               "DAB1", "DIP2B", "DMD", "DMPK", "FMR1", "FOXL2", "FXN", "GIPC1", "GLS", "HOXA13",
               "HOXA13", "HOXA13", "HOXD13", "HTT", "JPH3", "LRP12", "MARCH6", "NIPA1", "NOP56",
               "NOTCH2NLC", "NUTM2B-AS1", "PABPN1", "PHOX2B", "POLG", "PPP2R2B", "PRDM12", "RAPGEF2",
               "RFC1", "RILPL1", "RUNX2", "SAMD12", "SOX3", "STARD7", "TBP", "TBX1", "TCF4", "TNRC6A",
               "XYLT1", "YEATS2", "ZIC2", "ZIC3", "ZNF713", "ZNF9", "CNBP1", "CSTB", "EIF4A3", "PRNP",
               "VWA1", "ABCD3", "FGF14")

# Set a longer timeout for the Ensembl service
options(timeout = 240)

# Attempt to useMart with a longer timeout
try({
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
}, silent = TRUE)

# Reset the timeout option to its default value
options(timeout = NULL)

# Check if the mart object is created
if (exists("mart")) {
  # Retrieve the HGNC symbol and Gene Synonym for the genes in the list
  gene_info <- getBM(attributes = c("hgnc_symbol", "external_synonym"),
                     filters = "hgnc_symbol",
                     values = gene_list,
                     mart = mart)

} else {
  cat("Failed to create the mart object. Check your internet connection and try again.")
}

#if I want to copy this whole thing
#clipr::write_clip(gene_info)
# Replace missing values in 'external_synonym' with 'hgnc_symbol'
gene_info$external_synonym[is.na(gene_info$external_synonym) | gene_info$external_synonym == ""] <- gene_info$hgnc_symbol[is.na(gene_info$external_synonym) | gene_info$external_synonym == ""]

perform_pubmed_query <- function(gene_info) {
  publications <- list()

  for (i in 1:nrow(gene_info)) {
    gene_name <- gene_info$hgnc_symbol[i]
    gene_synonyms <- unlist(strsplit(gene_info$external_synonym[i], ", "))

    # Combine hgnc_symbol and external_synonym for the OR-based query
    all_terms <- c(gene_name, gene_synonyms)
    all_terms <- unique(unlist(all_terms))

    # Combine all the terms into a single OR statement without adding [Title/Abstract]
    or_terms <- paste(all_terms, collapse = ' OR ')

    # Construct the query with 'AND' logic among [Title/Abstract] segments
    query <- paste0('("repeat expansion"[Title/Abstract] OR "tandem repeat"[Title/Abstract] OR "repeat expansions"[Title/Abstract] OR "tandem repeats"[Title/Abstract]) AND (', or_terms, ')[Title/Abstract] AND "English"[Language] AND ("disease"[Title/Abstract] OR "disorder"[Title/Abstract] OR "diseases"[Title/Abstract] OR "disorders"[Title/Abstract] OR "syndrome"[Title/Abstract] OR "syndromes"[Title/Abstract]) AND "journal article"[Publication Type] NOT "review"[Publication Type]')

    # Clean up any unnecessary slashes from the query
    query <- gsub("\"", "", query, fixed = TRUE)
    print(query)

    search_results <- entrez_search(db = "pubmed", term = query, retmax = 500)
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


pubmed_results <- perform_pubmed_query(gene_info)



earliest_pub_dates_df <- data.frame(
  Gene = character(0),
  Earliest_Pub_Date = numeric(0),
  NumPub = numeric(0),
  Titles = I(list())
)

for (i in 1:length(pubmed_results)) {
  gene <- names(pubmed_results)[i]

  # Extract publication years for each gene using custom_grep
  publication_years <- custom_grep(pubmed_results[[i]], "PubDate", "char")
  publication_years <- str_extract(publication_years, "(?<=<Year>)\\d{4}(?=</Year>)")

  # Convert years to numeric format
  publication_years <- as.numeric(publication_years)

  # Find the minimum publication year
  min_year <- min(publication_years, na.rm = TRUE)

  # Extract titles
  titles <- custom_grep(pubmed_results[[i]], "ArticleTitle", "char")
  print(titles)
  # Count the number of titles
  titles_count <- length(titles)

  # Combine titles into a single string
  titles_string <- paste(titles, collapse = ", ")

  # Add values to the data frame
  row_index <- nrow(earliest_pub_dates_df) + 1
  earliest_pub_dates_df[row_index, ] <- list(gene, min_year, titles_count, I(titles_string))
}


write.table(earliest_pub_dates_df, "/Users/quinlan/Documents/Git/STRchive/scripts/earliest_pub_dates.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(earliest_pub_dates_df[, 1:3], "/Users/quinlan/Documents/Git/STRchive/scripts/earliest_pub_dates_only.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
