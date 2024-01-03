library(biomaRt)
#install.packages("rentrez")
library(rentrez)
library(dplyr)
#install.packages("easyPubMed")
library(easyPubMed)
library(stringr)

setwd('/Users/quinlan/Documents/Git/STRchive/data')
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

excluded_synonym_list <- c("B37", "MHP", "MED", "DM", "DM1", "FA", "GAC", "SPD",
                           "PRP", "A1", "CCD", "PHP", "VCF")

gene_info <- gene_info %>%
  filter(!grepl(paste(excluded_synonym_list, collapse = '|'), external_synonym))


gene_info <- gene_info %>%
  mutate(across(c(external_synonym, hgnc_symbol), ~sprintf('"%s"', .)))

gene_names <- unique(gene_info$hgnc_symbol)
publications <- list()
consolidated_strings <- gene_info %>%
  group_by(hgnc_symbol) %>%
  summarize(consolidated_strings = paste(unique(c(hgnc_symbol, external_synonym)), collapse = ' OR ')) %>%
  pull(consolidated_strings)
consolidated_strings <- gsub("BMD", "Becker muscular dystrophy", consolidated_strings)
base_directory <- '/Users/quinlan/Documents/Git/STRchive/data/'

perform_pubmed_query <- function(gene_info) {
  file_paths <- list()  # Initialize the list to store all publications
  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    cat("Processing gene:", gene_name, "\n")
    # Use str_detect to check if gene_name is present in each consolidated string
    idx <- str_detect(consolidated_strings, regex(gene_name, ignore_case = TRUE))

    # Filter out the relevant consolidated strings
    or_terms <- paste(consolidated_strings[idx], collapse = ' OR ')

    # Construct the query with 'AND' logic among [Title/Abstract] segments
    query <- paste0('("repeat expansion"[Title/Abstract] OR "tandem repeat"[Title/Abstract] OR "repeat expansions"[Title/Abstract] OR "tandem repeats"[Title/Abstract]) OR "repeat sequence"[Title/Abstract] OR "repeat sequences"[Title/Abstract] OR "repeat length"[Title/Abstract] OR "repeat lengths"[Title/Abstract] AND (', or_terms, ')[Title/Abstract] AND "English"[Language] AND ("disease"[Title/Abstract] OR "disorder"[Title/Abstract] OR "diseases"[Title/Abstract] OR "disorders"[Title/Abstract] OR "syndrome"[Title/Abstract] OR "syndromes"[Title/Abstract] OR "patient"[Title/Abstract] OR "patients"[Title/Abstract] OR "proband"[Title/Abstract] OR "probands"[Title/Abstract]) AND "journal article"[Publication Type] NOT "review"[Publication Type]')

    # Clean up any unnecessary slashes from the query
    query <- gsub("  ", " ", query)  # Remove double spaces
    print(query)
    gene_name <- gsub('"', '', gene_name)
    out_file <- paste0(base_directory, gene_name)

    # Include a separator ("/") between base_directory and gene_name
    # Modify dest_file_prefix to include the full file path
    out.A <- batch_pubmed_download(pubmed_query_string = query,
                                   format = "medline",
                                   batch_size = 10000,
                                   dest_file_prefix = out_file,
                                   encoding = "ASCII")
    # the function adds 01.txt so, gotta fix that here
    out_file <- paste0(base_directory, gene_name, "01.txt")
    print(out_file)

    # Check if the file was created successfully
    cat("Full file path:", out_file, "\n")
    if (file.exists(out_file)) {
      cat("File exists.\n")
      file_paths[[gene_name]] <- out_file
    } else {
      cat("Error: File not found -", out_file, "\n")
    }
  }

  return(file_paths)
}

file_paths <- perform_pubmed_query(gene_info)

all_publications <- list()

# Assuming file_paths is a list of file paths
for (gene_name in names(file_paths)) {
  # Append "01.txt" to the file path
  file_path <- paste0(file_paths[[gene_name]])

  tryCatch({
    # Read the file into a character vector
    current_publications <- readLines(file_path)

    # Append to the overall list
    all_publications[[gene_name]] <- current_publications
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n")
    # Append an empty character vector in case of an error
    all_publications[[gene_name]] <- character(0)
  })
}

extract_pub_info <- function(medline_data_list, gene_name) {
  # Combine the list of XML strings into a single string
  medline_string <- paste(medline_data_list, collapse = "")

  # Use regular expressions to extract PMID and publication years
  # Extract PMIDs
  pmids <- str_extract_all(medline_string, "(?<=PMID- )\\d+")[[1]]

  # Extract Publication Years
  publication_years <- str_extract_all(medline_string, "(?<=EDAT- )\\d{4}")[[1]]

  # Ensure both vectors have the same length
  length_diff <- length(pmids) - length(publication_years)
  if (length_diff > 0) {
    publication_years <- c(publication_years, rep(NA, length_diff))
  } else if (length_diff < 0) {
    pmids <- c(pmids, rep(NA, -length_diff))
  }

  # Create a dataframe with gene_name, PMID, and PublicationYear
  pub_info_df <- data.frame(GeneName = rep(gene_name, length(pmids)),
                            PMID = pmids,
                            PublicationYear = publication_years,
                            stringsAsFactors = FALSE)

  return(pub_info_df)
}


# Initialize an empty list to store the results
pub_info_list <- list()

# Loop through each gene_name in all_publications
for (gene_name in names(all_publications)) {
  # Get the list of XML data for the current gene_name
  medline_data_list <- all_publications[[gene_name]]
  print(gene_name)
  # Extract publication information using the function
  pub_info_df <- extract_pub_info(medline_data_list, gene_name)

  # Append the results to the list
  pub_info_list[[gene_name]] <- pub_info_df
}

# Combine all the dataframes into a single dataframe
all_pub_info_df <- do.call(rbind, pub_info_list)

# Get the current date
current_date <- format(Sys.Date(), "%Y%m%d")

# Concatenate the date to the file name
file_name <- paste0("/Users/quinlan/Documents/Git/STRchive/data/all_pub_info_", current_date, ".tsv")

# Write the table with the updated file name
write.table(all_pub_info_df, file_name, sep = "\t", quote = FALSE, row.names = FALSE)
