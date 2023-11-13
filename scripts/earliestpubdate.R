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

gene_names <- unique(gene_info$hgnc_symbol)
publications <- list()
consolidated_strings <- gene_info %>%
  group_by(hgnc_symbol) %>%
  summarize(consolidated_strings = paste(unique(c(hgnc_symbol, external_synonym)), collapse = ' OR ')) %>%
  pull(consolidated_strings)




base_directory <- '/Users/quinlan/Documents/Git/STRchive/data/'

perform_pubmed_query <- function(gene_info) {
  gene_names <- unique(gene_info$hgnc_symbol)
  file_paths <- list()  # Initialize the list to store all publications
  consolidated_strings <- gene_info %>%
    group_by(hgnc_symbol) %>%
    summarize(consolidated_strings = paste(unique(c(hgnc_symbol, external_synonym)), collapse = ' OR ')) %>%
    pull(consolidated_strings)

  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    cat("Processing gene:", gene_name, "\n")

    # Use str_detect to check if gene_name is present in each consolidated string
    idx <- str_detect(consolidated_strings, regex(gene_name, ignore_case = TRUE))

    # Filter out the relevant consolidated strings
    or_terms <- paste(consolidated_strings[idx], collapse = ' OR ')

    # Construct the query with 'AND' logic among [Title/Abstract] segments
    query <- paste0('("repeat expansion"[Title/Abstract] OR "tandem repeat"[Title/Abstract] OR "repeat expansions"[Title/Abstract] OR "tandem repeats"[Title/Abstract]) AND (', or_terms, ')[Title/Abstract] AND "English"[Language] AND ("disease"[Title/Abstract] OR "disorder"[Title/Abstract] OR "diseases"[Title/Abstract] OR "disorders"[Title/Abstract] OR "syndrome"[Title/Abstract] OR "syndromes"[Title/Abstract]) AND "journal article"[Publication Type] NOT "review"[Publication Type]')

    # Clean up any unnecessary slashes from the query
    query <- gsub("\"", "", query, fixed = TRUE)
    query <- gsub("  ", " ", query)  # Remove double spaces
    print(query)

    out_file <- paste0(base_directory, gene_name)

    # Include a separator ("/") between base_directory and gene_name
    # Modify dest_file_prefix to include the full file path
    out.A <- batch_pubmed_download(pubmed_query_string = query,
                                   format = "xml",
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


publications_df <- data.frame(
  gene_name = names(all_publications),
  publications = I(all_publications),
  stringsAsFactors = FALSE
)

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


gene_names <- names(all_publications)
publications_df <- data.frame(GeneName = character(), TypeCount = numeric())  # Initialize dataframe for storing results

for (gene_name in all_publications) {
  current_publications <- all_publications[[gene_name]]

  # Extract type counts from character strings and create a new dataframe
  type_counts <- sapply(str_extract_all(current_publications, "\\[\\d+\\]"), function(x) as.numeric(str_extract(x, "\\d+")))

  # Find the maximum length to ensure consistent lengths
  max_length <- max(lengths(type_counts))

  # Pad type_counts to match the maximum length
  type_counts <- lapply(type_counts, function(x) c(x, rep(NA, max_length - length(x))))

  # Combine the results into a matrix
  type_counts_matrix <- do.call(rbind, type_counts)

  # Create a dataframe from the matrix and add GeneName column
  type_counts_df <- data.frame(GeneName = gene_name, TypeCount = colSums(type_counts_matrix, na.rm = TRUE))

  # Append to the overall dataframe
  publications_df <- rbind(publications_df, type_counts_df)
}

# Now, publications_df contains the desired information

for (i in 1:length(pubmed_results)) {
  gene <- unique(gene_info$hgnc_symbol)[i]

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
