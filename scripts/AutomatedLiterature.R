library(jsonlite)
library(dplyr)
library(biomaRt)
library(rentrez)
# solution for potential error based on library versions
#devtools::install_version("dbplyr", version = "2.3.4")
#library(dbplyr)
library(easyPubMed)
library(stringr)
library(purrr)

### Data Setup
args <- commandArgs(trailingOnly = TRUE)

data <- fromJSON("/Users/quinlan/Documents/Git/STRchive-1/STRchive/data/STRchive-database.json")

if (length(args) < 3) {
  stop("Need input json, base directory, and output json")
}

data <- fromJSON(args[1])

# get unique references used in the json for each entry
data <- data %>%
  mutate(references = apply(data[, c("age_onset", "mechanism_detail", "details", "source",
                                     "prevalence_details", "year", "disease_description")], 1, function(row) {
                                       # Extract content inside square brackets
                                       matches <- unlist(regmatches(row, gregexpr("(?<=\\[)[^\\]]+(?=\\])", row, perl = TRUE)))

                                       # Replace "; " with "," for multiple citations within one
                                       matches <- gsub("; ", ",", matches)

                                       # Collapse matches into a single string, remove duplicates, and trim whitespaces
                                       unique_matches <- unique(trimws(matches))

                                       # Return the cleaned references
                                       paste(unique_matches, collapse = ",")
                                     }))


### Established loci lit  retrieval
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

# Retrieve the HGNC symbol and Gene Synonym for the genes in the list
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

# Replace missing values in 'external_synonym' with 'hgnc_symbol'
# This is because some gene_names have no synonyms, so it's easier to group this way
gene_info$external_synonym[is.na(gene_info$external_synonym) |
                             gene_info$external_synonym == ""] <- gene_info$hgnc_symbol[is.na(
                               gene_info$external_synonym) | gene_info$external_synonym == ""]

# curated list of synonyms to exclude based on non-specific results
excluded_synonym_list <- c("B37", "MHP", "MED", "DM", "DM1", "FA", "GAC", "SPD",
                           "PRP", "A1", "CCD", "PHP", "VCF", "PEM", "MCD", "EMA")

#remove from gene_info
gene_info <- gene_info %>%
  filter(!grepl(paste(excluded_synonym_list, collapse = '|'), external_synonym))

# adding quotation marks to terms for query
gene_info <- gene_info %>%
  mutate(across(c(external_synonym, hgnc_symbol), ~sprintf('"%s"', .)))

# adding gene synonym missed by biomaRt
gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"FMR1\"", external_synonym = "\"FMR-1\""))

#because pubmed hates slashes and these use slashes
gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"NUTM2B-AS1\"", external_synonym = "\"LOC642361/NUTM2B-AS1\""))
gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"ATXN10\"", external_synonym = "\"SCA10\""))

# collapse to gene names
gene_names <- unique(gene_info$hgnc_symbol)

# strings for queries
consolidated_strings <- gene_info %>%
  group_by(hgnc_symbol) %>%
  summarize(consolidated_strings = paste(unique(c(hgnc_symbol, external_synonym)), collapse = ' OR ')) %>%
  pull(consolidated_strings)

#substitute to avoid unrelated PMIDs
consolidated_strings <- gsub("BMD", "Becker muscular dystrophy", consolidated_strings)

#where results will be stored
base_directory <- args[2]

# function to perform the pubmed query
# Function printout includes gene name and if there are results, confirms
# that a file has been created
perform_pubmed_query <- function(gene_info) {
  file_paths <- list()  # Initialize the list to store all publications
  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    cat("Processing gene:", gene_name, "\n")
    # Use str_detect to check if gene_name is present in each consolidated string
    idx <- str_detect(consolidated_strings, regex(gene_name, ignore_case = TRUE))

    # Put together the relevant consolidated strings
    or_terms <- paste(consolidated_strings[idx], collapse = ' OR ')
    # Splitting the or_terms into individual terms
    individual_terms <- unlist(strsplit(or_terms, " OR "))

    # Adding [Title/Abstract] to each term
    joined_terms <- paste0(individual_terms, "[Title/Abstract]")
    #joined_terms <- paste0('(', paste(or_terms, collapse = '[Title/Abstract] OR '), ')[Title/Abstract]')
    # Construct the query with organized or_terms
    query <- paste0('("repeat expansion"[Title/Abstract] OR "tandem repeat"[Title/Abstract] OR "repeat expansions"[Title/Abstract] OR "tandem repeats"[Title/Abstract] OR "repeat sequence"[Title/Abstract] OR "repeat sequences"[Title/Abstract] OR "repeat length"[Title/Abstract] OR "repeat lengths"[Title/Abstract] OR "expansion"[Title] OR "expansions"[Title] OR "repeats"[Title]) AND (', paste(joined_terms, collapse = " OR "),') AND "English"[Language] AND ("disease"[Title/Abstract] OR "disorder"[Title/Abstract] OR "diseases"[Title/Abstract] OR "disorders"[Title/Abstract] OR "syndrome"[Title/Abstract] OR "syndromes"[Title/Abstract] OR "patient"[Title/Abstract] OR "patients"[Title/Abstract] OR "proband"[Title/Abstract] OR "probands"[Title/Abstract]) AND ("journal article"[Publication Type] OR "letter"[Publication Type] or "Case Reports"[Publication Type]) NOT "review"[Publication Type]')


    # Clean up any unnecessary slashes from the query
    query <- gsub("  ", " ", query)  # Remove double spaces
    print(query)
    gene_name <- gsub('"', '', gene_name)
    out_file <- paste0(base_directory, "/", gene_name)

    # Include a separator ("/") between base_directory and gene_name
    # Modify dest_file_prefix to include the full file path
    out.A <- batch_pubmed_download(pubmed_query_string = query,
                                   format = "medline",
                                   batch_size = 10000,
                                   dest_file_prefix = out_file,
                                   encoding = "ASCII")
    # the function adds 01.txt so, gotta fix that here
    out_file <- paste0(base_directory, "/", gene_name, "01.txt")
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

#creating list for publications
all_publications <- list()

# Assuming file_paths is a list of file paths
#let's get the files into a dataframe in R
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

# empty list for the publication info
pub_info_list <- list()

extract_citation_info <- function(medline_data_list, gene_name) {
  # Combine the list of XML strings into a single string
  medline_string <- paste(medline_data_list, collapse = "")

  # Use regular expressions to extract PMID, publication years, and titles
  # Extract PMIDs
  pmids <- str_extract_all(medline_string, "(?<=PMID- )\\d+")[[1]]
  #print(pmids)
  # Extract Publication Dates
  publication_dates <- str_extract_all(medline_string, "(?<=DP  - )\\d+")[[1]]
  #print(publication_dates)
  # Extract Titles
  title <- str_extract_all(medline_string, "(?<=TI  - ).+?(?=\\.|\\?)")[[1]]
  #print(title)

  # Ensure all vectors have the same length
  length_diff <- length(pmids) - length(publication_dates)
  if (length_diff > 0) {
    publication_dates <- c(publication_dates, rep(NA, length_diff))
  } else if (length_diff < 0) {
    pmids <- c(pmids, rep(NA, -length_diff))
  }

  length_diff <- length(pmids) - length(title)
  if (length_diff > 0) {
    title <- c(title, rep(NA, length_diff))
  } else if (length_diff < 0) {
    pmids <- c(pmids, rep(NA, -length_diff))
  }

  # Create a dataframe with gene_name, PMID, PublicationYear, and Title
  pub_info_df <- data.frame(gene = rep(gene_name, length(pmids)),
                            PMID = pmids,
                            PublicationDate = publication_dates,
                            Title = title,
                            stringsAsFactors = FALSE)

  return(pub_info_df)
}

for (gene_name in names(all_publications)) {
  # Get the list of XML data for the current gene_name
  medline_data_list <- all_publications[[gene_name]]
  #print(gene_name)
  # Extract publication information using the function
  pub_info_df <- extract_citation_info(medline_data_list, gene_name)

  # Append the results to the list
  pub_info_list[[gene_name]] <- pub_info_df
}

# Combine all the dataframes into a single dataframe
all_pub_info_df <- do.call(rbind, pub_info_list)

# add pubmed search results to literature field for entry
data <- data %>%
  mutate(additional_literature = map_chr(gene, function(g) {
    # Find matching PMIDs from all_pub_info_df for each gene
    matching_pmids <- all_pub_info_df %>%
      filter(gene == g) %>%
      pull(PMID) %>%
      unique() %>%
      paste0("@pmid:", ., collapse = ",")

    # If there are no matches, return an empty string
    if (length(matching_pmids) == 0) {
      return("")
    }
    return(matching_pmids)
  }))

# remove any redundant pmids from additional_literature that are in references
data <- data %>%
  mutate(additional_literature = mapply(function(lit, refs) {
    # Split the strings by commas to create lists
    lit_list <- unlist(str_split(lit, ",\\s*"))  # Split and remove extra spaces
    refs_list <- unlist(str_split(refs, ",\\s*"))  # Split and remove extra spaces

    # Find common elements between additional_literature and references
    common_elements <- intersect(lit_list, refs_list)

    # Print common elements being removed (if any)
    if (length(common_elements) > 0) {
      cat("Removing the following elements from additional_literature: ", paste(common_elements, collapse = ", "), "\n")
    }

    # Remove common elements from additional_literature
    lit_list <- setdiff(lit_list, common_elements)

    # Join the remaining elements back into a string
    cleaned_lit <- paste(lit_list, collapse = ",")

    # Return the cleaned additional_literature
    return(cleaned_lit)
  }, data$additional_literature, data$references))


#### New locus lit retrieval
#new locus query found from reviewing pertinent terms in discovery papers
perform_new_pubmed_query <- function() {
  file_path <- list()  # Initialize the list to store all publications
  #joined_terms <- paste0('(', paste(or_terms, collapse = '[Title/Abstract] OR '), ')[Title/Abstract]')
  # Construct the query with organized or_terms
  query <- paste0('("repeat expansion"[Title/Abstract] OR "tandem repeat"[Title/Abstract]) AND ("discovered"[Title/Abstract] OR "identified"[Title/Abstract] OR "causative"[Title/Abstract] OR "underlie"[Title/Abstract] OR "basis"[Title/Abstract]) AND "English"[Language] AND ("disease"[Title/Abstract] OR "disorder"[Title/Abstract] OR "syndrome"[Title/Abstract] OR "condition*"[Title/Abstract]) AND ("journal article"[Publication Type] OR "letter"[Publication Type] OR "Case Reports"[Publication Type]) NOT "review"[Publication Type])')

  # Clean up any unnecessary slashes from the query
  query <- gsub("  ", " ", query)  # Remove double spaces
  print(query)
  out_file <- paste0(base_directory, "new_loci")

  # Include a separator ("/") between base_directory and gene_name
  # Modify dest_file_prefix to include the full file path
  out.A <- batch_pubmed_download(pubmed_query_string = query,
                                 format = "medline",
                                 batch_size = 10000,
                                 dest_file_prefix = out_file,
                                 encoding = "ASCII")
  # the function adds 01.txt so, gotta fix that here
  out_file <- paste0(base_directory, "/new_loci", "01.txt")
  print(out_file)

  # Check if the file was created successfully
  cat("Full file path:", out_file, "\n")
  if (file.exists(out_file)) {
    cat("File exists.\n")
    file_path <- out_file
  } else {
    cat("Error: File not found -", out_file, "\n")
  }


  return(file_path)
}



perform_new_pubmed_query()
#
# ### Let's get all the citations to run manubot on
extract_citations <- function(column) {
  column %>%
    str_split(",") %>%                  # Split on commas
    unlist() %>%                        # Flatten list to vector
    str_trim() %>%                      # Trim whitespace
    str_subset("^@\\S+") %>%            # Ensure only terms starting with @ are included
    str_remove("^@") %>%                # Remove the @ symbol from each term
    unique()                            # Remove duplicates
}

# Extract citations from `references` and `additional_literature`
citations_references <- extract_citations(data$references)
citations_additional <- extract_citations(data$additional_literature)

# Combine citations into a single list and remove duplicates
all_citations <- unique(c(citations_references, citations_additional))

myfile = toJSON(all_citations)

write(myfile, args[3])
#
# bash_script <- file.path("/Users/quinlan/Documents/Git/STRchive/data/literature/run_manubot.sh")  # Adjust path if script is in a different directory
# # Run the script using its full path
# citation_chunks <- split(all_citations, ceiling(seq_along(all_citations) / 50)) # Process 50 citations at a time
#
# # Initialize a list to collect all JSON results
# all_json_results <- list()
#
# # Process each chunk
# for (chunk in citation_chunks) {
#   # Convert the chunk into a space-separated string
#   citation_string <- paste(chunk, collapse = " ")
#
#   # Run the bash script and capture output
#   json_output <- system(paste(bash_script, citation_string), intern = TRUE)
#
#   # Parse JSON output in R
#   json_results <- jsonlite::fromJSON(paste(json_output, collapse = "\n"))
#
#   # Append results to the list
#   all_json_results <- c(all_json_results, list(json_results))
# }
#
# # Combine all JSON results into one consolidated list
# consolidated_json <- do.call(c, all_json_results)
