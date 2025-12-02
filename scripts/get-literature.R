# Usage:
# Rscript get-literature.R STRchive-loci.json literature-directory out-citations.json STRchive-loci-literature.json

# More useful error messages
options(error=traceback)

## Set default CRAN mirror
local({r <- getOption("repos")
       r["CRAN"] <- "https://cran.r-project.org"
       options(repos=r)
})

# Install biomaRt if not already installed
tryCatch({
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", ask = FALSE, update = FALSE)
    }
    BiocManager::install("biomaRt", ask = FALSE)
  }
}, error = function(e) {
  cat("Error installing biomaRt package.\n", file=stderr())
  quit(status = 1)
})

tryCatch({
suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(biomaRt)
  library(rentrez)
  library(stringr)
  library(purrr)
})
}, error = function(e) {
  cat("Error loading required packages.\n", file=stderr())
  quit(status = 1)
})

### Data Setup
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Need input json, base directory, output json for all citations and output json for citations per locus")
}

cat("Arguments: ", args, "\n", file=stderr())

# Manually downloaded from biomart
synonym_file = 'data/literature/gene_synonyms.tsv' # should make an arg with default?

data <- fromJSON(args[1])

# get unique references used in the json for each entry
data <- data %>%
  mutate(references = apply(data[, c("age_onset", "mechanism_detail", "details",
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
cat("Searching for genes:", gene_list, "\n", file=stderr())

# Filter out NA values, if there are any, from the 'gene' column
filtered_data <- data[!is.na(data$gene), ]

# Extract STRchive gene names into a list
gene_list <- as.character(unique(filtered_data$gene))

### Use biomaRt to get gene synonyms
# This is necessary because gene names have evolved over time and can differ
# by publication
hosts <- c("https://useast.ensembl.org", "https://asia.ensembl.org", "https://www.ensembl.org")
for (host in hosts) {
  tryCatch({
    mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
                    host = host, verbose = TRUE)
    break
  }, error = function(e) {
    cat("Failed to create the mart object for host:", host, ". Trying another.\n", file=stderr())
  })
}

# Retrieve the HGNC symbol and Gene Synonym for the genes in the list
# Check if the mart object is created
if (exists("mart")) {
  # Retrieve the HGNC symbol and Gene Synonym for the genes in the list
  gene_info <- getBM(attributes = c("hgnc_symbol", "external_synonym"),
                     filters = "hgnc_symbol",
                     values = gene_list,
                     mart = mart)

} else {
  cat(paste0("Failed to create the mart object. Using previous gene synonyms from file: ", synonym_file), file=stderr())
  gene_info = read.csv(synonym_file, sep = '\t')
  gene_info = merge(gene_info, data.frame(hgnc_symbol = gene_list, external_synonym = NA), all = TRUE)
  gene_info = unique(gene_info)
}

# Sort gene_info by 'hgnc_symbol' to ensure consistency
gene_info <- gene_info[order(gene_info$hgnc_symbol), ]

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
if ("FMR1" %in% gene_list) {
  gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"FMR1\"", external_synonym = "\"FMR-1\""))
}

#because pubmed hates slashes and these use slashes
if ("NUTM2B-AS1" %in% gene_list) {
  gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"NUTM2B-AS1\"", external_synonym = "\"LOC642361/NUTM2B-AS1\""))
}
if ("ATXN10" %in% gene_list) {
  gene_info <- rbind(gene_info, data.frame(hgnc_symbol = "\"ATXN10\"", external_synonym = "\"SCA10\""))
}

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

# Optional: set NCBI email from environment (recommended by NCBI)
# if (nzchar(Sys.getenv("ENTREZ_EMAIL"))) {
#   set_entrez_email(Sys.getenv("ENTREZ_EMAIL"))
# } else {
#   cat("Note: Set ENTREZ_EMAIL env var to comply with NCBI policies.\n", file=stderr())
# }

# Helper to fetch MEDLINE text in batches and write to file prefix
fetch_medline_and_write <- function(pmids, outfile_prefix, batch_size = 200) {
  if (length(pmids) == 0) return(FALSE)
  medline_texts <- character(0)
  for (i in seq(1, length(pmids), by = batch_size)) {
    chunk <- pmids[i:min(i + batch_size - 1, length(pmids))]
    txt <- tryCatch({
      entrez_fetch(db = "pubmed", id = chunk, rettype = "medline", retmode = "text")
    }, error = function(e) {
      cat("entrez_fetch error:", conditionMessage(e), "\n", file=stderr())
      return(NULL)
    })
    if (!is.null(txt)) medline_texts <- c(medline_texts, txt)
  }
  if (!length(medline_texts)) return(FALSE)
  out_file <- paste0(outfile_prefix, "_batch_01.txt")
  writeLines(medline_texts, con = out_file, useBytes = TRUE)
  return(file.exists(out_file))
}

# function to perform the pubmed query
# Function printout includes gene name and if there are results, confirms
# that a file has been created
perform_pubmed_query <- function(gene_info) {
  file_paths <- list()  # Initialize the list to store all publications
  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    cat("Processing gene:", gene_name, "\n", file=stderr())
    # Use str_detect to check if gene_name is present in each consolidated string
    idx <- str_detect(consolidated_strings, regex(gene_name, ignore_case = TRUE))

    # Put together the relevant consolidated strings
    or_terms <- paste(consolidated_strings[idx], collapse = ' OR ')
    # Splitting the or_terms into individual terms
    individual_terms <- unlist(strsplit(or_terms, " OR "))

    # Adding [Title/Abstract] to each term
    joined_terms <- paste0(individual_terms, "[Title/Abstract]")
    # Construct the query with organized or_terms
    terms_repeat <- c(
      '"repeat expansion"[Title/Abstract]',
      '"repeat expansions"[Title/Abstract]',
      '"tandem repeat"[Title/Abstract]',
      '"tandem repeats"[Title/Abstract]',
      '"repeat sequence"[Title/Abstract]',
      '"repeat sequences"[Title/Abstract]',
      '"repeat length"[Title/Abstract]',
      '"repeat lengths"[Title/Abstract]',
      '"expansion"[Title]',
      '"expansions"[Title]',
      '"repeats"[Title]'
    )

    terms_gene <- joined_terms  # e.g., "FMR1"[Title/Abstract] OR "FMR-1"[Title/Abstract] ...
    terms_language <- '"English"[Language]'
    terms_disease <- c(
      "disease*[Title/Abstract]",
      "disorder*[Title/Abstract]",
      "syndrome*[Title/Abstract]",
      "patient*[Title/Abstract]",
      "proband*[Title/Abstract]"
    )
    terms_pubtype <- c(
      '"journal article"[Publication Type]',
      '"letter"[Publication Type]',
      '"Case Reports"[Publication Type]'
    )
    terms_exclude <- '"review"[Publication Type]'

    query <- paste(
      "(" , paste(terms_repeat, collapse = " OR "), ")",
      "AND (", paste(terms_gene, collapse = " OR "), ")",
      "AND", terms_language,
      "AND (", paste(terms_disease, collapse = " OR "), ")",
      "AND (", paste(terms_pubtype, collapse = " OR "), ")",
      "NOT", terms_exclude
    )

    # Clean up double spaces
    query <- gsub("\\s{2,}", " ", query)  # Remove double spaces
    gene_name <- gsub('"', '', gene_name)
    out_prefix <- paste0(base_directory, "/", gene_name)

    # print query for debugging
    # cat("PubMed query for gene", gene_name, ":\n", query, "\n")

    # Use rentrez to search and fetch
    tryCatch({
      srch0 <- entrez_search(db = "pubmed", term = query, retmax = 0, use_history = FALSE)
      count <- ifelse(is.null(srch0$count), 0L, as.integer(srch0$count))
      cat("Found", count, "articles for gene:", gene_name, "\n", file=stderr())
      if (count == 0L) {
        cat("Skipping fetch for:", gene_name, "\n", file=stderr())
        next
      }
      retmax <- min(count, 10000L)
      srch <- entrez_search(db = "pubmed", term = query, retmax = retmax)
      ok <- fetch_medline_and_write(srch$ids, out_prefix)
      if (!ok) stop("Failed to write MEDLINE file")
    }, error = function(e) {
      cat("batch pubmed download error: ", conditionMessage(e), "\n", file=stderr())
      quit(status = 1)
    })

    # output file name
    out_file <- paste0(out_prefix, "_batch_01.txt")
    cat(out_file, "\n", file=stderr())

    # Check if the file was created successfully
    cat("Full file path:", out_file, "\n", file=stderr())
    if (file.exists(out_file)) {
      file_paths[[gene_name]] <- out_file
    } else {
      cat("Error: File not found -", out_file, "\n", file=stderr())
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
  cat("Reading file:", file_path, "\n", file=stderr())

  tryCatch({
    # Read the file into a character vector
    current_publications <- readLines(file_path)

    # Append to the overall list
    all_publications[[gene_name]] <- current_publications
  }, error = function(e) {
    cat("Warning: error reading file:", file_path, "\n", file=stderr())
    # Append an empty character vector in case of an error
    all_publications[[gene_name]] <- character(0)
  })
}

# empty list for the publication info
pub_info_list <- list()

extract_citation_info <- function(medline_data_list, gene_name) {
  # Combine the MEDLINE records into a single string
  medline_string <- paste(medline_data_list, collapse = "")

  # Extract PMIDs, dates, and titles
  pmids <- str_extract_all(medline_string, "(?<=PMID- )\\d+")[[1]]
  publication_dates <- str_extract_all(medline_string, "(?<=DP  - )\\d+")[[1]]
  title <- str_extract_all(medline_string, "(?<=TI  - ).+?(?=\\.|\\?)")[[1]]

  # Align vector lengths
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

  # Data frame output
  data.frame(
    gene = rep(gene_name, length(pmids)),
    PMID = pmids,
    PublicationDate = publication_dates,
    Title = title,
    stringsAsFactors = FALSE
  )
}

for (gene_name in names(all_publications)) {
  medline_data_list <- all_publications[[gene_name]]
  if (length(medline_data_list) == 0) next  # skip empty files
  pub_info_df <- extract_citation_info(medline_data_list, gene_name)
  pub_info_list[[gene_name]] <- pub_info_df
}

# Combine all the dataframes into a single dataframe
if (length(pub_info_list) == 0) {
  all_pub_info_df <- data.frame(
    gene = character(),
    PMID = character(),
    PublicationDate = character(),
    Title = character(),
    stringsAsFactors = FALSE
  )
} else {
  all_pub_info_df <- do.call(rbind, pub_info_list)
}

# add pubmed search results to literature field for entry
data <- data %>%
  mutate(additional_literature = map_chr(gene, function(g) {
    # Guard if no publications parsed
    if (nrow(all_pub_info_df) == 0) return("")
    matching_pmids <- all_pub_info_df %>%
      filter(gene == g) %>%
      pull(PMID) %>%
      unique() %>%
      paste0("@pmid:", ., collapse = ",")
    if (length(matching_pmids) == 0) "" else matching_pmids
  }))

# remove any redundant pmids from additional_literature that are in references
data <- data %>%
  mutate(additional_literature = mapply(function(lit, refs) {
    # Split the strings by commas to create lists
    lit_list <- unique(unlist(str_split(lit, ",\\s*")))  # Split and remove extra spaces
    refs_list <- unique(unlist(str_split(refs, ",\\s*")))  # Split and remove extra spaces

    # Find common elements between additional_literature and references
    common_elements <- intersect(lit_list, refs_list)

    # Print common elements being removed (if any)
    if (length(common_elements) > 0) {
      cat("Removing the following elements from additional_literature: ", paste(common_elements, collapse = ", "), "\n", file=stderr())
    }

    # Remove common elements from additional_literature
    lit_list <- setdiff(lit_list, common_elements)

    # Join the remaining elements back into a string
    cleaned_lit <- paste(lit_list, collapse = ",")

    # Return the cleaned additional_literature
    return(cleaned_lit)
  }, data$additional_literature, data$references))

# Write out a json with just the locus id and the two literature fields
lit_data <- data[, c("id", "additional_literature", "references")]
write_json(lit_data, args[4])

#### New locus lit retrieval
#new locus query found from reviewing pertinent terms in discovery papers
perform_new_pubmed_query <- function() {
  file_path <- list()  # Initialize the list to store all publications

  terms_repeat <- c(
    '"repeat expansion"[Title/Abstract]',
    '"tandem repeat"[Title/Abstract]'
  )
  terms_discovery <- c(
    '"discovered"[Title/Abstract]',
    '"identified"[Title/Abstract]',
    '"causative"[Title/Abstract]',
    '"underlie"[Title/Abstract]',
    '"basis"[Title/Abstract]'
  )
  terms_language <- '"English"[Language]'
  terms_disease <- c(
    '"disease"[Title/Abstract]',
    '"disorder"[Title/Abstract]',
    '"syndrome"[Title/Abstract]',
    '"condition*"[Title/Abstract]'
  )
  terms_pubtype <- c(
    '"journal article"[Publication Type]',
    '"letter"[Publication Type]',
    '"Case Reports"[Publication Type]'
  )
  terms_exclude <- '"review"[Publication Type]'

  query <- paste(
    "(" , paste(terms_repeat, collapse = " OR "), ")",
    "AND (", paste(terms_discovery, collapse = " OR "), ")",
    "AND", terms_language,
    "AND (", paste(terms_disease, collapse = " OR "), ")",
    "AND (", paste(terms_pubtype, collapse = " OR "), ")",
    "NOT", terms_exclude
  )

  query <- gsub("\\s{2,}", " ", query)
  out_prefix <- paste0(base_directory, "/new_loci")

  tryCatch({
    srch0 <- entrez_search(db = "pubmed", term = query, retmax = 0, use_history = FALSE)
    count <- ifelse(is.null(srch0$count), 0L, as.integer(srch0$count))
    if (count == 0L) {
      cat("Skipping fetch for new loci - no articles found.\n", file=stderr())
      return(NULL)
    }
    retmax <- min(count, 10000L)
    srch <- entrez_search(db = "pubmed", term = query, retmax = retmax)
    ok <- fetch_medline_and_write(srch$ids, out_prefix)
    if (!ok) stop("Failed to write MEDLINE file")
  }, error = function(e) {
    cat("batch pubmed download error: ", conditionMessage(e), "\n", file=stderr())
    quit(status = 1)
  })

  out_file <- paste0(out_prefix, "_batch_01.txt")
  cat("Full file path:", out_file, "\n", file=stderr())
  if (file.exists(out_file)) {
    file_path <- out_file
  } else {
    cat("Error: File not found -", out_file, "\n", file=stderr())
  }

  return(file_path)
}

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

cat(paste0("Wrote all citations to ", args[3], "\n"), file=stderr())
quit(status = 0)
