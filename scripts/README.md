# Script details

- `Snakefile`: runs all scripts below using the above `snakemake` command
- `check-loci.py`: checks and repairs errors in the STRchive loci JSON file
- `make-catalog.py`: creates TRGT, atarva, LongTR, STRkit, ExpansionHunter, STRaglr, STRanger, UCSC BED16, or extended BED catalogs from all STRchive loci
- `get-literature.py`: searches pubmed for relevent literature related to TR disease loci
- `run-manubot.py`: fetches metadata for all citations into `STRchive-citations.json`
  - `pmid:`, `genereviews:` and `omim:` IDs are fetched in batches from NCBI E-utilities (`citation_utils.py`). GeneReviews chapters are looked up in PubMed, falling back to NCBI Bookshelf for accessions not indexed there. OMIM titles come from NCBI's copy of OMIM, since omim.org blocks automated requests. Takes about a minute to run (on a few thousand citations).
  - `mondo:` IDs are looked up in the EBI Ontology Lookup Service API (`citation_utils.py`).
  - `gnomad:` IDs are generated locally from the ID (no network lookup).
  - All other types (DOIs, Orphanet, URLs, ...) are looked up one at a time with Manubot.
  - With `--cache <previous citations JSON>`, only new or previously failed input IDs are looked up. Other cached entries are kept, even if they're not in the input. Cached PubMed, GeneReviews, OMIM and MONDO entries that were made by Manubot are looked up again once. Used for pull requests (`old-refs`) and new literature (`new-refs`).
  - With `--cache` and `--refresh`, all citations are looked up again, and the previous entry is only used where a lookup fails (e.g. a site is down). This is used by the monthly scheduled run (`all`).
  - A lookup counts as failed if the title is empty or a known bot-check/placeholder page (e.g. "Checking your browser - reCAPTCHA"); failed entries are retried on the next run if their IDs are in the input.
  - Optional environment variables: `NCBI_API_KEY` (raises the NCBI rate limit from 3 to 10 requests/second) and `NCBI_EMAIL` (contact address NCBI asks automated clients to provide).
- `citation_utils.py`: NCBI E-utilities and EBI Ontology Lookup Service lookups used by `run-manubot.py`