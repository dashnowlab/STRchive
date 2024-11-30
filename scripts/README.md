# Script details

- `Snakefile`: runs all scripts below using the above `snakemake` command
- `check-loci.py`: checks and repairs errors in the STRchive loci JSON file
- `make-catalog.py`: creates TRGT catalog or general extended BED format files from all STRchive loci
- `get-literature.R`: searches pubmed for relevent literature related to TR disease loci
- `run-manubot.py`: fetched metadata for all citations (currently takes several hours)