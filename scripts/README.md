Note: all commands are assumed to be run from this directory (`scripts`)

# Install dependencies

New install:  
`conda env create --file environment.yml`

Update existing installation:  
```
conda activate strchive
conda env update --file environment.yml --prune
```

## Install biomaRt

I can't figure out how to get biomaRt to play nicely with conda, so installing it manually within the R environment:

```
$ conda activate strchive

$ R

> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install("biomaRt")
```

# Run all scripts to update STRchive

From the `scripts` directory, run `snakemake`

# Script details

- `Snakefile`: runs all scripts below using the above `snakemake` command
- `check-loci.py`: checks and repairs errors in the STRchive loci JSON file
- `make-catalog.py`: creates TRGT catalog or general extended BED format files from all STRchive loci
- `get-literature.R`: searches pubmed for relevent literature related to TR disease loci
- `run-manubot.py`: fetched metadata for all citations (currently takes several hours)