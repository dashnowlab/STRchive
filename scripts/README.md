Note: all commands are assumed to be run from this directory (`scripts`)

# Install dependencies

`conda env update --file environment.yml --prune`

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

`snakemake`