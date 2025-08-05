# STRchive

Short Tandem Repeats (STRs) are a type of genetic variation that are associated with many rare diseases. Information about pathogenic STRs is often out-of-date and scattered across different databases, making it difficult to find and interpret STR variants. STRchive ("ess tee archive") aims to solve this problem by providing a central community resource.

[⭐️ View the data at strchive.org ⭐️](http://strchive.org/)

If you use STRchive in your research, please cite:
Hiatt, L., Weisburd, B., Dolzhenko, E., Rubinetti, V., Avvaru, A.K., VanNoy, G.E., Kurtas, N.E., Rehm, H.L., Quinlan, A. and Dashnow, H.✉, 2025. STRchive: a dynamic resource detailing population-level and locus-specific insights at tandem repeat disease loci. Genome medicine doi: [https://doi.org/10.1186/s13073-025-01454-4](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-025-01454-4).

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="http://strchive.org/">STRchive</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://github.com/hdashnow">Harriet Dashnow</a> is licensed under <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>

## Contributors

- Harriet Dashnow
- Laurel Hiatt
- Akshay Avvaru
- Vincent Rubinetti
- Macayla Weiner

## Contributing

If you notice an error, omission, or update, feel free to leave a comment or create a pull request.

To make a change to the STRchive data itself, please edit `data/STRchive-loci.json`

Then run the "linting" script and fix any errors:  
`python scripts/check-loci.py data/STRchive-loci.json`

## Development

### Run all scripts to update STRchive

From the root directory, run:  
`snakemake`

Or to skip retrieve and manubot stages, which will speed things up substantially:  
`snakemake --config stages="skip-refs"`

### Update TRGT genotyping catalogs

```
python scripts/make-catalog.py -g hg38 -f TRGT data/STRchive-loci.json data/STRchive-disease-loci.hg38.TRGT.bed
python scripts/make-catalog.py -g T2T -f TRGT data/STRchive-loci.json data/STRchive-disease-loci.T2T-chm13.TRGT.bed
python scripts/make-catalog.py -g hg19 -f TRGT data/STRchive-loci.json data/STRchive-disease-loci.hg19.TRGT.bed
```

### Update extended BED files

```
python scripts/make-catalog.py -f bed -g hg38 data/STRchive-loci.json data/STRchive-disease-loci.hg38.bed
python scripts/make-catalog.py -f bed -g T2T data/STRchive-loci.json data/STRchive-disease-loci.T2T-chm13.bed
python scripts/make-catalog.py -f bed -g hg19 data/STRchive-loci.json data/STRchive-disease-loci.hg19.bed
```

### Install dependencies

New install:  
```
conda env create --file scripts/environment.yml
conda activate strchive
```

Update existing installation:  
```
conda activate strchive
conda env update --file scripts/environment.yml --prune
conda activate strchive
```


Note: biomaRt isn't playing nicely with conda, so installing it within the R script where it is used.
