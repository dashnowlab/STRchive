# STRchive Catalogs

This directory contains various genotyping and annotation catalogs and based on STRchive tandem repeat loci.

**CAVEATS:**
- Some of these files are still in active development and should be used with care. The specific coordinates and motifs chosen can affect genotyping accuracy.
- Information about the overall pathogenicity of a locus and about specific ranges and motifs is provided as our best estimate. All information should be verified.

## Reference Genomes

- hg38 is the "default" reference from which the others are derived
- hg19
- CHM13-T2T

## Genotypers and File Descriptions

File format:
`STRchive-disease-loci.[reference genome].[software].[file extension(s) e.g. bed, json, bed.gz]`

### TRGT
- `STRchive-disease-loci.hg19.TRGT.bed`

### STRanger
- `STRchive-disease-loci.hg38.stranger.json`

This file is designed to work with the [wf-human-variation workflow](https://github.com/epi2me-labs/wf-human-variation/tree/master). It is modeled after this file: [variant_catalog_hg38.json](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/variant_catalog_hg38.json). It should be used with a matching STRagler bed file where the IDs and start coordinates match.

**WARNING:**  

STRanger requires values for "NormalMax" and "PathologicMin". For some loci these values may be missing from STRchive because they could not be verified from the literature. In cases where one of these values is missing it will be inferred as such:  
PathologicMin = NormalMax + 1  
NormalMax = PathologicMin - 1  

If both values are missing from STRchive the locus will not be included in this file (e.g. where pathogenicity is caused by motif change, not allele size).

### Atarva

- `STRchive-disease-loci.hg38.atarva.bed.gz`
- `STRchive-disease-loci.hg38.atarva.bed.gz.tbi`

Additionally this file is included as the source and to track changes, but is not used by atarva:  
- `STRchive-disease-loci.hg38.atarva.bed`

### Expansion Hunter
- `STRchive-disease-loci.hg19.expansionhunter.json`
