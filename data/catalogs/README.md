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

### Why Catalog Counts Differ

Most catalogs contain one record per STRchive locus. Atarva, however, emits an
additional BED record for each flanking repeat, so its number of BED lines can
exceed its number of loci.

STRanger and STRaglr require an allele-size threshold. They retain loci when
either the normal maximum or pathogenic minimum is available, inferring the
missing threshold where necessary. Loci without both values are omitted because
they cannot be genotyped or annotated with a size-based threshold.

### [TRGT](https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md)
- `STRchive-disease-loci.hg38.TRGT.bed`
- `STRchive-disease-loci.hg19.TRGT.bed`
- `STRchive-disease-loci.T2T-chm13.TRGT.bed`

### [STRanger](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/variant_catalog_hg38.json)
- `STRchive-disease-loci.hg38.stranger.json`
- `STRchive-disease-loci.hg19.stranger.json`
- `STRchive-disease-loci.T2T-chm13.stranger.json`

This file is designed to work with the [wf-human-variation workflow](https://github.com/epi2me-labs/wf-human-variation/tree/master). It is modeled after this file: [variant_catalog_hg38.json](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/variant_catalog_hg38.json). It should be used with a matching STRagler bed file where the IDs and start coordinates match.

**WARNING:**  

STRanger requires values for "NormalMax" and "PathologicMin". For some loci these values may be missing from STRchive because they could not be verified from the literature. In cases where one of these values is missing it will be inferred as such:  
PathologicMin = NormalMax + 1  
NormalMax = PathologicMin - 1  

If both values are missing from STRchive the locus will not be included in this file (e.g. where pathogenicity is caused by motif change, not allele size).

### [Atarva](https://github.com/dashnowlab/ATaRVa#region-file)

- `STRchive-disease-loci.hg38.atarva.bed.gz` and `.tbi`
- `STRchive-disease-loci.hg19.atarva.bed.gz` and `.tbi`
- `STRchive-disease-loci.T2T-chm13.atarva.bed.gz` and `.tbi`

Uncompressed BED files are also included as the source files used to create the
bgzip-compressed and tabix-indexed catalogs:
- `STRchive-disease-loci.hg38.atarva.bed`
- `STRchive-disease-loci.hg19.atarva.bed`
- `STRchive-disease-loci.T2T-chm13.atarva.bed`

### [STRaglr](https://github.com/BirolLab/straglr#usage)
- `STRchive-disease-loci.hg38.straglr.bed`
- `STRchive-disease-loci.hg19.straglr.bed`
- `STRchive-disease-loci.T2T-chm13.straglr.bed`

### [LongTR](https://github.com/gymrek-lab/longtr#tr-region-bed-file)
- `STRchive-disease-loci.hg38.longTR.bed`
- `STRchive-disease-loci.hg19.longTR.bed`
- `STRchive-disease-loci.T2T-chm13.longTR.bed`

### [General BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
- `STRchive-disease-loci.hg38.general.bed`
- `STRchive-disease-loci.hg19.general.bed`
- `STRchive-disease-loci.T2T-chm13.general.bed`

### [STRkit](https://github.com/davidlougheed/strkit/blob/master/docs/caller_catalog.md)
- `STRchive-disease-loci.hg38.strkit.bed`
- `STRchive-disease-loci.hg19.strkit.bed`
- `STRchive-disease-loci.T2T-chm13.strkit.bed`

These headerless, 0-based, half-open BED-like files assign each locus a stable
STRchive ID and use the first pathogenic motif for STRkit calling.

### [ExpansionHunter](https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md)
- `STRchive-disease-loci.hg38.expansionhunter.json`
- `STRchive-disease-loci.hg19.expansionhunter.json`
- `STRchive-disease-loci.T2T-chm13.expansionhunter.json`

### [UCSC Genome Browser](https://genome.ucsc.edu/goldenPath/help/bigBed.html)
- `STRchive-disease-loci.hg38.ucsc.bed`
- `STRchive-disease-loci.hg19.ucsc.bed`
- `STRchive-disease-loci.T2T-chm13.ucsc.bed`
- `STRchive-disease-loci.hg38.ucsc.bb`
- `STRchive-disease-loci.hg19.ucsc.bb`
- `STRchive-disease-loci.T2T-chm13.ucsc.bb`

The BED16 files can be imported as UCSC custom tracks. The bigBed files can be
loaded in UCSC and other genome browsers, including IGV.

Item colors encode evidence strength: green for Definitive or Strong, blue for
Moderate, amber for Limited, gray for Provisional, and red for Disputed or
Refuted loci.

The BED file has no column header: `bedToBigBed` accepts only BED records.
UCSC Table Browser adds a `#chrom` header when it exports a track for download.
For direct custom-track loading, an optional `track` line can be prepended:

```text
track name="STRchive" description="STRchive disease-associated STR loci" type=bed itemRgb=On
```

Do not include that `track` line when creating a bigBed; configure its display
settings in the track hub's `trackDb.txt` instead.

Snakemake downloads the matching UCSC chromosome-size files and creates all
three bigBed files with `bedToBigBed -as=strchive.as -type=bed9+7`.
