import { FaBookMedical, FaLaptopMedical } from "react-icons/fa6";

/** get identifiers for datum */
export const getIdentifiers = (d) =>
  /** list of identifier types */
  [
    {
      name: "Clinical Links",
      Icon: FaBookMedical,
      identifiers: [
        {
          key: "gard",
          name: "GARD",
          link: "https://rarediseases.info.nih.gov/diseases/$ID/index",
          tooltip: "Easy-to-understand rare disease information and resources",
          info: "https://rarediseases.info.nih.gov",
        },
        {
          key: "genereviews",
          name: "GeneReviews",
          link: "https://www.ncbi.nlm.nih.gov/books/$ID",
          tooltip:
            "Point-of-care clinical resource for diagnosis, management, and genetic counseling",
          info: "https://www.ncbi.nlm.nih.gov/books/NBK1116",
        },
        {
          key: "malacard",
          name: "MalaCard",
          link: "https://www.malacards.org/card/$ID",
          tooltip:
            "Detailed disease database with integrated data from 75 web sources",
          info: "https://www.malacards.org",
        },
        {
          key: "medgen",
          name: "MedGen",
          link: "https://www.ncbi.nlm.nih.gov/medgen/?term=$ID",
          tooltip:
            "Aggregated data on conditions and phenotypes within medical genetics",
          info: "https://www.ncbi.nlm.nih.gov/medgen",
        },

        {
          key: "mondo",
          name: "Mondo",
          link: "https://purl.obolibrary.org/obo/MONDO_$ID",
          tooltip:
            "Disease definitions, ontologies, and harmonized terminologies",
          info: "http://obofoundry.org/ontology/mondo.html",
        },
        {
          key: "omim",
          name: "OMIM",
          link: "https://omim.org/entry/$ID",
          tooltip:
            "Compendium focused on genotype-phenotype relationships in Mendelian disorders",
          info: "https://omim.org",
        },
        {
          key: "orphanet",
          name: "Orphanet",
          link: "https://www.orpha.net/en/disease/detail/$ID",
          tooltip: "High-quality information on rare diseases",
          info: "https://www.orpha.net",
        },
      ],
    },
    {
      name: "Bioinformatical Links",
      Icon: FaLaptopMedical,
      identifiers: [
        {
          key: "gnomad",
          name: "gnomAD",
          link: "https://gnomad.broadinstitute.org/short-tandem-repeat/$ID?dataset=gnomad_r4",
          tooltip:
            "TR genotype distributions from 18.5k individuals with ancestry, sex, and read visualization data",
          info: "https://gnomad.broadinstitute.org/news/2022-01-the-addition-of-short-tandem-repeat-calls-to-gnomad",
        },
        {
          key: "gnomad",
          name: "TR Explorer",
          link: "https://trexplorer.broadinstitute.org/#sc=isPathogenic&sd=DESC&showRs=1&searchQuery=$ID&showColumns=0i1i2i3i4i7i21i17",
          tooltip:
            "An online portal for exploring genome-wide tandem repeat (TR) catalogs",
          info: "https://trexplorer.broadinstitute.org/",
        },
        {
          key: "stripy",
          name: "STRipy",
          link: "https://stripy.org/database/$ID",
          tooltip:
            "Pathogenic STR database with population-wide repeat frequencies",
          info: "https://stripy.org",
        },
        {
          key: "tr_atlas",
          name: "TR-Atlas",
          link: "https://wlcb.oit.uci.edu/TRatlas/Distribution.php?index_id=$ID",
          tooltip:
            "Biobank-scale reference of 0.86 million TRs derived from 338,963 diverse samples",
          info: "https://wlcb.oit.uci.edu/TRatlas",
        },

        {
          key: "webstr_hg38",
          name: "WebSTR hg38",
          link: "https://webstr.ucsd.edu/locus?genome=hg38&repeat_id=$ID",
          tooltip:
            "Genome-wide STR data across populations, with details on variation, expression, constraint, and imputation metrics",
          info: "https://webstr.ucsd.edu",
        },
        {
          key: "webstr_hg19",
          name: "WebSTR hg19",
          link: "https://webstr.ucsd.edu/locus?repeat_id=$ID",
          tooltip:
            "Genome-wide STR data across populations, with details on variation, expression, constraint, and imputation metrics",
          info: "https://webstr.ucsd.edu",
        },
      ],
    },
  ]
    .map(({ identifiers, ...rest }) => ({
      ...rest,
      identifiers: identifiers
        /** for each reference type */
        .map(({ key, link, ...rest }) => ({
          ...rest,
          links:
            /** access key field, normalize to array */
            [d[key] ?? []]
              .flat()
              /** normalize to string */
              .map((id) => String(id))
              /** make links */
              .map((id) => ({
                /** link text */
                label: id,
                /** link target, with id inserted */
                link: link.replace("$ID", id),
              })) || [],
        }))
        /** remove identifiers with no links */
        .filter(({ links }) => links?.length),
    }))
    /** remove types with no identifiers */
    .filter(({ identifiers }) => identifiers.length);
