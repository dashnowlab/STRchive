import { mapValues } from "lodash-es";

/** derive/compute some props from existing props on datum  */
export const deriveDatum = (d) =>
  mapValues(
    {
      ...d,
      position_hg38: `${d.chrom}:${d.start_hg38}-${d.stop_hg38}`,
      tags: [
        ...(d.tags ?? []),
        d.details?.match(/conflict/i) && "conflicting",
        new Date().getFullYear() - d.Year < 3 && "new",
      ].filter(Boolean),
    },
    (value, key) =>
      [
        "OMIM",
        "Mondo",
        "MedGen",
        "Orphanet",
        "GARD",
        "GeneReviews",
        "gnomAD_gene",
        "STRipy_gene",
        "WebSTR_hg38",
        "WebSTR_hg19",
      ].includes(key) && typeof value === "string"
        ? value.split("; ")
        : value,
  );

/** get unique inheritance values */
export const getUniqueInheritance = (data) =>
  [...new Set(data.map((d) => d.Inheritance))].filter(Boolean);
