import { mapValues, uniq } from "lodash-es";

/** derive/compute some props from existing props on datum  */
export const deriveDatum = (d) => ({
  ...d,
  position_hg38: `${d.chrom}:${d.start_hg38}-${d.stop_hg38}`,
  locus_tags: [
    ...(d.locus_tags ?? []),
    d.details?.match(/conflict/i) && "conflicting",
    new Date().getFullYear() - d.year < 3 && "new",
  ].filter(Boolean),
});

/** get unique values */
export const getUnique = (data) => uniq(data).filter(Boolean);
