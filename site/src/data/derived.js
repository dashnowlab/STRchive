import { uniq } from "lodash-es";

/** years old before not "new" anymore */
export const newThreshold = 2;

/** derive/compute some props from existing props on datum  */
export const deriveDatum = (d) => {
  /** keep existing data */
  d = { ...d };

  /** construct full position strings */
  for (const assembly of ["hg19", "hg38", "t2t"])
    d[`position_${assembly}`] =
      `${d.chrom}:${d[`start_${assembly}`]}-${d[`stop_${assembly}`]}`;

  /** make prevalence fraction */
  if (d.prevalence)
    d.prevalence = d.prevalence
      .split("/")
      .map(Number)
      .map((value) => (value || 0).toLocaleString());

  /** init tags */
  d.locus_tags ??= [];

  /** add tags */
  if (new Date().getFullYear() - parseInt(d.year) <= newThreshold)
    d.locus_tags.push("new");

  /** clean up tags */
  d.locus_tags = uniq(d.locus_tags.filter(Boolean));

  return d;
};
