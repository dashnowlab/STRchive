/** years old before not "new" anymore */
export const newThreshold = 2;

/** derive/compute some props from existing props on datum  */
export const deriveDatum = (d) => ({
  ...d,
  position_hg38: `${d.chrom}:${d.start_hg38}-${d.stop_hg38}`,
  locus_tags: [
    ...(d.locus_tags ?? []),
    d.details?.match(/conflict/i) && "conflicting",
    new Date().getFullYear() - d.year <= newThreshold && "new",
  ].filter(Boolean),
});
