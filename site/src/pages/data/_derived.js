/** derive/compute some props from existing props on datum  */
export const deriveDatum = (d) => ({
  ...d,
  position_hg38: `${d.chrom}:${d.start_hg38}-${d.stop_hg38}`,
  tags: [
    ...(d.tags ?? []),
    d.details?.match(/conflict/i) && "conflicting",
    new Date().getFullYear() - d.Year < 3 && "new",
  ].filter(Boolean),
});

/** get unique inheritance values */
export const getUniqueInheritance = (data) =>
  [...new Set(data.map((d) => d.Inheritance))].filter(Boolean);
