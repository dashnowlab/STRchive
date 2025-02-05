import { sortBy, uniq } from "lodash-es";

/** years old before not "new" anymore */
export const newThreshold = 2;

/** derive/compute some props from existing props on locus  */
export const deriveLocus = (locus, loci, citations) => {
  /** keep existing data */
  locus = { ...locus };

  /** nice looking id, for labels */
  locus.nice_id = locus.id.replace("_", " ");

  /** construct full position strings */
  for (const assembly of ["hg19", "hg38", "t2t"])
    locus[`position_base0_${assembly}`] =
      `${locus.chrom}:${locus[`start_${assembly}`]}-${locus[`stop_${assembly}`]}`;
  for (const assembly of ["hg19", "hg38", "t2t"])
    locus[`position_base1_${assembly}`] =
      `${locus.chrom}:${locus[`start_${assembly}`] + 1}-${locus[`stop_${assembly}`]}`;

  /** make prevalence fraction */
  if (locus.prevalence)
    locus.prevalence = locus.prevalence
      .split("/")
      .map(Number)
      .map((value) => (value || 0).toLocaleString());

  /** init tags */
  locus.locus_tags ??= [];

  /** add tags */
  if (new Date().getFullYear() - parseInt(locus.year) <= newThreshold)
    locus.locus_tags.push("new");

  /** clean up tags */
  locus.locus_tags = uniq(locus.locus_tags.filter(Boolean));

  if (loci) {
    /** find other loci that are the same gene */
    locus.genes = loci.filter(({ gene }) => gene === locus.gene);
  }

  if (citations) {
    /** map for quick lookup of citation by id */
    const citationLookup = Object.fromEntries(
      citations.map((citation) => [citation.id, { ...citation }]),
    );

    /** incrementing counter for order of citations */
    let number = 0;

    /** fields that may contain in-text citations, in order they appear on locus page (so number increases as you go down page) */
    const inTextFields = [
      "disease_description",
      "prevalence_details",
      "age_onset",
      "details",
      "mechanism_detail",
      "year",
    ];

    /** extract in-text citations */
    for (const key of inTextFields)
      if (locus[key])
        locus[key] = extractCitations(locus[key]).map(
          ({ text, references }) => ({
            text,
            references: references.map((id) => {
              const citation = citationLookup[id];
              if (citation && !citation.number) citation.number = ++number;
              return { id, ...(citation || {}) };
            }),
          }),
        );

    /** look up citation by id */
    const lookupCitation = (id) => ({ id, ...(citationLookup[id] ?? {}) });

    /** fill in citation details, and sort by number (if any) */
    const mapReferences = (refs) =>
      sortBy((refs || []).map(lookupCitation), "number");

    /** map references */
    locus.references = mapReferences(locus.references || []);

    /** map additional literature */
    locus.additional_literature = mapReferences(
      locus.additional_literature || [],
    );
  }

  return locus;
};

/** split field that may contain in-text references into parts */
const extractCitations = (text) => {
  if (typeof text !== "string") return text;

  /** find sub-string indices in text that match citation format */
  const indices = uniq([
    0,
    ...[...(text.matchAll(/\[@\w+:.+?\]/g) || [])]
      .map((match) => [match.index, match.index + match[0].length])
      .flat(),
    text.length,
  ]);

  return (
    Array(indices.length - 1)
      .fill(null)
      /** extract sub-strings from pairs of indices */
      .map((_, index) =>
        text.substring(indices[index], indices[index + 1]).trim(),
      )
      .map((text) => {
        /** if reference */
        if (text.startsWith("["))
          return {
            text: "",
            references: text
              /** widdle down to just reference id */
              .replaceAll("@", "")
              .replaceAll("[", "")
              .replaceAll("]", "")
              /** split multiple */
              .split(";")
              .map((reference) => reference.trim()),
          };
        /** else, return plain text */ else return { text, references: [] };
      })
  );
};
