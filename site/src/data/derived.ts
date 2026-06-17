import type loci from "~/STRchive-loci.json";
import { cloneDeep, sortBy, uniq } from "lodash-es";

type Loci = typeof loci;
type Locus = Loci[number];
/** citations file too big for ts to infer type */
type Citation = {
  id?: string;
  link?: string;
  title?: string;
  type?: string;
  doi?: string;
  authors?: string[][];
  publisher?: string;
  issn?: string;
  date?: string;
  abstract?: string;
  language?: string;
  note?: string;
};

/** years old before not "new" anymore */
export const newThreshold = 2;

/** derive/compute some props from existing props on locus  */
export const deriveLocus = (
  locus: Locus,
  loci: Loci,
  citations?: Citation[],
) => {
  /** get position string */
  const position = (assembly: "hg19" | "hg38" | "t2t", base = 0) =>
    `${locus.chrom}:${locus[`start_${assembly}`] + base}-${locus[`stop_${assembly}`]}`;

  /** map for quick lookup of citation by id */
  const citationLookup = Object.fromEntries(
    citations?.map((citation) => [citation.id, { ...citation }]) ?? [],
  );

  /** look up citation by id */
  const lookupCitation = (id: string) => ({
    id,
    ...(citationLookup[id] ?? {}),
  });

  /** fill in citation details, and sort by number (if any) */
  const mapReferences = (refs: string[]) =>
    sortBy((refs || []).map(lookupCitation), "number");

  /** incrementing counter for order of citations */
  let citationNumber = 0;

  /** split field that may contain in-text references into parts */
  const extractCitations = (text: string | null | undefined) => {
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
        .map(({ text, references }) => ({
          text,
          references: references.map((id) => {
            const citation = citationLookup[id];
            if (citation && !citation.number)
              citation.number = ++citationNumber;
            return { id, ...(citation || {}) };
          }),
        }))
    );
  };

  return {
    /** keep existing data */
    ...cloneDeep(locus),

    /** nice looking id, for labels */
    nice_id: locus.id.replace("_", " "),

    /** construct full position strings */
    position_base0_hg19: position("hg19"),
    position_base0_hg38: position("hg38"),
    position_base0_t2t: position("t2t"),
    position_base1_hg19: position("hg19", 1),
    position_base1_hg38: position("hg38", 1),
    position_base1_t2t: position("t2t", 1),

    /** make prevalence fraction */
    prevalence: locus.prevalence
      ?.split("/")
      .map(Number)
      .map((value) => (value || 0).toLocaleString()),

    /** tags */
    locus_tags: uniq([
      ...locus.locus_tags,
      locus.year &&
      new Date().getFullYear() - parseInt(locus.year) <= newThreshold
        ? "new"
        : "",
    ]).filter(Boolean),

    /** find other loci that are the same gene */
    genes: loci.filter(({ gene }) => gene === locus.gene),

    /** process fields with in-text citations */
    disease_description: extractCitations(locus.disease_description),
    prevalence_details: extractCitations(locus.prevalence_details),
    age_onset: extractCitations(locus.age_onset),
    details: extractCitations(locus.details),
    mechanism_detail: extractCitations(locus.mechanism_detail),
    year: extractCitations(locus.year),

    /** map references */
    references: mapReferences(locus.references || []),

    /** map additional literature */
    additional_literature: mapReferences(locus.additional_literature || []),
  };
};
