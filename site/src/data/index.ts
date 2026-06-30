import { parse } from "@/util/markdown";
import { cloneDeep, sortBy, uniq } from "lodash-es";
import rawCurations from "~/criTRia-curations.json";
import rawCitations from "~/STRchive-citations.json";
import rawLoci from "~/STRchive-loci.json";

/** years old before not "new" anymore */
export const newThreshold = 2;

/** map for quick lookup of citation by id */
const citationLookup: Record<string, Citation> = Object.fromEntries(
  (rawCitations as Citations).map((citation) => [citation.id, citation]),
);

/** look up citation by id */
const lookupCitation = (id: string) => ({
  id,
  ...(citationLookup[id] ?? {}),
});

/** fill in citation details, and sort by number (if any) */
const mapReferences = (refs: string[], number = false) =>
  sortBy(
    (refs || [])
      .map(lookupCitation)
      .map((citation, index) =>
        number ? { ...citation, number: index + 1 } : citation,
      ),
    "number",
  );

/** split field that may contain in-text references into parts */
const extractCitations = (text: string | null | undefined, list: string[]) => {
  if (typeof text !== "string") return text;

  text = parse(text);

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
          if (!list.includes(id)) list.push(id);
          const number = list.indexOf(id) + 1;
          return { number, id, ...(citation ?? {}) };
        }),
      }))
  );
};

/** raw loci data, plus extra derived/computed props */
export const loci = rawLoci.map((locus) => {
  /** avoid mutation bugs */
  locus = cloneDeep(locus);

  /** get position string */
  const position = (assembly: "hg19" | "hg38" | "t2t", base = 0) =>
    `${locus.chrom}:${locus[`start_${assembly}`] + base}-${locus[`stop_${assembly}`]}`;

  /** incrementing reference list */
  const references: string[] = [];

  return {
    /** keep existing data */
    ...locus,

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
    genes: rawLoci.filter(({ gene }) => gene === locus.gene),

    /** process fields with in-text citations */
    disease_description: extractCitations(
      locus.disease_description,
      references,
    ),
    prevalence_details: extractCitations(locus.prevalence_details, references),
    age_onset: extractCitations(locus.age_onset, references),
    details: extractCitations(locus.details, references),
    mechanism_detail: extractCitations(locus.mechanism_detail, references),
    detection: extractCitations(locus.detection, references),
    year: extractCitations(locus.year, references),

    /** map references */
    references: mapReferences(references, true),

    /** map additional literature */
    additional_literature: mapReferences(locus.additional_literature),
  };
});

/** raw curation data, plus extra derived/computed props */
export const curations = rawCurations.map((curation) => {
  /** avoid mutation */
  curation = cloneDeep(curation);

  /** incrementing reference list */
  const references: string[] = [];

  return {
    /** keep existing data */
    ...curation,

    /** process fields with in-text citations */
    Description: extractCitations(curation.Description, references),
    genetic_evidence_details: curation.genetic_evidence_details.map(
      (evidence) => ({
        ...evidence,
        "Evidence detail": extractCitations(
          evidence["Evidence detail"],
          references,
        ),
      }),
    ),
    experimental_evidence_details: curation.experimental_evidence_details.map(
      (evidence) => ({
        ...evidence,
        "Evidence detail": extractCitations(
          evidence["Evidence detail"],
          references,
        ),
      }),
    ),

    /** map references */
    references: mapReferences(references, true),
  };
});

export type Loci = typeof loci;
export type Locus = Loci[number];

/** citations file too big for ts to infer type */
export type Citation = {
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
export type Citations = Citation[];

export type Curations = typeof curations;
export type Curation = Curations[number];
