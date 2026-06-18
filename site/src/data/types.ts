import type curations from "~/criTRia-curations.json";
import type loci from "~/STRchive-loci.json";

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

export type Curations = typeof curations;
export type Curation = Curations[number];
