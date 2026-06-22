import { newThreshold } from ".";
import { classifications } from "@/data/curations";
import {
  IconCheck,
  IconExclamationMark,
  IconMinus,
  IconPlus,
  IconPoint,
  IconQuestionMark,
  IconSlash,
  IconSparkles,
  IconX,
} from "@tabler/icons-react";

/** https://tailwindcss.com/docs/colors */

/** order and properties of tags */
export const tagOptions = [
  /** general */
  {
    /** key on locus item where tag can be found */
    key: "locus_tags",
    /** value of tag in locus[key] array */
    value: "new",
    /** human readable label */
    label: "New",
    /** longer description, for tooltips and etc. */
    description: `Less than ~${newThreshold} years old`,
    /** icon */
    Icon: IconSparkles,
    /** classes on icon */
    className: "bg-primary text-white",
    /** whether to show in loci page table */
    filter: true,
  },
  {
    key: "gene_strand",
    value: "+",
    Icon: IconPlus,
    className: "bg-primary text-white",
  },
  {
    key: "gene_strand",
    value: "-",
    Icon: IconMinus,
    className: "bg-secondary text-white",
  },

  /** locus tags */
  {
    key: "locus_tags",
    value: "somatic_instability",
    label: "Som. Inst.",
    description: "Somatic instability",
  },
  {
    key: "locus_tags",
    value: "anticipation",
    label: "Anticip.",
    description: "Anticipation",
  },
  {
    key: "locus_tags",
    value: "contraction",
    label: "Contrac.",
    description: "Contraction",
  },
  {
    key: "locus_tags",
    value: "maternal_expansion",
    label: "Mat. Exp.",
    description: "Maternal expansion",
  },
  {
    key: "locus_tags",
    value: "paternal_expansion",
    label: "Pat. Exp.",
    description: "Paternal expansion",
  },
  {
    key: "locus_tags",
    value: "length_affects_severity",
    label: "Len. → Sev.",
    description: "Length affects severity",
  },
  {
    key: "locus_tags",
    value: "length_affects_onset",
    label: "Len. → Onset",
    description: "Length affects onset",
  },
  {
    key: "locus_tags",
    value: "length_affects_phenotype",
    label: "Len. → Pheno.",
    description: "Length affects phenotype",
  },
  {
    key: "locus_tags",
    value: "length_affects_penetrance",
    label: "Len. → Pen.",
    description: "Length affects penetrance",
  },
  {
    key: "locus_tags",
    value: "motif_affects_severity",
    label: "Mot. → Sev.",
    description: "Motif affects severity",
  },
  {
    key: "locus_tags",
    value: "motif_affects_onset",
    label: "Mot. → Onset",
    description: "Motif affects onset",
  },
  {
    key: "locus_tags",
    value: "motif_affects_phenotype",
    label: "Mot. → Pheno.",
    description: "Motif affects phenotype",
  },
  {
    key: "locus_tags",
    value: "motif_affects_penetrance",
    label: "Mot. → Pen.",
    description: "Motif affects penetrance",
  },
  {
    key: "locus_tags",
    value: "motif_affects_instability",
    label: "Mot. → Inst.",
    description: "Motif affects instability",
  },
  {
    key: "locus_tags",
    value: "proposed_modifier",
    label: "Prop. Mod.",
    description: "Proposed modifier",
  },

  /** disease tags */
  {
    key: "disease_tags",
    value: "oculopharyngodistal_myopathy",
    label: "OPDM",
    description: "Oculopharyngodistal myopathy",
  },
  {
    key: "disease_tags",
    value: "phenotypic_spectrum",
    label: "Pheno. Spec.",
    description: "Phenotypic spectrum",
  },
  {
    key: "disease_tags",
    value: "ataxia",
    label: "Ataxia",
    description: "Ataxia",
  },
  {
    key: "disease_tags",
    value: "epilepsy",
    label: "Epilepsy",
    description: "Epilepsy",
  },
  {
    key: "disease_tags",
    value: "spinocerebellar_ataxia",
    label: "SCA",
    description: "Spinocerebellar ataxia",
  },
  {
    key: "disease_tags",
    value: "myotonic_dystrophy",
    label: "Myo. Dys.",
    description: "Myotonic dystrophy",
  },
  {
    key: "disease_tags",
    value: "hand_foot_genital_syndrome",
    label: "HFGS",
    description: "Hand-foot-genital syndrome",
  },
  {
    key: "disease_tags",
    value: "familial_adult_myoclonic_epilepsy",
    label: "FAME",
    description: "Familial adult myoclonic epilepsy",
  },
  {
    key: "disease_tags",
    value: "phenotypic_spectrum",
    label: "Pheno. Spec.",
    description: "Phenotypic spectrum",
  },

  /** evidence */
  {
    key: "evidence",
    ...classifications.Definitive,
    Icon: IconCheck,
    className: "bg-best text-white",
    filter: true,
  },
  {
    key: "evidence",
    ...classifications.Strong,
    Icon: IconPlus,
    className: "bg-good/50 text-black",
    filter: true,
  },
  {
    key: "evidence",
    ...classifications.Moderate,
    Icon: IconPoint,
    className: "bg-okay/25 text-black",
    filter: true,
  },
  {
    key: "evidence",
    ...classifications.Limited,
    Icon: IconMinus,
    className: "bg-warn/25 text-black",
    filter: true,
  },
  {
    key: "evidence",
    ...classifications.Disputed,
    Icon: IconExclamationMark,
    className: "bg-bad/50 text-black",
    filter: true,
  },
  {
    key: "evidence",
    ...classifications.Refuted,
    Icon: IconX,
    className: "bg-worst text-white",
    filter: true,
  },
  {
    key: "evidence",
    ...classifications["No Known Relationship"],
    Icon: IconSlash,
    className: "bg-warn/25 text-black",
    filter: true,
  },
  {
    key: "evidence",
    ...classifications.Provisional,
    Icon: IconQuestionMark,
    className: "bg-gray text-white",
    filter: true,
  },
];
