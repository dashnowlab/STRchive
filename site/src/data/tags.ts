import { BsStars } from "react-icons/bs";
import {
  FaCheck,
  FaExclamation,
  FaMinus,
  FaPlus,
  FaPlusMinus,
  FaQuestion,
  FaSlash,
  FaXmark,
} from "react-icons/fa6";
import { newThreshold } from "./derived";

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
    Icon: BsStars,
    /** classes on icon */
    className: "bg-primary text-white",
    /** whether to show in loci page table */
    filter: true,
  },
  {
    key: "gene_strand",
    value: "+",
    Icon: FaPlus,
    className: "bg-primary text-white",
  },
  {
    key: "gene_strand",
    value: "-",
    Icon: FaMinus,
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
    value: "Definitive",
    label: "Definitive",
    Icon: FaCheck,
    description:
      "There is definitive evidence for this locus-disease relationship",
    className: "bg-best text-white",
    filter: true,
  },
  {
    key: "evidence",
    value: "Strong",
    label: "Strong",
    Icon: FaPlus,
    description: "There is strong evidence for this locus-disease relationship",
    className: "bg-good/66 text-black",
    filter: true,
  },
  {
    key: "evidence",
    value: "Moderate",
    label: "Moderate",
    Icon: FaPlusMinus,
    description:
      "There is moderate evidence for this locus-disease relationship",
    className: "bg-okay/33 text-black",
    filter: true,
  },
  {
    key: "evidence",
    value: "Limited",
    label: "Limited",
    Icon: FaMinus,
    description:
      "There is limited evidence for this locus-disease relationship",
    className: "bg-warn/33 text-black",
    filter: true,
  },
  {
    key: "evidence",
    value: "Disputed",
    label: "Disputed",
    Icon: FaExclamation,
    description: "Evidence disputes this gene-disease relationship",
    className: "bg-bad/66 text-black",
    filter: true,
  },
  {
    key: "evidence",
    value: "Refuted",
    label: "Refuted",
    Icon: FaXmark,
    description: "Evidence refutes this gene-disease relationship",
    className: "bg-worst text-white",
    filter: true,
  },
  {
    key: "evidence",
    value: "No Known Relationship",
    label: "No Known",
    Icon: FaSlash,
    description: "No known evidence for this locus-disease relationship",
    className: "bg-warn/33 text-black",
    filter: true,
  },
  {
    key: "evidence",
    value: "Provisional",
    label: "Provisional",
    Icon: FaQuestion,
    description:
      "Provisional: this locus-disease relationship has been proposed but has not yet been evaluated",
    className: "bg-gray text-white",
    filter: true,
  },
];
