import { BsStars } from "react-icons/bs";
import { FaCircleExclamation, FaCircleQuestion } from "react-icons/fa6";
import { newThreshold } from "./derived";

/** top-level tag types */
export const tagOptions = [
  /** "important" tags */
  {
    value: "new",
    label: "New",
    tooltip: `Relatively new (less than ~${newThreshold} years old)`,
    Icon: BsStars,
    color: `var(--primary)`,
    important: true,
  },
  {
    value: "conflicting_evidence",
    label: "Conflicting",
    tooltip: "Conflicting evidence",
    Icon: FaCircleExclamation,
    color: "var(--secondary)",
    important: true,
  },
  {
    value: "sparse_evidence",
    label: "Sparse",
    tooltip: "Sparse evidence (< 5 independent cases after 5 years)",
    Icon: FaCircleQuestion,
    color: `var(--tertiary)`,
    important: true,
  },

  /** locus tags */
  {
    value: "somatic_instability",
    label: "Som. Inst.",
    tooltip: "Somatic instability",
  },
  {
    value: "anticipation",
    label: "Anticip.",
    tooltip: "Anticipation",
  },
  {
    value: "contraction",
    label: "Contrac.",
    tooltip: "Contraction",
  },
  {
    value: "maternal_expansion",
    label: "Mat. Exp.",
    tooltip: "Maternal expansion",
  },
  {
    value: "paternal_expansion",
    label: "Pat. Exp.",
    tooltip: "Paternal expansion",
  },
  {
    value: "length_affects_severity",
    label: "Len. → Sev.",
    tooltip: "Length affects severity",
  },
  {
    value: "length_affects_onset",
    label: "Len. → Onset",
    tooltip: "Length affects onset",
  },
  {
    value: "length_affects_phenotype",
    label: "Len. → Pheno.",
    tooltip: "Length affects phenotype",
  },
  {
    value: "length_affects_penetrance",
    label: "Len. → Pen.",
    tooltip: "Length affects penetrance",
  },
  {
    value: "motif_affects_severity",
    label: "Mot. → Sev.",
    tooltip: "Motif affects severity",
  },
  {
    value: "motif_affects_onset",
    label: "Mot. → Onset",
    tooltip: "Motif affects onset",
  },
  {
    value: "motif_affects_phenotype",
    label: "Mot. → Pheno.",
    tooltip: "Motif affects phenotype",
  },
  {
    value: "motif_affects_penetrance",
    label: "Mot. → Pen.",
    tooltip: "Motif affects penetrance",
  },
  {
    value: "motif_affects_instability",
    label: "Mot. → Inst.",
    tooltip: "Motif affects instability",
  },
  {
    value: "proposed_modifier",
    label: "Prop. Mod.",
    tooltip: "Proposed modifier",
  },

  /** disease tags */
  {
    value: "oculopharyngodistal_myopathy",
    label: "OPDM",
    tooltip: "Oculopharyngodistal myopathy",
  },
  {
    value: "phenotypic_spectrum",
    label: "Pheno. Spec.",
    tooltip: "Phenotypic spectrum",
  },
  {
    value: "ataxia",
    label: "Ataxia",
    tooltip: "Ataxia",
  },
  {
    value: "epilepsy",
    label: "Epilepsy",
    tooltip: "Epilepsy",
  },
  {
    value: "spinocerebellar_ataxia",
    label: "SCA",
    tooltip: "Spinocerebellar ataxia",
  },
  {
    value: "myotonic_dystrophy",
    label: "Myo. Dys.",
    tooltip: "Myotonic dystrophy",
  },
  {
    value: "hand_foot_genital_syndrome",
    label: "HFGS",
    tooltip: "Hand-foot-genital syndrome",
  },
  {
    value: "familial_adult_myoclonic_epilepsy",
    label: "FAME",
    tooltip: "Familial adult myoclonic epilepsy",
  },
  {
    value: "phenotypic_spectrum",
    label: "Pheno. Spec.",
    tooltip: "Phenotypic spectrum",
  },
];
