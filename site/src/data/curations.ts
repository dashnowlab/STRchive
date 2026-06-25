export const classifications = {
  Definitive: {
    value: "Definitive",
    label: "Definitive",
    score: "12 to 18",
    description:
      "There is definitive evidence for this locus-disease relationship",
  },
  Strong: {
    value: "Strong",
    label: "Strong",
    score: "12 to 18",
    description: "There is strong evidence for this locus-disease relationship",
  },
  Moderate: {
    value: "Moderate",
    label: "Moderate",
    score: "7 to 11",
    description:
      "There is moderate evidence for this locus-disease relationship",
  },
  Limited: {
    value: "Limited",
    label: "Limited",
    score: "0.1 to 6",
    description:
      "There is limited evidence for this locus-disease relationship",
  },
  Disputed: {
    value: "Disputed",
    label: "Disputed",
    score: "any",
    description: "Evidence disputes this locus-disease relationship",
  },
  Refuted: {
    value: "Refuted",
    label: "Refuted",
    score: "any",
    description: "Evidence refutes this locus-disease relationship",
  },
  "No Known Relationship": {
    value: "No Known Relationship",
    label: "No Known",
    score: "0",
    description: "No known evidence for this locus-disease relationship",
  },
  Provisional: {
    value: "Provisional",
    label: "Provisional",
    score: "",
    description:
      "This locus-disease relationship has been proposed but has not yet been evaluated",
  },
} as const;
