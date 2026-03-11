/** order and color hues of each classification type */
export const classifications = [
  {
    match: /no known/i,
    backgroundColor: "var(--gray)",
    color: "var(--black)",
  },
  {
    match: /refuted/i,
    backgroundColor: "color-mix(in srgb, var(--secondary), var(--white) 0%)",
    color: "var(--white)",
  },
  {
    match: /disputed/i,
    backgroundColor: "color-mix(in srgb, var(--secondary), var(--white) 50%)",
    color: "var(--black)",
  },
  {
    match: /limited/i,
    backgroundColor: "color-mix(in srgb, var(--secondary), var(--white) 75%)",
    color: "var(--black)",
  },
  {
    match: /moderate/i,
    backgroundColor: "color-mix(in srgb, var(--primary), var(--white) 75%)",
    color: "var(--black)",
  },
  {
    match: /strong/i,
    backgroundColor: "color-mix(in srgb, var(--primary), var(--white) 50%)",
    color: "var(--black)",
  },
  {
    match: /definitive/i,
    backgroundColor: "color-mix(in srgb, var(--primary), var(--white) 0%)",
    color: "var(--white)",
  },
];

/** get coloring styles for classification */
export const getStyles = (classification) => {
  const match = classifications.find(({ match }) =>
    classification.match(match),
  );
  return match ?? {};
};
