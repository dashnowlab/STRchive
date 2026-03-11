/** order and color hues of each classification type */
export const classifications = [
  {
    match: /no known/i,
    backgroundColor: "var(--gray)",
    color: "var(--black)",
  },
  {
    match: /refuted/i,
    backgroundColor: "#8C2D2D",
    color: "var(--white)",
  },
  {
    match: /disputed/i,
    backgroundColor: "#D66A4E",
    color: "var(--white)",
  },
  {
    match: /limited/i,
    backgroundColor: "#E9B872",
    color: "var(--black)",
  },
  {
    match: /moderate/i,
    backgroundColor: "#8FBBD9",
    color: "var(--black)",
  },
  {
    match: /strong/i,
    backgroundColor: "#2A9D8F",
    color: "var(--white)",
  },
  {
    match: /definitive/i,
    backgroundColor: "#0B6E6E",
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
