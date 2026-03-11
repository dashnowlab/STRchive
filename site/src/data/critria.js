/** order and color hues of each classification type */
export const classifications = [
  {
    match: /no known/i,
    backgroundColor: "var(--gray)",
    color: "var(--black)",
  },
  {
    match: /refuted/i,
    backgroundColor: "rgb(155, 44, 44)",
    color: "var(--white)",
  },
  {
    match: /disputed/i,
    backgroundColor: "rgb(229, 62, 62)",
    color: "var(--black)",
  },
  {
    match: /limited/i,
    backgroundColor: "rgb(252, 129, 130)",
    color: "var(--black)",
  },
  {
    match: /moderate/i,
    backgroundColor: "rgb(104, 211, 145)",
    color: "var(--black)",
  },
  {
    match: /strong/i,
    backgroundColor: "rgb(55, 161, 105)",
    color: "var(--black)",
  },
  {
    match: /definitive/i,
    backgroundColor: "rgb(39, 103, 73)",
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
