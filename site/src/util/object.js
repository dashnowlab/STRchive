/** get values of nested object
 * e.g. getValues({ a: 1, b: [2, { b: 3 } ] }) => [1, 2, 3] */
export const getValues = (object) =>
  object && typeof object === "object"
    ? Object.values(object).map(getValues).flat()
    : [object];
