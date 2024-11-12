/** get values of nested object */
export const getValues = (object) =>
  object && typeof object === "object"
    ? Object.values(object).map(getValues).flat()
    : [object];
