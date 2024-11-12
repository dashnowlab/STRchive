/** capitalize first char, but unlike lodash's, don't alter rest of string */
export const capitalize = (string) =>
  string.charAt(0).toUpperCase() + string.substring(1);
