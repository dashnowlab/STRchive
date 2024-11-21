/** wait */
export const sleep = (ms = 0) =>
  new Promise((resolve) => window.setTimeout(resolve, ms));
