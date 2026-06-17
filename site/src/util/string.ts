/** shorten url text */
export const shortenUrl = (value: string) => {
  try {
    const url = new URL(value);
    return (url.hostname + url.pathname).replace(/\/+$/, "");
  } catch (error) {
    return value;
  }
};

/** make string url-safe */
export const slugify = (value: string) =>
  value
    .toLowerCase()
    .replaceAll(/[^a-z0-9]+/g, " ")
    .replaceAll(/\s+/g, " ")
    .trim()
    .replaceAll(" ", "-");
