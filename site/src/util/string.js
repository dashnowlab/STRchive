/** shorten url text */
export const shortenUrl = (value) => {
  try {
    const url = new URL(value);
    return url.hostname + url.pathname;
  } catch (error) {
    return value;
  }
};

/** make string url-safe */
export const slugify = (value) =>
  value
    .toLowerCase()
    .replaceAll(/[^a-z0-9]+/g, " ")
    .replaceAll(/\s+/g, " ")
    .trim()
    .replaceAll(" ", "-");
