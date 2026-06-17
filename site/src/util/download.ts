type Filename = string | string[];

/** download url as file */
const download = (
  /** url to download */
  url: string,
  /** single filename string or filename "parts" */
  filename: Filename,
  /** extension, with or without dot */
  ext: string,
) => {
  const link = document.createElement("a");
  link.href = url;
  link.download =
    [import.meta.env.PUBLIC_TITLE, filename]
      .flat()
      .filter(Boolean)
      .map((part) =>
        part
          /** make path safe */
          .replace(/[^A-Za-z0-9]+/g, "-")
          /** remove leading/trailing dashes */
          .replace(/(^-+)|(-+$)/g, ""),
      )
      .filter(Boolean)
      .join("_")
      /** remove extension if already included */
      .replace(new RegExp("." + ext + "$"), "") +
    "." +
    ext;
  link.click();
  window.URL.revokeObjectURL(url);
};

/** make url from blob */
export const getUrl = (
  /** blob data to download */
  data: Blob | string,
  /** mime type */
  type: string,
) =>
  typeof data === "string" && data.startsWith("data:")
    ? data
    : window.URL.createObjectURL(new Blob([data], { type }));

/** download data as json */
export const downloadJson = (data: unknown, filename: Filename) =>
  download(getUrl(JSON.stringify(data), "application/json"), filename, "json");
