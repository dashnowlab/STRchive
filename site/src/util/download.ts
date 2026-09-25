import papaparse from "papaparse";

const getFilename = (filename: string) =>
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
    .join("_");

/** download url as file */
const downloadFile = (
  /** url to download */
  url: string,
  /** single filename string */
  filename: string,
  /** extension, with dot */
  ext: string,
) => {
  let download = getFilename(filename);

  /** add extension */
  if (!download.endsWith(ext)) download += ext;

  /** trigger download */
  const link = document.createElement("a");
  link.href = url;
  link.download = download;
  link.click();
  window.URL.revokeObjectURL(url);
};

/** make url from data */
const getUrl = (
  /** data to download */
  data: string | BlobPart | Blob,
  /** mime type */
  type?: string,
) =>
  typeof data === "string" && data.startsWith("data:")
    ? data
    : window.URL.createObjectURL(
        data instanceof Blob ? data : new Blob([data], { type }),
      );

/** download data as json file */
export const downloadJson = (data: unknown, filename: string) =>
  downloadFile(
    getUrl(JSON.stringify(data, null, 2), "application/json;charset=utf-8"),
    filename,
    ".json",
  );

export type Tabular = (string | number | boolean | null | undefined)[][];

/** assemble csv/tsv from arrays */
const stringifyTable = (table: Tabular, delimiter = "\t") =>
  papaparse.BYTE_ORDER_MARK + papaparse.unparse(table, { delimiter });

/** download data as csv file */
export const downloadCsv = (data: Tabular, filename: string) =>
  downloadFile(
    getUrl(stringifyTable(data, ","), "text/csv;charset=utf-8"),
    filename,
    ".csv",
  );

/** download data as tsv file */
export const downloadTsv = (data: Tabular, filename: string) =>
  downloadFile(
    getUrl(
      stringifyTable(data, "\t"),
      "text/tab-separated-values;charset=utf-8",
    ),
    filename,
    ".tsv",
  );
