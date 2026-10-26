import { readFileSync } from "fs";
import { maxBy } from "lodash-es";
import { simpleGit } from "simple-git";

/** git client */
const git = simpleGit();

/** get git blame info of json file entry. makes assumptions about structure and formatting. */
export const getJsonBlame = async (file: string, regex: string) => {
  /** split file into lines */
  const lines = readFileSync(file, "utf-8").split("\n");

  /** essentially make line numbers 1-indexed */
  lines.unshift("");

  /** find line number of id */
  const line = lines.findIndex((line) => line.match(new RegExp(regex)));

  if (line === -1) return {};

  /** find start/end lines of entry */
  let start = line;
  while (!lines[start].startsWith("{")) start--;
  let end = line;
  while (!lines[end].startsWith("}")) end++;

  /** get git blame between start/end lines */
  const blame = (await git.raw(["blame", "-L", `${start},${end}`, file]))
    .split("\n")
    .map((line) => {
      const [, hash, author, date] =
        /** extract info */
        line.match(/([A-Za-z0-9]+) \((.+) (\d\d\d\d-\d\d-\d\d)/i) ?? [];
      return { hash, author, date: new Date(date) };
    });

  /** get most recent blame line */
  const { hash, author, date } = maxBy(blame, "date") ?? {};

  /** format date */
  const formattedDate = new Date(date ?? "").toLocaleString(undefined, {
    dateStyle: "medium",
  });

  /** get project version */
  const version = getVersion(await getFile("../CITATION.cff", hash));

  return { start, end, line, hash, author, date: formattedDate, version };
};

/** extract version field from yaml string */
const getVersion = (contents?: string) =>
  contents?.match(/^version:\s*(.+)$/im)?.[1] ?? "";

/** get contents of file w/ optional hash version */
const getFile = async (file: string, hash?: string) => {
  try {
    if (!hash) throw Error("no hash");
    const contents = await git.show([`${hash}:${file}`]);
    if (!contents.trim()) throw Error("empty");
    return contents;
  } catch {
    return readFileSync(file, "utf-8");
  }
};
