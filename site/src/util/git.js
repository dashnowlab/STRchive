import { execSync } from "child_process";
import { readFileSync } from "fs";
import { maxBy } from "lodash-es";

/** get git blame info of json file entry. makes assumptions about structure and formatting. */
export const getJsonBlame = (file, regex) => {
  /** split file into lines */
  const lines = readFileSync(file, "utf-8").split("\n");

  /** essentially make line numbers 1-indexed */
  lines.unshift("");

  /** find line number of id */
  const line = lines.findIndex((line) => line.match(new RegExp(regex)));

  if (!line) return {};

  /** find start/end lines of entry */
  let start = line;
  while (!lines[start].startsWith("{")) start--;
  let end = line;
  while (!lines[end].startsWith("}")) end++;

  /** get git blame between start/end lines */
  const blame = execSync(`git blame -L ${start},${end} ${file}`)
    .toString()
    .split("\n")
    .map((line) => {
      let [, hash, author, date] =
        /** extract info */
        line.match(/([A-Za-z0-9]+) \((.+) (\d\d\d\d-\d\d-\d\d)/i) ?? [];
      return { hash, author, date: new Date(date) };
    });

  /** get most recent blame line */
  let { hash, author, date } = maxBy(blame, "date");

  /** format date */
  date = new Date(date).toLocaleString(undefined, { dateStyle: "medium" });

  /** get project version */
  const version = getProjectVersion(hash);

  return { start, end, line, hash, author, date, version };
};

/** get version string from citation file */
const getProjectVersion = (hash) => {
  const [, version] =
    getFileVersion("../CITATION.cff", hash).match(/version: (.+)/i) ?? [];
  return version;
};

/** get specific version of file */
const getFileVersion = (file, hash) =>
  execSync(hash ? `git show ${hash}:${file}` : `cat ${file}`).toString();
