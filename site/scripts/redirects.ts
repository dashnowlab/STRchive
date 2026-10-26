import { appendFileSync, readFileSync } from "fs";
import { isEqual, uniq } from "lodash-es";
import { simpleGit } from "simple-git";
import { parse, stringify } from "smol-toml";

/** file paths, relative to /site */
const lociFile = "../data/STRchive-loci.json";
const redirectsFile = "../netlify.toml";

type Locus = { id: string } & Record<string, unknown>;

/** git client */
const git = simpleGit();

/** commit hashes on main that touched loci file, oldest first */
const hashes = (await git.log(["origin/main", "--", lociFile])).all
  .map(({ hash }) => hash)
  .reverse();

/** get loci list as of particular commit */
const getLoci = async (hash: string) => {
  try {
    const contents = await git.show(`${hash}:${lociFile}`).catch(() => "null");
    return JSON.parse(contents) as Locus[];
  } catch {
    return undefined;
  }
};

/** fraction of values that are equal between two loci */
const compare = (a: Locus, b: Locus) => {
  const keys = uniq(Object.keys(a).concat(Object.keys(b)));
  return keys.filter((key) => isEqual(a[key], b[key])).length / keys.length;
};

/** track loci renames */
const renames: Record<string, string> = {};

/** loci lists as of each commit, oldest first */
const snapshots = await Promise.all(hashes.map(getLoci));

/** add current working tree (PR changes) as final snapshot */
snapshots.push(JSON.parse(readFileSync(lociFile, "utf8")) as Locus[]);
/** approximates what final state will be after squash merging into main */

/** track loci list as of previous snapshot */
let before: Locus[] | undefined;

for (const after of snapshots) {
  if (before && after) {
    /** loci whose ids disappeared/appeared */
    const beforeIds = new Set(before.map(({ id }) => id));
    const afterIds = new Set(after.map(({ id }) => id));
    const removedLoci = before.filter(({ id }) => !afterIds.has(id));
    const addedLoci = after.filter(({ id }) => !beforeIds.has(id));

    for (const removedLocus of removedLoci) {
      /** find most closely matching locus between removed and added */
      const { similarity, locus: addedLocus } =
        addedLoci
          .map((locus) => ({ similarity: compare(removedLocus, locus), locus }))
          .sort((a, b) => b.similarity - a.similarity)[0] ?? {};

      /** consider before/after entries to represent same locus if similarity high enough */
      if (addedLocus && similarity > 0.9)
        renames[removedLocus.id] = addedLocus.id;
    }
  }

  /** track previous snapshot */
  before = after;
}

/** pairs of redirects */
const redirects: [string, string][] = [];

/** for each start of id rename chain */
for (let [from, to] of Object.entries(renames)) {
  /** follow redirects to end of chain (latest id) */
  const maxTries = 20;
  /** limit tries to catch possible infinite loops */
  for (let tries = 0; renames[to]; tries++) {
    if (tries >= maxTries)
      throw Error("Possible infinite loop in rename chain");
    to = renames[to];
  }
  redirects.push([from, to]);
}

/** existing redirects in config */
const existing = parse(readFileSync(redirectsFile, "utf8"));
type Redirect = { from?: string; to?: string };
const existingRedirects = (existing.redirects ?? []) as Redirect[];

/** new redirect entries */
const newRedirects = redirects
  /** netlify redirect pattern */
  .map(([from, to]) => ({
    from: `/loci/${from}`,
    to: `/loci/${to}`,
    status: 301,
  }))
  /** skip any that already exist */
  .filter(
    (redirect) =>
      !existingRedirects.some(
        ({ from, to }) => from === redirect.from && to === redirect.to,
      ),
  );

/** format */
const newRedirectsText = stringify({ redirects: newRedirects })
  .split("\n")
  .map((line) => (line && !line.startsWith("[") ? "  " + line : line))
  .join("\n");

/** append to bottom, leave rest of file untouched */
if (newRedirects.length) appendFileSync(redirectsFile, "\n" + newRedirectsText);
