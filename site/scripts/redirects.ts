import { writeFileSync } from "fs";
import { groupBy, isEqual, max, orderBy, uniq, uniqWith } from "lodash-es";
import { simpleGit } from "simple-git";
import loci from "~/STRchive-loci.json";

/** pairs of redirects */
const redirects: Record<string, string> = {};

/** add redirect */
const add = (from: string, to: string) => (redirects[from] = to);

/** legacy pages */
add("/index.html", "/");
add("/about.html", "/about");
add("/contribute.html", "/about");
add("/resources.html", "/resources");
add("/plots.html", "/loci#plots");
add("/table.html", "/loci");
add("/database/index.html", "/loci");

/** legacy gene pages */
for (const [gene, ids] of Object.entries(groupBy(loci, ({ gene }) => gene))) {
  const from = `/database/${gene}.html`;
  const first = ids[0]?.id;
  const to = first ? `/loci/${first}` : "/loci";
  add(from, to);
}

/** direct gene links */
for (const [gene, ids] of Object.entries(groupBy(loci, ({ gene }) => gene))) {
  const from = `/loci/${gene}`;
  const first = ids[0]?.id;
  const to = first ? `/loci/${first}` : "/loci";
  add(from, to);
}

/** common misspellings */
add("/criteria", "/critria");
add("/criteria/*", "/critria/:splat");
add("/criTRia", "/critria");
add("/criTRia/*", "/critria/:splat");
add("/locus", "/loci");
add("/locus/*", "/loci/:splat");

/** file paths, relative to /site */
const lociFile = "../data/STRchive-loci.json";
const curationsFile = "../data/criTRia-curations.json";
const redirectsFile = "./public/_redirects";

type Locus = { id: string } & Record<string, unknown>;
type Curation = { Locus_ID: string } & Record<string, unknown>;

/** git client */
const git = simpleGit();

/** get latest main with full history */
await git.fetch(["--no-tags", "--depth=2147483647", "origin", "main"]);

/** get pairs of entries */
const pairs = <Entry>(array: Entry[]): [Entry, Entry][] =>
  array.slice(1).map((after, index) => [array[index], after]);

/** fraction of values that are equal between two entries */
const compare = <Entry extends Record<string, unknown>>(a: Entry, b: Entry) => {
  const keys = uniq(Object.keys(a).concat(Object.keys(b)));
  return keys.filter((key) => isEqual(a[key], b[key])).length / keys.length;
};

/** get full list of commit hashes that touched file */
const getHashes = async (file: string) => {
  /** commit hashes on main that touched loci file, oldest first */
  const hashes = (await git.log(["origin/main", "--", file])).all
    .map(({ hash }) => `${hash}:${file}`)
    .reverse();
  /** add merge base hash (file as it will look after merging this branch into main) */
  hashes.push(
    `${(await git.raw("merge-base", "origin/main", "HEAD")).trim()}:${file}`,
  );
  /** add final, current hash (file as it looks right now on current branch) */
  hashes.push(`${(await git.revparse(["HEAD"])).trim()}:${file}`);
  return hashes;
};

/** get file contents as of each commit, oldest first (skips commits with invalid json) */
const getSnapshots = async <Type>(hashes: string[]) =>
  (
    await Promise.all(
      hashes.map(async (hash) => {
        try {
          const contents = await git.show(hash);
          return JSON.parse(contents) as Type[];
        } catch (error) {
          console.warn(`Skipping invalid JSON for ${hash}`, error);
          return undefined;
        }
      }),
    )
  ).filter((snapshot) => !!snapshot);

/** get renames of primary id field on list of entries */
const getRenames = <
  Key extends string,
  Entry extends Record<string, unknown> & { [key in Key]: string },
>(
  snapshots: Entry[][],
  key: Key,
) => {
  const renames: Record<string, string> = {};
  for (const [before, after] of pairs(snapshots)) {
    /** file entries whose ids disappeared/appeared */
    const beforeIds = new Set(before.map((entry) => entry[key]));
    const afterIds = new Set(after.map((entry) => entry[key]));
    const removedEntries = before.filter((entry) => !afterIds.has(entry[key]));
    const addedEntries = after.filter((entry) => !beforeIds.has(entry[key]));

    for (const removedEntry of removedEntries) {
      /** find most closely matching entry between removed and added */
      const { similarity, entry: addedEntry } =
        addedEntries
          .map((entry) => ({ similarity: compare(removedEntry, entry), entry }))
          .sort((a, b) => b.similarity - a.similarity)[0] ?? {};

      /** consider before/after entries to represent same locus if similarity high enough */
      if (addedEntry && similarity > 0.25)
        renames[removedEntry[key]] = addedEntry[key];
    }
  }
  return Object.entries(renames);
};

/** loci id renames */
const lociRenames = getRenames(
  await getSnapshots<Locus>(await getHashes(lociFile)),
  "id",
);
/** curation id renames */
const curationRenames = getRenames(
  await getSnapshots<Curation>(await getHashes(curationsFile)),
  "Locus_ID",
);

/** unique id renames */
const renames = orderBy(
  uniqWith(lociRenames.concat(curationRenames), isEqual),
  ([from]) => from,
);

/** add id rename redirects */
for (const [from, to] of renames) add(`/loci/${from}`, `/loci/${to}`);
for (const [from, to] of renames) add(`/loci/${from}/*`, `/loci/${to}/:splat`);
for (const [from, to] of renames) add(`/critria/${from}`, `/critria/${to}`);
for (const [from, to] of renames)
  add(`/critria/${from}/*`, `/critria/${to}/:splat`);

/** for each start of id rename chain */
for (let [from, to] of Object.entries(redirects)) {
  /** follow redirects to end of chain (latest id) */
  const maxTries = 20;
  /** limit tries to catch possible infinite loops */
  for (let tries = 0; redirects[to]; tries++) {
    if (tries >= maxTries)
      throw Error("Possible infinite loop in rename chain");
    to = redirects[to];
  }
  redirects[from] = to;
}

/** get longest strings widths */
const longestFrom = max(Object.keys(redirects).map((from) => from.length)) || 0;
const longestTo = max(Object.values(redirects).map((to) => to.length)) || 0;

/** redirects text to write */
const lines = Object.entries(redirects).map(([from, to]) =>
  [
    globalThis.encodeURI(from).padEnd(longestFrom),
    globalThis.encodeURI(to).padEnd(longestTo),
    "301",
  ].join("  "),
);
/** add warning */
lines.unshift("# AUTO-GENERATED, DO NOT EDIT");
/** add EOF newline */
lines.push("");

/** save redirects */
writeFileSync(redirectsFile, lines.join("\n"));
