import { readFileSync, writeFileSync } from "node:fs";

const read = (path) =>
  JSON.parse(readFileSync(new URL(path, import.meta.url), "utf8"));

/**
 * make netlify redirect rules from each locus's previous ids
 * https://docs.netlify.com/manage/routing/redirects/overview
 */
export const makeRedirects = (loci, curations) => {
  const curated = new Set(curations.map(({ Locus_ID }) => Locus_ID));
  const rules = [];
  for (const { id, previous_ids } of loci) {
    for (const old of previous_ids ?? []) {
      const from = encodeURI(old);
      const to = encodeURI(id);
      rules.push([`/loci/${from}`, `/loci/${to}`]);
      rules.push([`/loci/${from}/*`, `/loci/${to}/:splat`]);
      if (curated.has(id)) rules.push([`/critria/${from}`, `/critria/${to}`]);
    }
  }
  return [
    "# Generated from previous_ids in data/STRchive-loci.json. Do not edit.",
    ...rules.map(([from, to]) => `${from} ${to} 301`),
    "",
  ].join("\n");
};

/** astro integration that writes netlify _redirects file to build output */
export const redirects = () => ({
  name: "redirects",
  hooks: {
    "astro:build:done": ({ dir, logger }) => {
      const text = makeRedirects(
        read("../data/STRchive-loci.json"),
        read("../data/criTRia-curations.json"),
      );
      writeFileSync(new URL("_redirects", dir), text);
      logger.info(`wrote ${text.split("\n").length - 2} redirect rules`);
    },
  },
});
