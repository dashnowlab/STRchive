import { uniq } from "lodash-es";

/** render citations in free text */
const Cited = ({ text }) => {
  const matches = uniq([
    0,
    ...[...text.matchAll(/\[@\w+:.+?\]/g)]
      .map((match) => [match.index, match.index + match[0].length])
      .flat(),
    text.length,
  ]);

  const split = Array(matches.length - 1)
    .fill(null)
    .map((_, index) => text.substring(matches[index], matches[index + 1]))
    .map((text) => {
      if (text.includes("@"))
        return {
          references: text
            .replaceAll("@", "")
            .replaceAll("[", "")
            .replaceAll("]", "")
            .split(";")
            .map((ref) => ref.trim()),
        };
      else return { text };
    });

  console.log(matches, split);

  return text;
};

export default Cited;
