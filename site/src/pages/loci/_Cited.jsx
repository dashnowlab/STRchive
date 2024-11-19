import { Fragment } from "react";
import { uniq } from "lodash-es";
import Link from "@/components/Link";

/** render free text with citations */
const Cited = ({ text, literature }) => {
  if (!text) return "";

  /** find sub-string indices in text that match citation format */
  const indices = uniq([
    0,
    ...[...text.matchAll(/\[@\w+:.+?\]/g)]
      .map((match) => [match.index, match.index + match[0].length])
      .flat(),
    text.length,
  ]);

  const split = Array(indices.length - 1)
    .fill(null)
    /** extract sub-strings from pairs of indices */
    .map((_, index) => text.substring(indices[index], indices[index + 1]))
    .map((text) => {
      /** if citation */
      if (text.includes("@"))
        return {
          citations: text
            /** widdle down to just citation id */
            .replaceAll("@", "")
            .replaceAll("[", "")
            .replaceAll("]", "")
            /** split multiple */
            .split(";")
            .map(
              /** look up full citation details from list of literature */
              (ref) => literature.find((lit) => lit.id === ref.trim()) ?? {},
            ),
        };
      /** else, return plain text */ else return { text };
    });

  return split.map(({ text, citations }, index) =>
    text ? (
      /** plain text */
      <span key={index}>{text}</span>
    ) : (
      <sup key={index}>
        {citations.map(({ number, title }, index) => (
          <Fragment key={index}>
            {/* link to citation id */}
            <Link to={`#citation-${number}`} data-tooltip={`${title ?? "-"}`}>
              {number}
            </Link>
            {index < citations.length - 1 && ","}
          </Fragment>
        ))}
      </sup>
    ),
  );
};

export default Cited;
