import { Fragment } from "react";
import { uniq } from "lodash-es";
import Link from "@/components/Link";

/** render free text with citations */
const Cited = ({ text, citations }) => {
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
    .map((_, index) =>
      text.substring(indices[index], indices[index + 1]).trim(),
    )
    .map((text) => {
      /** if citation */
      if (text.startsWith("["))
        return {
          citations: text
            /** widdle down to just citation id */
            .replaceAll("@", "")
            .replaceAll("[", "")
            .replaceAll("]", "")
            /** split multiple */
            .split(";")
            .map((citation) => citation.trim())
            .map(
              (citation) =>
                /** look up full citation details  */
                citations.find(
                  ({ id }) => id === citation,
                ) /** if citation doesn't exist, fall-back to plain text */ ?? {
                  text: citation,
                },
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
        {citations?.map(
          ({ text, number, title, authors, publisher }, index) => (
            <Fragment key={index}>
              {text ? (
                /** plain text fallback */
                <span data-tooltip={text}>#</span>
              ) : (
                /** link to citation id below */
                <Link
                  to={`#citation-${number}`}
                  data-tooltip={[title, authors?.join(", "), publisher]
                    .filter(Boolean)
                    .map((line) => `<div class="truncate-lines">${line}</div>`)
                    .join("")}
                >
                  {number}
                </Link>
              )}
              {index < citations.length - 1 && ","}
            </Fragment>
          ),
        )}{" "}
      </sup>
    ),
  );
};

export default Cited;
