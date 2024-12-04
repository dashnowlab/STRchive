import { Fragment } from "react";
import Link from "@/components/Link";

/** render field value that has in-text references */
const Cited = ({ value }) => {
  return value.map(({ text, references }, index) =>
    text ? (
      /** plain text */
      <span key={index}>{text}</span>
    ) : (
      <sup key={index}>
        {references.map(
          ({ text, id, number, title, authors, publisher }, index) => (
            <Fragment key={index}>
              {text ? (
                /** plain text fallback */
                <span data-tooltip={text}>?</span>
              ) : (
                /** link to citation id below */
                <Link
                  to={`#${id}`}
                  onClick={(event) => event.stopPropagation()}
                  data-tooltip={[title, authors?.join(", "), publisher]
                    .filter(Boolean)
                    .map((line) => `<div class="truncate-lines">${line}</div>`)
                    .join("")}
                >
                  {number}
                </Link>
              )}
              {index < references.length - 1 && ","}
            </Fragment>
          ),
        )}{" "}
      </sup>
    ),
  );
};

export default Cited;
