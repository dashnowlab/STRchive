import type { Citation } from "@/data/types";
import { Fragment } from "react";
import Link from "@/components/Link";

type Props = {
  value: {
    text: string;
    references: (Citation & { number: number })[];
  }[];
};

/** render field value that has in-text references */
export default function Cited({ value }: Props) {
  return value.map(({ text, references }, index) =>
    text ? (
      /** plain text */
      <span key={index}>{text}</span>
    ) : (
      <sup key={index}>
        {references.map(({ id, number, title, authors, publisher }, index) => (
          <Fragment key={index}>
            ( /** link to citation id below */
            <Link
              to={`#${id}`}
              onClick={(event) => event.stopPropagation()}
              data-tooltip={[title, authors?.join(", "), publisher]
                .filter(Boolean)
                .map(
                  (line) =>
                    `<div class="overflow-hidden [display:-webkit-box] [-webkit-box-orient:vertical] [-webkit-line-clamp:var(--lines,1)]">${line}</div>`,
                )
                .join("")}
            >
              {number}
            </Link>
            ){index < references.length - 1 && ","}
          </Fragment>
        ))}{" "}
      </sup>
    ),
  );
}
