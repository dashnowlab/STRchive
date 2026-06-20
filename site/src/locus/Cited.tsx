import type { Citation } from "@/data";
import { Fragment } from "react";
import Link from "@/components/Link";
import Popover from "@/components/Popover";

type Props = {
  value?:
    | {
        text: string;
        references: (Citation & { number: number })[];
      }[]
    | null;
};

/** render field value that has in-text references */
export default function Cited({ value }: Props) {
  if (!value) return "–";
  return value.map(({ text, references }, index) =>
    text ? (
      /** plain text */
      <span key={index}>{text}</span>
    ) : (
      <sup key={index}>
        {references.map(({ id, number, title, authors, publisher }, index) => (
          <Fragment key={index}>
            {/* link to citation id below */}
            <Popover
              content={
                <>
                  {title && <div className="line-clamp-2">{title}</div>}
                  {authors && (
                    <div className="line-clamp-1">{authors.join(" ")}</div>
                  )}
                  {publisher && <div className="line-clamp-1">{publisher}</div>}
                  {id}
                </>
              }
              button={false}
            >
              <Link to={`#${id}`}>{number}</Link>
            </Popover>
            {index < references.length - 1 && ","}
          </Fragment>
        ))}{" "}
      </sup>
    ),
  );
}
