import type { Citation } from "@/data/types";
import Link from "@/components/Link";
import ShowMoreLines from "@/components/ShowMoreLines";
import { IconBookmark } from "@tabler/icons-react";

type Props = Citation & { number?: number };

/** full citation details for a source */
export default function Citation({
  number,
  id,
  title,
  authors,
  publisher,
  date,
  link,
}: Props) {
  const details = [publisher, date].filter(Boolean);
  return (
    <div className="box" id={id}>
      {number && (
        <strong>
          <IconBookmark />
          {number}
        </strong>
      )}
      {title && <ShowMoreLines lines={1}>{title}</ShowMoreLines>}
      {!!authors?.length && (
        <ShowMoreLines lines={1}>{authors?.join(", ")}</ShowMoreLines>
      )}
      {!!details.length && <div>{details.join(" · ")}</div>}
      {id && (
        <Link to={link ?? `https://www.google.com/search?q=${id}`}>{id}</Link>
      )}
    </div>
  );
}
