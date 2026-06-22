import type { Citation } from "@/data";
import Link from "@/components/Link";
import ShowMoreLines from "@/components/ShowMoreLines";

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
    <div id={id} className="relative pl-4">
      {number && (
        <strong className="absolute top-1 right-full grid size-6 place-content-center rounded-full bg-primary/10">
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
