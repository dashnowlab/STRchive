import type { Citation } from "@/data/types";
import Link from "@/components/Link";
import ShowMoreLines from "@/components/ShowMoreLines";

type Props = Citation & { number: number };

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
    <div
      className="relative flex max-w-full flex-col items-start gap-1 px-4 py-2"
      id={id}
    >
      {number && (
        <div className="absolute top-2 left-2 flex size-6 max-w-full flex-wrap items-center justify-center gap-x-4 gap-y-2 rounded-full bg-primary/25">
          {number}
        </div>
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
