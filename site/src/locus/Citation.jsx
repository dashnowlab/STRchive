import clsx from "clsx";
import Link from "@/components/Link";
import ShowMoreLines from "@/components/ShowMoreLines";

/** full citation details for a source */
const Citation = ({ number, id, title, authors, publisher, date, link }) => {
  const details = [publisher, date].filter(Boolean);
  return (
    <div
      className="relative flex max-w-full flex-col items-start gap-[5px] px-5 py-2.5  has-[.citation-number]:pl-[45px]"
      id={id}
    >
      {number && (
        <div
          className="citation-number absolute left-2.5 top-2.5 flex h-[25px] w-[25px] max-w-full flex-wrap items-center justify-center gap-x-5 gap-y-2 rounded-full bg-[color-mix(in_srgb,var(--color-primary),var(--color-white)_75%)]"
        >
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
};

export default Citation;
