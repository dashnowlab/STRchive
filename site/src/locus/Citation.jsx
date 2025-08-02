import clsx from "clsx";
import Link from "@/components/Link";
import ShowMoreLines from "@/components/ShowMoreLines";
import classes from "./Citation.module.css";

/** full citation details for a source */
const Citation = ({ number, id, title, authors, publisher, date, link }) => {
  const details = [publisher, date].filter(Boolean);
  return (
    <div className={clsx("col", classes.citation)} id={id}>
      {number && <div className={clsx("row", classes.number)}>{number}</div>}
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
