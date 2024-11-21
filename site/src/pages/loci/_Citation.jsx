import ShowMoreLines from "@/components/ShowMoreLines";
import classes from "./_Citation.module.css";

/** full citation details for a source */
const Citation = ({ number, id, title, authors, publisher, date }) => (
  <div
    className={classes.citation}
    id={number ? `citation-${number}` : undefined}
  >
    {number && <strong className={classes.number}>{number}</strong>}
    <ShowMoreLines lines={1}>{title ?? "-"}</ShowMoreLines>
    <ShowMoreLines lines={1}>{authors?.join(", ") || "-"}</ShowMoreLines>
    <div>{[publisher, date].filter(Boolean).join(" · ")}</div>
    <div>{id}</div>
  </div>
);

export default Citation;
