import Link from "@/components/Link";
import TableComponent from "@/components/Table";
import classes from "./Table.module.css";

/** order and color hues of each classification type */
const classifications = [
  [/no known/i, 300],
  [/refuted/i, 30],
  [/disputed/i, 80],
  [/limited/i, 100],
  [/moderate/i, 120],
  [/strong/i, 150],
  [/definitive/i, 180],
];

/** hue to full color. use oklch instead of hsl for more even brightness. */
const color = (hue) => `oklch(75% 0.5 ${hue} / 0.25)`;

/** column definitions */
const cols = [
  {
    key: "Gene",
    name: "Gene",
    render: (value, row) => (
      <Link to={`/loci/${row.Locus_ID}`} arrow={false}>
        {value}
      </Link>
    ),
  },
  {
    key: "Disease_ID",
    name: "Disease",
  },
  {
    key: "Inheritance",
    name: "Inheritance",
  },
  {
    key: "total_score",
    name: "Total Score",
    align: "center",
  },
  {
    /** sort by classificationValue */
    key: "classificationValue",
    name: "Classification",
    style: {
      padding: 0,
    },
    render: (_, row) => {
      /** display normal classification string */
      const value = row.classification;
      const hue =
        classifications.find(([regex]) => value.match(regex))?.[1] ?? 0;
      return (
        <div
          className={classes.classification}
          style={{ backgroundColor: color(hue) }}
        >
          {value}
        </div>
      );
    },
  },
  {
    key: "Date",
    name: "Date",
    render: (value) =>
      new Date(value).toLocaleDateString(undefined, {
        year: "numeric",
        month: "short",
        day: "numeric",
      }),
  },
];

/** table for main critria page */
const Table = ({ curations }) => {
  const mappedCurations = curations.map((curation) => ({
    ...curation,
    classificationValue: classifications.findIndex(([regex]) =>
      curation.classification.match(regex),
    ),
  }));

  return <TableComponent cols={cols} rows={mappedCurations} />;
};

export default Table;
