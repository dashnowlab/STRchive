import Link from "@/components/Link";
import Table from "@/components/Table";
import { classifications, getStyles } from "@/data/critria";
import classes from "./CurationTable.module.css";

/** column definitions */
const cols = [
  {
    key: "Gene",
    name: "Gene",
    render: (cell, row) => (
      <Link to={`/critria/${row.Locus_ID}`} arrow={false}>
        {cell}
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
    /** use classification number value so col sorted by that instead of alphabetically */
    key: "classificationValue",
    name: "Classification",
    style: {
      padding: 0,
    },
    render: (cell, row) => {
      /** display normal classification string */
      const value = row.classification;
      return (
        <div className={classes.classification} style={getStyles(value)}>
          {value}
        </div>
      );
    },
  },
  {
    key: "Date",
    name: "Date",
    render: (cell) =>
      new Date(cell).toLocaleDateString(undefined, {
        year: "numeric",
        month: "short",
        day: "numeric",
      }),
  },
];

/** table for main critria page */
const CurationTable = ({ curations }) => {
  const mappedCurations = curations.map((curation) => ({
    ...curation,
    classificationValue: classifications.findIndex(({ match }) =>
      curation.classification.match(match),
    ),
  }));

  return <Table cols={cols} rows={mappedCurations} />;
};

export default CurationTable;
