import Link from "@/components/Link";
import Table from "@/components/Table";
import { evidenceOptions } from "@/data/evidence";
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
    /** use number value so col sorted by that instead of alphabetically */
    key: "classification_index",
    name: "Classification",
    style: {
      padding: 0,
    },
    render: (cell, row) => {
      /** display normal classification string */
      const value = row.classification;
      return (
        <div
          className={classes.classification}
          style={{ backgroundColor: row.bg, color: row.text }}
        >
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
  {
    key: "Source",
    name: "Source",
  },
];

/** table for main critria page */
const CurationTable = ({ curations }) => {
  const mappedCurations = curations.map((curation) => {
    const index = evidenceOptions.findIndex(
      (tag) => curation.classification === tag.value,
    );
    const tag = evidenceOptions[index];
    return {
      ...curation,
      bg: tag?.bg,
      text: tag?.text,
      classification_index: index,
    };
  });

  return (
    <Table
      cols={cols}
      rows={mappedCurations}
      sort={[{ id: "0", desc: false }]}
    />
  );
};

export default CurationTable;
