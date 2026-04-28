import Link from "@/components/Link";
import Table from "@/components/Table";
import Tag from "@/components/Tag";
import { tagOptions } from "@/data/tags";

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
    render: (cell, row) => <Tag value={row.classification} />,
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
  const mappedCurations = curations.map((curation) => ({
    ...curation,
    classification_index: tagOptions.findIndex(
      (tag) => curation.classification === tag.value,
    ),
  }));

  return (
    <Table
      cols={cols}
      rows={mappedCurations}
      sort={[{ id: "0", desc: false }]}
    />
  );
};

export default CurationTable;
