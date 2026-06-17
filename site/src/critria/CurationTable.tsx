import type curations from "~/criTRia-curations.json";
import Link from "@/components/Link";
import Table, { defineData } from "@/components/Table";
import Tag from "@/components/Tag";
import { tagOptions } from "@/data/tags";

type Curations = typeof curations;

type Props = {
  curations: Curations;
};

/** table for main critria page */
export default function CurationTable({ curations }: Props) {
  const mappedCurations = curations.map((curation) => ({
    ...curation,
    classification_index: tagOptions.findIndex(
      (tag) => curation.classification === tag.value,
    ),
  }));

  return (
    <Table
      {...defineData(mappedCurations, (column) => [
        column({
          key: "Gene",
          name: "Gene",
          render: (cell, row) => (
            <Link to={`/critria/${row.Locus_ID}`} arrow={false}>
              {cell}
            </Link>
          ),
        }),
        column({
          key: "Disease_ID",
          name: "Disease",
        }),
        column({
          key: "Inheritance",
          name: "Inheritance",
        }),
        column({
          key: "total_score",
          name: "Total Score",
          style: { textAlign: "center" },
        }),
        column({
          /** use number value so column sorted by that instead of alphabetically */
          key: "classification_index",
          name: "Classification",
          style: { padding: 0 },
          render: (cell, row) => <Tag value={row.classification} />,
        }),
        column({
          key: "Date",
          name: "Date",
          render: (cell) =>
            new Date(cell).toLocaleDateString(undefined, {
              year: "numeric",
              month: "short",
              day: "numeric",
            }),
        }),
        column({
          key: "Source",
          name: "Source",
        }),
      ])}
      sort={[{ id: "Gene", desc: false }]}
    />
  );
}
