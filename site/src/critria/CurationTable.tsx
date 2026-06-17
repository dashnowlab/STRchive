import type curations from "~/criTRia-curations.json";
import Link from "@/components/Link";
import Table, { defineData } from "@/components/Table";
import Tag from "@/components/Tag";
import { tagOptions } from "@/data/tags";

type Props = {
  curations: typeof curations;
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
      {...defineData(mappedCurations, (col) => [
        col({
          key: "Gene",
          name: "Gene",
          render: (cell, row) => (
            <Link to={`/critria/${row.Locus_ID}`} arrow={false}>
              {cell}
            </Link>
          ),
        }),
        col({
          key: "Disease_ID",
          name: "Disease",
        }),
        col({
          key: "Inheritance",
          name: "Inheritance",
        }),
        col({
          key: "total_score",
          name: "Total Score",
          style: { textAlign: "center" },
        }),
        col({
          /** use number value so col sorted by that instead of alphabetically */
          key: "classification_index",
          name: "Classification",
          style: { padding: 0 },
          render: (cell, row) => <Tag value={row.classification} />,
        }),
        col({
          key: "Date",
          name: "Date",
          render: (cell) =>
            new Date(cell).toLocaleDateString(undefined, {
              year: "numeric",
              month: "short",
              day: "numeric",
            }),
        }),
        col({
          key: "Source",
          name: "Source",
        }),
      ])}
      sort={[{ id: "Gene", desc: false }]}
    />
  );
}
