import type curations from "~/criTRia-curations.json";
import Button from "@/components/Button";
import Table, { defineData } from "@/components/Table";
import Tag from "@/components/Tag";
import { tagOptions } from "@/data/tags";
import { IconArrowRight } from "@tabler/icons-react";

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
          key: "Locus_ID",
          render: (cell) => (
            <Button
              to={`/critria/${cell}`}
              className="p-0!"
              design="bubble"
              data-tooltip="Go to curation page"
            >
              <IconArrowRight />
            </Button>
          ),
          sortable: false,
        }),
        column({
          key: "Gene",
          name: "Gene",
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
        }),
        column({
          /** use number value so column sorted by that instead of alphabetically */
          key: "classification_index",
          name: "Classification",
          className: "py-0!",
          render: (cell, row) => (
            <Tag value={row.classification} className="w-full" />
          ),
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
      sort={[{ id: "1", desc: false }]}
    />
  );
}
