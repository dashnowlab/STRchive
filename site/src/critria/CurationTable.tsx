import Button from "@/components/Button";
import Popover from "@/components/Popover";
import Table from "@/components/Table";
import Tag from "@/components/Tag";
import { curations } from "@/data";
import { tagOptions } from "@/data/tags";
import { IconArrowRight } from "@tabler/icons-react";

/** table for main critria page */
export default function CurationTable() {
  const mappedCurations = curations.map((curation) => ({
    ...curation,
    classification_index: tagOptions.findIndex(
      (tag) => curation.classification === tag.value,
    ),
  }));

  type Datum = (typeof mappedCurations)[number];

  return (
    <Table
      itemNames="curations"
      rows={mappedCurations}
      columns={[
        {
          key: "Locus_ID",
          render: (cell: Datum["Locus_ID"]) => (
            <Popover content="Go to curation page" button={false}>
              <Button className="p-0!" design="bubble" to={`/critria/${cell}`}>
                <IconArrowRight />
              </Button>
            </Popover>
          ),
          sortable: false,
        },
        {
          key: "Gene",
          name: "Gene",
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
          name: "Score",
        },
        {
          /** use number value so column sorted by that instead of alphabetically */
          key: "classification_index",
          name: "Classification",
          className: "py-0!",
          render: (cell: Datum["classification_index"], row: Datum) => (
            <Tag value={row.classification} className="w-full" />
          ),
        },
        {
          key: "Date",
          name: "Date",
          render: (cell: Datum["Date"]) =>
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
      ]}
      sort={[{ id: "1", desc: false }]}
    />
  );
}
