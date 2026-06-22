import Button from "@/components/Button";
import Popover from "@/components/Popover";
import Table, { defineData } from "@/components/Table";
import Tag from "@/components/Tag";
import { curations } from "@/data";
import { tagOptions } from "@/data/tags";
import { downloadJson } from "@/util/download";
import { IconArrowRight, IconDownload } from "@tabler/icons-react";

/** table for main critria page */
export default function CurationTable() {
  const mappedCurations = curations.map((curation) => ({
    ...curation,
    classification_index: tagOptions.findIndex(
      (tag) => curation.classification === tag.value,
    ),
  }));

  return (
    <>
      <Table
        {...defineData(mappedCurations, (column) => [
          column({
            key: "Locus_ID",
            render: (cell) => (
              <Popover content="Go to curation page" button={false}>
                <Button
                  className="p-0!"
                  design="bubble"
                  to={`/critria/${cell}`}
                >
                  <IconArrowRight />
                </Button>
              </Popover>
            ),
            sortable: false,
          }),
          column({
            key: "Gene",
            name: "Gene",
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
            name: "Score",
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
      <Button
        design="plain"
        className="self-center"
        onClick={() => downloadJson(curations, "curations")}
      >
        <IconDownload />
        Download
      </Button>
    </>
  );
}
