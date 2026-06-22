import Button from "@/components/Button";
import Table, { defineData } from "@/components/Table";
import Tag from "@/components/Tag";
import { classifications } from "@/data/curations";
import { downloadJson } from "@/util/download";
import { IconDownload } from "@tabler/icons-react";

export default function LegendTable() {
  return (
    <>
      <Button
        design="plain"
        className="self-end"
        onClick={() => downloadJson(classifications, "classifications")}
      >
        <IconDownload />
        Download
      </Button>
      <Table
        className="w-full"
        {...defineData(Object.values(classifications), (column) => [
          column({
            key: "value",
            name: "Classification",
            render: (cell) => <Tag value={cell} className="w-full" />,
          }),
          column({
            key: "score",
            name: "Score",
          }),
          column({
            key: "description",
            name: "Description",
            className: "justify-start text-left",
          }),
        ])}
        showControls={false}
      />
    </>
  );
}
