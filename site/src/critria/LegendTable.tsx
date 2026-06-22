import Table, { defineData } from "@/components/Table";
import Tag from "@/components/Tag";
import { classifications } from "@/data/curations";

export default function LegendTable() {
  return (
    <Table
      className="w-full"
      itemNames="classifications"
      pageControls={false}
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
    />
  );
}
