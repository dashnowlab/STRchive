import Table from "@/components/Table";
import Tag from "@/components/Tag";
import { classifications } from "@/data/curations";

export default function LegendTable() {
  return (
    <Table
      rows={Object.values(classifications)}
      columns={[
        {
          key: "value",
          name: "Classification",
          render: (cell) => <Tag value={cell} className="w-full" />,
        },
        {
          key: "score",
          name: "Score",
        },
        {
          key: "description",
          name: "Description",
          className: "justify-start text-left",
        },
      ]}
      className="w-full"
      itemNames="classifications"
      pageControls={false}
    />
  );
}
