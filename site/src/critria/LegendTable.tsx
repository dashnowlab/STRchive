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
          sortable: false,
          render: (cell) => <Tag value={cell} className="w-full" />,
        },
        {
          key: "score",
          name: "Score",
          sortable: false,
        },
        {
          key: "description",
          name: "Description",
          sortable: false,
          className: "justify-start text-left",
        },
      ]}
      className="w-full"
      itemNames="classifications"
      pageControls={false}
    />
  );
}
