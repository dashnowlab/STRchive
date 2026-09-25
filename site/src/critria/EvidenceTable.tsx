import type { Curation } from "@/data";
import { Fragment } from "react/jsx-runtime";
import Link from "@/components/Link";
import Table from "@/components/Table";
import Cited from "@/locus/Cited";

type Props = {
  name: string;
  evidence:
    | Curation["genetic_evidence_details"]
    | Curation["experimental_evidence_details"];
};

/** evidence table on individual critria page */
export default function EvidenceTable({ name, evidence }: Props) {
  type Datum = Props["evidence"][number];

  return (
    <Table
      itemNames={name}
      rows={evidence}
      columns={[
        {
          key: "evidence_category",
          name: "Category",
        },
        {
          key: "Evidence type",
          name: "Type",
        },
        {
          key: "Citation",
          name: "Citation",
          render: (cell: Datum["Citation"]) =>
            [...cell.matchAll(/pmid:\s*(\d+)/gi)]
              .map((match) => match[1])
              .map((pmid, index, array) => (
                <Fragment key={pmid}>
                  <Link
                    to={`https://pubmed.ncbi.nlm.nih.gov/${pmid}`}
                    arrow={false}
                  >
                    {`PMID:${pmid}`}
                  </Link>
                  {index < array.length - 1 ? " " : ""}
                </Fragment>
              )),
        },
        {
          key: "Score",
          name: "Score",
        },
        {
          key: "Evidence detail",
          name: "Details",
          className: "min-w-100 justify-start text-left",
          render: (cell: Datum["Evidence detail"]) => (
            <p>
              <Cited value={cell} />
            </p>
          ),
          download: (cell: Datum["Evidence detail"]) =>
            cell?.map((part) => part.text).join(" ") ?? "",
        },
      ]}
    />
  );
}
