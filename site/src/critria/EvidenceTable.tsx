import type { Curation } from "@/data";
import { Fragment } from "react/jsx-runtime";
import Link from "@/components/Link";
import Table, { defineData } from "@/components/Table";
import Cited from "@/locus/Cited";

type Props = {
  name: string;
  evidence:
    | Curation["genetic_evidence_details"]
    | Curation["experimental_evidence_details"];
};

/** evidence table on individual critria page */
export default function EvidenceTable({ name, evidence }: Props) {
  return (
    <Table
      itemNames={name}
      {...defineData(evidence, (column) => [
        column({
          key: "evidence_category",
          name: "Category",
        }),
        column({
          key: "Evidence type",
          name: "Type",
        }),
        column({
          key: "Citation",
          name: "Citation",
          render: (cell) =>
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
        }),
        column({
          key: "Score",
          name: "Score",
        }),
        column({
          key: "Evidence detail",
          name: "Details",
          className: "min-w-100 justify-start text-left",
          render: (cell) => (
            <p>
              <Cited value={cell} />
            </p>
          ),
        }),
      ])}
    />
  );
}
