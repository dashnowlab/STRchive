import type curations from "~/criTRia-curations.json";
import { Fragment } from "react/jsx-runtime";
import Link from "@/components/Link";
import Table, { defineData } from "@/components/Table";

type Curation = (typeof curations)[number];

type Props = {
  evidence:
    | Curation["genetic_evidence_details"]
    | Curation["experimental_evidence_details"];
};

/** evidence table on individual critria page */
export default function EvidenceTable({ evidence }: Props) {
  return (
    <Table
      {...defineData(evidence, (col) => [
        col({ key: "evidence_category", name: "Evidence Category" }),
        col({ key: "Evidence type", name: "Evidence Type" }),
        col({
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
        col({ key: "Score", name: "Score" }),
        col({ key: "Evidence detail", name: "Evidence Detail" }),
      ])}
    />
  );
}
