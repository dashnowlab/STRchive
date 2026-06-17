import type curations from "~/criTRia-curations.json";
import { Fragment } from "react/jsx-runtime";
import Link from "@/components/Link";
import Table, { defineData } from "@/components/Table";

type Curations = typeof curations;
type Curation = Curations[number];

type Props = {
  evidence:
    | Curation["genetic_evidence_details"]
    | Curation["experimental_evidence_details"];
};

/** evidence table on individual critria page */
export default function EvidenceTable({ evidence }: Props) {
  return (
    <Table
      {...defineData(evidence, (column) => [
        column({
          key: "evidence_category",
          name: "Evidence Category",
        }),
        column({
          key: "Evidence type",
          name: "Evidence Type",
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
          name: "Evidence Detail",
        }),
      ])}
    />
  );
}
