import { Fragment } from "react/jsx-runtime";
import Link from "@/components/Link";
import Table from "@/components/Table";

/** column definitions */
const cols = [
  { key: "evidence_category", name: "Evidence Category" },
  { key: "Evidence type", name: "Evidence Type" },
  {
    key: "Citation",
    name: "Citation",
    render: (cell) =>
      getPmidLinks(cell).map((pmid, index, array) => (
        <Fragment key={pmid}>
          <Link to={`https://pubmed.ncbi.nlm.nih.gov/${pmid}`} arrow={false}>
            {`PMID:${pmid}`}
          </Link>
          {index < array.length - 1 ? " " : ""}
        </Fragment>
      )),
  },
  { key: "Score", name: "Score" },
  { key: "Values", name: "Values" },
  { key: "Evidence detail", name: "Evidence Detail" },
  { key: "Notes", name: "Notes" },
];

const getPmidLinks = (citation) => {
  if (typeof citation !== "string") return [];
  return [...citation.matchAll(/pmid:\s*(\d+)/gi)].map((match) => match[1]);
};

/** evidence table on individual critria page */
const EvidenceTable = ({ evidence }) => {
  return <Table cols={cols} rows={evidence} />;
};

export default EvidenceTable;
