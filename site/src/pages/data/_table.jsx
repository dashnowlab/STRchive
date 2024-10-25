import { Fragment } from "react";
import TableComponent from "@/components/Table";
import Link from "@/components/Link";
import Sparkle from "@/assets/sparkle.svg?react";
import Error from "@/assets/error.svg?react";

const cols = [
  {
    key: "gene",
    render: (cell) => (
      <Link to={`/data/${cell}`} className="button">
        View
      </Link>
    ),
    sortable: false,
  },
  {
    render: (cell, row) => (
      <>
        {row.new && <Sparkle className="success" aria-label="New" />}
        {row.conflict && (
          <Error className="error" aria-label="Conflicting evidence" />
        )}
      </>
    ),
    sortable: false,
  },
  {
    key: "gene",
    name: "Gene",
  },
  {
    key: "disease_id",
    name: "Disease",
  },
  {
    key: "disease",
    name: "Description",
  },
  {
    key: "position_hg38",
    name: "Position hg38",
    render: (cell) => (
      <Link
        to={`https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${cell}`}
      >
        {cell}
      </Link>
    ),
  },
  {
    key: "pathogenic_motif_reference_orientation",
    name: "Ref. Orient.",
    render: (cell) =>
      cell.match(/.{1,5}/g).map((part, index) => (
        <Fragment key={index}>
          {part}
          <wbr />
        </Fragment>
      )),
    attrs: { width: "100px" },
  },
  {
    key: "Inheritance",
    name: "Inheritance",
  },
];

const Table = ({ data }) => (
  <TableComponent
    cols={cols}
    rows={data.map((d) => ({
      ...d,
      new: new Date().getFullYear() - d.Year < 3,
      conflict: !!d.details?.match(/conflict/i),
      position_hg38: `${d.chrom}:${d.start_hg38}-${d.stop_hg38}`,
    }))}
  />
);

export default Table;
