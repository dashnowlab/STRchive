import { useState } from "react";
import { map, pick, startCase, uniq } from "lodash-es";
import CheckBox from "@/components/CheckBox";
import Link from "@/components/Link";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TableComponent from "@/components/Table";
import { getValues } from "@/util/object";
import { deriveDatum } from "./_derived";
import classes from "./_Table.module.css";
import { tagOptions } from "./_tags";

/** column definitions */
const cols = [
  {
    key: "id",
    render: (cell) => (
      <Link to={`/loci/${cell}`} className="button">
        View
      </Link>
    ),
    sortable: false,
  },
  {
    key: "locus_tags",
    name: "Tags",
    render: (cell) => {
      return tagOptions
        .filter(({ value }) => cell.includes(value))
        .map(({ Icon, color, tooltip }, index) => (
          <Icon key={index} style={{ color }} data-tooltip={tooltip} />
        ));
    },
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
    name: "Motif",
    render: (cell) => (
      <div
        data-tooltip={cell.join(", ")}
        style={{
          display: "inline-block",
          maxWidth: "80px",
          overflow: "hidden",
          whiteSpace: "nowrap",
          textOverflow: "ellipsis",
        }}
      >
        {cell.join(", ")}
      </div>
    ),
  },
  {
    key: "inheritance",
    name: "Inheritance",
    render: (cell) => cell?.join("/"),
  },
];

/** table for main loci page */
const Table = ({ data }) => {
  /** find longest motif length */
  const maxMotif = Math.max(
    ...data
      .map(
        (d) =>
          d?.pathogenic_motif_reference_orientation.map(
            (motif) => motif.length,
          ) || 0,
      )
      .flat()
      .filter(Boolean),
  );

  const derivedData = data.map(deriveDatum);

  /** options for inheritance filter */
  const inheritanceOptions = [{ value: "all", label: "All" }].concat(
    uniq(map(derivedData, "inheritance").flat()).map((inheritance) => ({
      value: inheritance,
      label: inheritance,
    })),
  );

  /** selected tag filters */
  const [tags, setTags] = useState(Array(tagOptions.length).fill(false));
  /** motif filter state */
  const [motif, setMotif] = useState(maxMotif);
  /** inheritance filter state */
  const [inheritance, setInheritance] = useState("all");
  /** search filter state */
  const [search, setSearch] = useState("");

  /** filter data */
  const filteredData = derivedData.filter(
    (d) =>
      /** free text search visible columns */
      getValues(pick(d, map(cols, "key")))
        .join(" ")
        .match(new RegExp(search.trim(), "i")) &&
      /** tags */
      (tags.filter(Boolean).length
        ? tagOptions
            .filter((_, index) => tags[index])
            .every(({ value }) => d.locus_tags.includes(value))
        : true) &&
      /** motif */
      d.pathogenic_motif_reference_orientation.every(
        (m) => m.length <= (motif || Infinity),
      ) &&
      /** inheritance */
      (inheritance === "all" || d.inheritance.includes(inheritance)),
  );

  return (
    <>
      {/* filters */}
      <div className={classes.filters}>
        <input
          placeholder="Search"
          value={search}
          onChange={(event) => setSearch(event.target.value)}
        />
        <div className={classes["filter-row"]}>
          {tagOptions.map(({ value, Icon, color, tooltip }, index) => (
            <CheckBox
              key={index}
              label={
                <>
                  <Icon style={{ color }} />
                  {startCase(value)}
                </>
              }
              value={tags[index]}
              onChange={(value) => {
                const newTags = [...tags];
                newTags[index] = value;
                setTags(newTags);
              }}
              tooltip={tooltip}
            />
          ))}
        </div>
        <NumberBox
          label="Motif max"
          value={motif}
          onChange={setMotif}
          min={1}
          max={maxMotif}
        />
        <div>
          <Select
            label="Inheritance"
            value={inheritance}
            onChange={setInheritance}
            options={inheritanceOptions}
          />
        </div>
      </div>

      {/* table */}
      <TableComponent cols={cols} rows={filteredData} />
    </>
  );
};

export default Table;
