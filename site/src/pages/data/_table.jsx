import { Fragment, useState } from "react";
import { BsStars } from "react-icons/bs";
import { FaCircleExclamation } from "react-icons/fa6";
import CheckBox from "@/components/CheckBox";
import Link from "@/components/Link";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TableComponent from "@/components/Table";
import { getValues } from "@/util/object";
import { capitalize } from "@/util/string";
import classes from "./_table.module.css";

const tagOptions = [
  {
    value: "conflicting",
    icon: (
      <FaCircleExclamation
        className="error"
        aria-label="Conflicting evidence"
      />
    ),
  },
  { value: "new", icon: <BsStars className="success" aria-label="New" /> },
];

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
    key: "tags",
    render: (cell) => {
      return tagOptions
        .filter(({ value }) => cell.includes(value))
        .map(({ icon }) => icon);
    },
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
    name: "Motif",
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

const Table = ({ data }) => {
  const [tags, setTags] = useState(Array(tagOptions.length).fill(false));
  const [motif, setMotif] = useState(
    Math.max(
      ...data
        .map((d) => d?.pathogenic_motif_reference_orientation.length || 0)
        .filter(Boolean),
    ),
  );

  const [inheritance, setInheritance] = useState("all");
  const [search, setSearch] = useState("");

  /** derive some props at load time */
  const derivedData = data.map((d) => ({
    ...d,
    tags: [
      ...(d.tags ?? []),
      new Date().getFullYear() - d.Year < 3 && "new",
    ].filter(Boolean),
  }));

  const filteredData = derivedData
    /** filter data */
    .filter(
      (d) =>
        /** free text search */
        getValues(d).join(" ").match(new RegExp(search.trim(), "i")) &&
        /** tags */
        (tags.filter(Boolean).length
          ? tagOptions
              .filter((_, index) => tags[index])
              .every(({ value }) => d.tags.includes(value))
          : true) &&
        /** motif */
        d.pathogenic_motif_reference_orientation.length <=
          (motif || Infinity) &&
        /** inheritance */
        (inheritance === "all" || d.Inheritance === inheritance),
    );

  const inheritanceOptions = [{ value: "all", label: "All" }].concat(
    /** unique inheritance values */
    [...new Set(derivedData.map((d) => d.Inheritance))]
      .filter(Boolean)
      .map((inheritance) => ({
        value: inheritance,
        label: inheritance,
      })),
  );

  return (
    <>
      <div className={classes.filters}>
        <input
          placeholder="Search"
          value={search}
          onChange={(event) => setSearch(event.target.value)}
        />
        <div className={classes["filter-row"]}>
          {tagOptions.map(({ value, icon }, index) => (
            <CheckBox
              key={index}
              label={
                <>
                  {icon}
                  {capitalize(value)}
                </>
              }
              value={tags[index]}
              onChange={(value) => {
                const newTags = [...tags];
                newTags[index] = value;
                setTags(newTags);
              }}
            />
          ))}
        </div>
        <NumberBox
          label="Motif max"
          value={motif}
          onChange={setMotif}
          min={1}
          max={100}
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

      <TableComponent cols={cols} rows={filteredData} />
    </>
  );
};

export default Table;
