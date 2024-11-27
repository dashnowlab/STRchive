import { useEffect, useState } from "react";
import { LuDownload } from "react-icons/lu";
import { map, pick, uniq } from "lodash-es";
import CheckBox from "@/components/CheckBox";
import Link from "@/components/Link";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TableComponent from "@/components/Table";
import TextBox from "@/components/TextBox";
import { deriveDatum } from "@/data/derived";
import { tagOptions } from "@/data/tags";
import { downloadJson } from "@/util/download";
import { getValues } from "@/util/object";
import classes from "./Table.module.css";

/** tags considered important enough to show in table and filters */
const importantTagOptions = tagOptions.filter((tag) => tag.important);

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
    render: (cell) => (
      <div className={classes["tags-cell"]}>
        {importantTagOptions
          .filter(({ value }) => cell.includes(value))
          .map(({ Icon, color, tooltip }, index) => (
            <Icon key={index} style={{ color }} data-tooltip={tooltip} />
          ))}
      </div>
    ),
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

/** normalize string for search comparison */
const normalize = (string) =>
  string
    .toLowerCase()
    .replaceAll(/[^A-Za-z0-9]/g, " ")
    .replaceAll(/\s+/g, " ")
    .trim();

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

  useEffect(() => {
    /** tag from url param */
    const tag = new URL(location).searchParams.get("tag") ?? "";

    if (!tag) return;

    /** set states. don't do as default useState values due to window obj and url param being blank on first SSR. */

    const checked = map(importantTagOptions, "value").map(
      (value) => value === tag,
    );
    if (checked.some(Boolean)) {
      /** set checkboxes */
      setTags(checked);
    } else {
      /** set search */
      setSearch(normalize(tag));
    }
  }, []);

  /** selected tag filters */
  const [tags, setTags] = useState(
    Array(importantTagOptions.length).fill(false),
  );
  /** motif filter state */
  const [motif, setMotif] = useState(maxMotif);
  /** inheritance filter state */
  const [inheritance, setInheritance] = useState("all");
  /** search filter state */
  const [search, setSearch] = useState("");

  const normalizedSearch = normalize(search);

  /** filter data */
  const filteredData = derivedData.filter(
    (d) =>
      /** free text search visible columns */
      normalize(
        getValues(
          pick(d, [...map(cols, "key"), "locus_tags", "disease_tags"]),
        ).join(" "),
      ).includes(normalizedSearch) &&
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
        <TextBox placeholder="Search" value={search} onChange={setSearch} />
        <div className={classes["filter-row"]}>
          {importantTagOptions.map(({ Icon, label, color, tooltip }, index) => (
            <CheckBox
              key={index}
              label={
                <>
                  <Icon style={{ color }} />
                  {label}
                </>
              }
              checked={tags[index]}
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

      {/* row count */}
      <div className={classes["filter-row"]}>
        <strong>{filteredData.length.toLocaleString()} loci</strong>
        <button
          className={classes.download}
          onClick={() =>
            /** download filtered data */
            downloadJson(filteredData, [
              filteredData.length < derivedData.length ? "filtered" : "",
              "loci",
            ])
          }
          data-tooltip="Download filtered loci"
        >
          Download <LuDownload />
        </button>
      </div>

      {/* table */}
      <TableComponent cols={cols} rows={filteredData} />
    </>
  );
};

export default Table;
