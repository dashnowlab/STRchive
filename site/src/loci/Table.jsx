import { useEffect, useState } from "react";
import { FaArrowRight } from "react-icons/fa";
import { LuDownload } from "react-icons/lu";
import clsx from "clsx";
import { map, pick, uniq } from "lodash-es";
import Button from "@/components/Button";
import CheckBox from "@/components/CheckBox";
import Link from "@/components/Link";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TableComponent from "@/components/Table";
import TextBox from "@/components/TextBox";
import { deriveLocus } from "@/data/derived";
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
      <Button
        to={`/loci/${cell}`}
        design="bubble"
        data-tooltip="Go to locus page"
      >
        <FaArrowRight />
      </Button>
    ),
    sortable: false,
  },
  {
    key: "locus_tags",
    name: "Tags",
    render: (cell) => (
      <div className={clsx("row", classes["tags-cell"])}>
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
    key: "position_base0_hg38",
    name: "Position hg38",
    render: (cell, row) => (
      <Link
        to={`https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${row.position_base1_hg38}`}
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
const Table = ({ loci }) => {
  /** find longest motif length */
  const maxMotif = Math.max(
    ...loci
      .map(
        (d) =>
          d?.pathogenic_motif_reference_orientation.map(
            (motif) => motif.length,
          ) || 0,
      )
      .flat()
      .filter(Boolean),
  );

  const derivedLoci = loci.map((locus) => deriveLocus(locus, loci));

  /** options for inheritance filter */
  const inheritanceOptions = [{ value: "all", label: "All" }].concat(
    uniq(map(derivedLoci, "inheritance").flat()).map((inheritance) => ({
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

  /** filter loci */
  const filteredLoci = derivedLoci.filter(
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
    <div className="col">
      {/* filters */}
      <div className="row">
        <div className="row">
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

        <div className="row">
          <TextBox placeholder="Search" value={search} onChange={setSearch} />

          <NumberBox
            label="Motif max"
            value={motif}
            onChange={setMotif}
            min={1}
            max={maxMotif}
          />

          <Select
            label="Inheritance"
            value={inheritance}
            onChange={setInheritance}
            options={inheritanceOptions}
          />
        </div>

        {/* row count */}
        <div className="row">
          <strong>{filteredLoci.length.toLocaleString()} loci</strong>
          <Button
            design="plain"
            onClick={() =>
              /** download filtered loci */
              downloadJson(filteredLoci, [
                filteredLoci.length < derivedLoci.length ? "filtered" : "",
                "loci",
              ])
            }
            data-tooltip="Download filtered loci"
          >
            Download <LuDownload />
          </Button>
        </div>
      </div>

      {/* table */}
      <TableComponent cols={cols} rows={filteredLoci} />
    </div>
  );
};

export default Table;
