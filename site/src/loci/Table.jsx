import { useEffect, useState } from "react";
import { FaArrowRight, FaCheck, FaXmark } from "react-icons/fa6";
import { LuDownload } from "react-icons/lu";
import clsx from "clsx";
import { countBy, map, max, min, pick, uniq } from "lodash-es";
import Button from "@/components/Button";
import CheckBox from "@/components/CheckBox";
import Link from "@/components/Link";
import NumberBox from "@/components/NumberBox";
import Popover from "@/components/Popover";
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
    name: "Motif (len)",
    render: (cell) => (
      <div className={classes.motif}>
        <div data-tooltip={cell.join(", ")} className={classes["motif-chars"]}>
          {cell.join(", ")}
        </div>
        <div>
          ({cell.map((motif) => motif.length.toLocaleString()).join(", ")})
        </div>
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
  /** find shortest/longest motif lengths */
  const motifLengths = loci
    .map(
      (d) =>
        d?.pathogenic_motif_reference_orientation.map(
          (motif) => motif.length,
        ) || 0,
    )
    .flat();
  const shortestMotif = min(motifLengths);
  const longestMotif = max(motifLengths);

  /** loci with extra derived props */
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

    /** if url tag not one of important tags (checkboxes), search for tag */
    if (!map(importantTagOptions, "value").includes(tag))
      setSearch(normalize(tag));
    else
      /** set checkboxes */
      setTags(
        map(importantTagOptions, "value").map((value) =>
          value === tag ? true : "mixed",
        ),
      );
  }, []);

  /** selected tag filters */
  const [tags, setTags] = useState(
    Array(importantTagOptions.length).fill("mixed"),
  );
  /** motif filter state */
  const [motifMin, setMotifMin] = useState(shortestMotif);
  const [motifMax, setMotifMax] = useState(longestMotif);
  /** inheritance filter state */
  const [inheritance, setInheritance] = useState("all");
  /** search filter state */
  const [search, setSearch] = useState("");

  const normalizedSearch = normalize(search);

  const filteredLoci = derivedLoci.filter(
    (d) =>
      /** filter by free text search visible columns */
      normalize(
        getValues(
          pick(d, [...map(cols, "key"), "locus_tags", "disease_tags"]),
        ).join(" "),
      ).includes(normalizedSearch) &&
      /** filter by tags */
      importantTagOptions.every(({ value }, index) => {
        const filter = tags[index];
        if (filter === undefined) return;
        /** if "mixed", keep locus whether it includes tag or not */
        if (filter === "mixed") return true;
        /** if true, keep locus if includes tag. if false, keep if doesn't. */ else
          return d.locus_tags.includes(value) === filter;
      }) &&
      /** filter by motif lengths */
      d.pathogenic_motif_reference_orientation.every(
        (m) => m.length >= motifMin && m.length <= motifMax,
      ) &&
      /** filter by inheritance type */
      (inheritance === "all" || d.inheritance.includes(inheritance)),
  );

  return (
    <div className={clsx("col", classes.table)}>
      {/* filters */}
      <div className="row">
        <div className="row">
          <Popover
            label="Tags"
            button={(() => {
              const counts = countBy(tags);
              if (counts.mixed === tags.length) return <>Any</>;

              return (
                <>
                  {counts.true && (
                    <>
                      {counts.true}
                      <FaCheck style={{ color: "var(--primary)" }} />
                    </>
                  )}
                  {counts.false && (
                    <>
                      {counts.false}
                      <FaXmark style={{ color: "var(--secondary)" }} />
                    </>
                  )}
                </>
              );
            })()}
          >
            {importantTagOptions.map(
              ({ Icon, label, color, tooltip }, index) => (
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
              ),
            )}
          </Popover>
        </div>

        <div className="row">
          <TextBox
            className={classes.search}
            placeholder="Search"
            value={search}
            onChange={setSearch}
          />

          <div className={classes["motif-length"]}>
            Motif length
            <NumberBox
              data-tooltip="Motif length min"
              value={motifMin}
              onChange={setMotifMin}
              min={shortestMotif}
              max={Math.min(longestMotif, motifMax)}
              snapValues={motifLengths}
            />
            &ndash;
            <NumberBox
              data-tooltip="Motif length max"
              value={motifMax}
              onChange={setMotifMax}
              min={Math.max(shortestMotif, motifMin)}
              max={longestMotif}
              snapValues={motifLengths}
            />
          </div>

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
            <span>Download</span>
            <LuDownload />
          </Button>
        </div>
      </div>

      {/* table */}
      <TableComponent cols={cols} rows={filteredLoci} />
    </div>
  );
};

export default Table;
