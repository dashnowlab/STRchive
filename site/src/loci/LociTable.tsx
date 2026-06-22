import { useEffect, useState } from "react";
import Button from "@/components/Button";
import CheckBox from "@/components/CheckBox";
import Link from "@/components/Link";
import NumberBox from "@/components/NumberBox";
import Popover from "@/components/Popover";
import Select from "@/components/Select";
import Table, { defineData } from "@/components/Table";
import Tag from "@/components/Tag";
import TextBox from "@/components/TextBox";
import { loci } from "@/data";
import { tagOptions } from "@/data/tags";
import { getValues } from "@/util/object";
import { IconArrowRight, IconCheck, IconX } from "@tabler/icons-react";
import { countBy, map, max, min, pick, uniq } from "lodash-es";

/** tags to show in table and filters */
const filterTags = tagOptions.filter((tag) => tag.filter);

/** normalize string for search comparison */
const normalize = (string: string) =>
  string
    .toLowerCase()
    .replaceAll(/[^A-Za-z0-9]/g, " ")
    .replaceAll(/\s+/g, " ")
    .trim();

/** table for main loci page */
export default function LociTable() {
  /** find shortest/longest motif lengths */
  const motifLengths = loci
    .map(
      (locus) =>
        locus?.pathogenic_motif_reference_orientation.map(
          (motif) => motif.length,
        ) || 0,
    )
    .flat();
  const shortestMotif = min(motifLengths) ?? 0;
  const longestMotif = max(motifLengths) ?? Infinity;

  /** options for inheritance filter */
  const inheritanceOptions = [{ value: "all", label: "All" }].concat(
    uniq(map(loci, "inheritance").flat()).map((inheritance) => ({
      value: inheritance,
      label: inheritance,
    })),
  );

  /** selected tag filters */
  const [tags, setTags] = useState(Array(filterTags.length).fill("mixed"));
  /** motif filter state */
  const [motifMin, setMotifMin] = useState(shortestMotif);
  const [motifMax, setMotifMax] = useState(longestMotif);
  /** inheritance filter state */
  const [inheritance, setInheritance] = useState("all");
  /** search filter state */
  const [search, setSearch] = useState("");

  useEffect(() => {
    /** tag from url param */
    const tag = new URL(window.location.href).searchParams.get("tag") ?? "";

    if (!tag) return;

    /** set states. don't do as default useState values due to window obj and url param being blank on first SSR. */

    /** if url tag not one of important tags (checkboxes), search for tag */
    // eslint-disable-next-line
    if (!map(filterTags, "value").includes(tag)) setSearch(normalize(tag));
    else
      /** set checkboxes */
      setTags(
        map(filterTags, "value").map((value) =>
          value === tag ? true : "mixed",
        ),
      );
  }, []);

  const normalizedSearch = normalize(search);

  /** access locus key as string[] */
  const includes = (locus: object, key: string, value: string) =>
    (
      [locus[key as keyof typeof locus]].flat().filter(Boolean) as string[]
    ).includes(value);

  const filteredLoci = loci
    .map((locus) => ({
      ...locus,
      /** for sorting */
      tag_sort: (["evidence", "locus_tags", "disease_tags"] as const).map(
        (key) =>
          filterTags.findIndex(({ value }) => includes(locus, key, value)),
      ),
    }))
    .filter(
      (locus) =>
        /** filter by free text search visible columns */
        normalize(
          getValues(
            pick(locus, [
              "gene",
              "disease_id",
              "disease",
              "position_base0_hg38",
              "pathogenic_motif_reference_orientation",
              "inheritance",
              "locus_tags",
              "disease_tags",
              "evidence",
            ]),
          ).join(" "),
        ).includes(normalizedSearch) &&
        /** filter by tags */
        filterTags.every(({ key, value }, index) => {
          const filter = tags[index];
          if (filter === undefined) return;
          /** if "mixed", keep locus whether it includes tag or not */
          if (filter === "mixed") return true;
          /** if true, keep locus if includes tag. if false, keep if doesn't. */ else
            return includes(locus, key, value) === filter;
        }) &&
        /** filter by motif lengths */
        locus.pathogenic_motif_reference_orientation.every(
          (motif) =>
            motif.length >= (motifMin ?? NaN) &&
            motif.length <= (motifMax ?? NaN),
        ) &&
        /** filter by inheritance type */
        (inheritance === "all" ||
          [locus.inheritance].flat().includes(inheritance)),
    );

  /** selected tags label */
  let selected = <></>;
  const counts = countBy(tags);
  if (counts.mixed === tags.length) selected = <>Any</>;
  else
    selected = (
      <>
        {counts.true && (
          <>
            {counts.true}
            <IconCheck className="text-primary" />
          </>
        )}
        {counts.false && (
          <>
            {counts.false}
            <IconX className="text-secondary" />
          </>
        )}
      </>
    );

  return (
    <>
      <div className="flex flex-wrap items-center justify-center gap-4 max-lg:flex-col">
        {/* filters */}
        <div className="box w-full">
          <strong className="flex justify-between gap-8">
            Filters
            <span className="font-regular">
              {filteredLoci.length === loci.length
                ? "all"
                : filteredLoci.length.toLocaleString()}{" "}
              results
            </span>
          </strong>

          <div className="flex flex-wrap items-center gap-x-8 gap-y-4">
            <label>
              Tags
              <Popover
                hover={false}
                content={
                  <div className="grid grid-cols-[auto_auto_auto] gap-2 *:contents">
                    {filterTags.map(({ value, label }, index) => (
                      <CheckBox
                        key={index}
                        label={
                          <>
                            <Tag value={value} small />
                            <span className="leading-normal">{label}</span>
                          </>
                        }
                        checked={tags[index]}
                        onChange={(value) => {
                          const newTags = [...tags];
                          newTags[index] = value;
                          setTags(newTags);
                        }}
                      />
                    ))}
                  </div>
                }
              >
                <Button design="plain">{selected}</Button>
              </Popover>
            </label>

            <TextBox
              className="min-w-30 grow"
              placeholder="Search"
              value={search}
              onChange={setSearch}
            />

            <div className="flex items-center gap-1">
              Motif length
              <NumberBox
                value={motifMin}
                onChange={(value) => setMotifMin(value ?? shortestMotif)}
                min={shortestMotif}
                max={Math.min(longestMotif, motifMax)}
                snapValues={motifLengths}
                aria-label="Motif length min"
              />
              &ndash;
              <NumberBox
                value={motifMax}
                onChange={(value) => setMotifMax(value ?? longestMotif)}
                min={Math.max(shortestMotif, motifMin)}
                max={longestMotif}
                snapValues={motifLengths}
                aria-label="Motif length max"
              />
            </div>

            <Select
              label="Inheritance"
              value={inheritance}
              onChange={setInheritance}
              options={inheritanceOptions}
            />
          </div>
        </div>
      </div>

      {/* table */}
      <Table
        itemNames={filteredLoci.length < loci.length ? "filtered loci" : "loci"}
        {...defineData(filteredLoci, (column) => [
          column({
            key: "id",
            render: (cell) => (
              <Popover content="Go to locus page" button={false}>
                <Button className="p-0!" design="bubble" to={`/loci/${cell}`}>
                  <IconArrowRight />
                </Button>
              </Popover>
            ),
            sortable: false,
          }),
          column({
            /** use number value so column sorted by that instead of alphabetically */
            key: "tag_sort",
            name: "Tags",
            className: "gap-1",
            render: (cell, row) =>
              tagOptions
                .filter(
                  ({ key, value, filter }) =>
                    filter && includes(row, key, value),
                )
                .map(({ value }, index) => (
                  <Tag key={index} value={value} small />
                )),
          }),
          column({
            key: "gene",
            name: "Gene",
          }),
          column({
            key: "disease_id",
            name: "Disease",
          }),
          column({
            key: "disease",
            name: "Description",
            className: "justify-start text-left",
          }),
          column({
            key: "position_base0_hg38",
            name: "Position hg38",
            render: (cell, row) => (
              <Link
                to={`https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${row.position_base1_hg38}`}
              >
                {cell}
              </Link>
            ),
          }),
          column({
            key: "pathogenic_motif_reference_orientation",
            name: "Motif (len)",
            render: (cell) => (
              <div className="flex flex-col items-center">
                <Popover
                  content={
                    <span className="wrap-anywhere">{cell.join(", ")}</span>
                  }
                  button={false}
                >
                  <span className="max-w-20 truncate">{cell.join(", ")}</span>
                </Popover>
                <div>
                  (
                  {cell
                    .map((motif) => motif.length.toLocaleString())
                    .join(", ")}
                  )
                </div>
              </div>
            ),
          }),
          column({
            key: "inheritance",
            name: "Inheritance",
            render: (cell) => cell?.join("/"),
          }),
        ])}
      />
    </>
  );
}
