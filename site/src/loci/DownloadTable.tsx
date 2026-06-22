import Button from "@/components/Button";
import Link from "@/components/Link";
import Table, { defineData } from "@/components/Table";
import { repoRaw, version } from "@/layouts/meta";
import { filter, map, uniq } from "lodash-es";

/** import catalog files */
const catalogs = Object.keys(
  import.meta.glob(`~/catalogs/STRchive-disease-loci*.*`, {
    query: "url",
    import: "default",
    eager: true,
  }),
).map((path) => {
  /** last part of path */
  path = path.split("/").pop() ?? "";
  /** filename parts */
  const [, genome, software, ...extensions] = path.split(".");
  /** insert version in download filename */
  const filename = path.replace(/(STRchive-disease-loci)/, `$1-v${version}`);
  return {
    path,
    filename,
    genome,
    software,
    extension: ["", ...extensions].join("."),
  };
});

/** rows, by software */
const rows = [
  {
    key: "general",
    primary: "General",
    secondary:
      "a general-purpose extended bed file for filtering and annotating loci",
  },

  {
    key: "TRGT",
    primary: "TRGT",
    secondary: "for genotyping full allele sequences in PacBio HiFi reads",
    link: "https://github.com/PacificBiosciences/trgt",
  },
  {
    key: "atarva",
    primary: "Atarva",
    secondary: "⚠️ for genotyping full allele sequences in long-read data",
    link: "https://github.com/SowpatiLab/ATaRVa",
  },
  {
    key: "longTR",
    primary: "LongTR",
    secondary: "⚠️ for genotyping full allele sequences in long-read data",
    link: "https://github.com/gymrek-lab/LongTR",
  },
  {
    key: "straglr",
    primary: "Straglr",
    secondary: "⚠️ for genotyping allele sizes in long read-data",
    link: "https://github.com/bcgsc/straglr",
  },
  {
    key: "stranger",
    primary: "Stranger",
    secondary:
      "⚠️ for annotating TRGT or ExpansionHunter allele sizes with pathologic implications.",
    link: "https://github.com/Clinical-Genomics/stranger",
  },
];

export default function DownloadTable() {
  return (
    <Table
      showControls={false}
      {...defineData(rows, (column) => [
        column({
          key: "key",
          name: "",
          className: "flex-col items-start text-left",
          sortable: false,
          render: (value, row) => (
            <>
              {row.link ? (
                <Link to={row.link}>
                  <strong>{row.primary}</strong>
                </Link>
              ) : (
                <div>
                  <strong>{row.primary}</strong>
                </div>
              )}

              <div>{row.secondary}</div>
            </>
          ),
        }),
        ...uniq(map(catalogs, "genome")).map((genome) =>
          column({
            key: "key",
            name: genome,
            sortable: false,
            className: "flex-wrap gap-1",
            render: (value, row) =>
              filter(
                catalogs,
                (catalog) =>
                  catalog.software === row.key && catalog.genome === genome,
              ).map((download, index) => (
                <Button
                  key={index}
                  design="plain"
                  to={`${repoRaw}/refs/heads/main/data/catalogs/${download.path}`}
                  download={download.filename}
                  arrow={false}
                >
                  {download.extension}
                </Button>
              )),
          }),
        ),
      ])}
    />
  );
}
