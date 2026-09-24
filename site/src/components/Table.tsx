import type { ReactNode } from "react";
import type { Tabular } from "@/util/download";
import type { NoInfer, SortingState } from "@tanstack/react-table";
import { useState } from "react";
import Button from "@/components/Button";
import Popover from "@/components/Popover";
import Select from "@/components/Select";
import { downloadCsv, downloadJson, downloadTsv } from "@/util/download";
import {
  IconArrowsSort,
  IconChevronLeft,
  IconChevronRight,
  IconChevronsLeft,
  IconChevronsRight,
  IconDownload,
  IconSortAscending,
  IconSortDescending,
} from "@tabler/icons-react";
import {
  createColumnHelper,
  createPaginatedRowModel,
  createSortedRowModel,
  metaHelper,
  rowPaginationFeature,
  rowSortingFeature,
  tableFeatures,
  useTable,
} from "@tanstack/react-table";
import clsx from "clsx";
import { clamp } from "lodash-es";

type Meta = Pick<Column, "className">;

const features = tableFeatures({
  rowSortingFeature,
  sortedRowModel: createSortedRowModel(),
  rowPaginationFeature,
  paginatedRowModel: createPaginatedRowModel(),
  columnMeta: metaHelper<Meta>(),
});

type Features = typeof features;

export type Column<
  Datum extends object = object,
  Key extends keyof Datum = keyof Datum,
> = {
  /** key of row object to access as cell value */
  key: Key;
  /** label for header */
  name?: ReactNode;
  /** is sortable (default true) */
  sortable?: boolean;
  /** class on cells */
  className?: string;
  /** custom render function for cell */
  render?: (cell: NoInfer<Datum[Key]>, row: Datum) => ReactNode;
};

/**
 * https://stackoverflow.com/questions/68274805/typescript-reference-type-of-property-by-other-property-of-same-object
 * https://github.com/vuejs/core/discussions/8851
 */
type _Column<Datum extends object> = {
  [Key in keyof Datum]: Column<Datum, Key extends string ? Key : never>;
}[keyof Datum];

export type Columns<Datum extends object> = _Column<Datum>[];

type Props<Datum extends object> = {
  columns: _Column<Datum>[];
  rows: Datum[];
  sort?: SortingState;
  pageControls?: boolean;
  actionControls?: boolean;
  itemNames?: string;
  className?: string;
};

/** options for per-page select */
const perPageOptions = [
  { value: "5", label: "5" },
  { value: "10", label: "10" },
  { value: "25", label: "25" },
  { value: "50", label: "50" },
  { value: "100", label: "100" },
  { value: "9999", label: "All" },
];

/** per-page option selected at start */
const defaultPerPage = perPageOptions.at(-1)!;

/** table component with sorting, filtering, and more */
export default function Table<Datum extends object>({
  columns,
  rows,
  sort,
  pageControls = true,
  actionControls = true,
  itemNames = "rows",
  className,
}: Props<Datum>) {
  /** current per-page selection */
  const [perPage] = useState(defaultPerPage.value);

  const columnHelper = createColumnHelper<Features, Datum>();

  const columnDefinitions = columnHelper.columns(
    columns.map((column, index) =>
      columnHelper.accessor((row) => row[column.key], {
        id: String(index),
        header: () => column.name,
        enableSorting: column.sortable ?? true,
        meta: {
          className: column.className,
        },
        /** render func for cell */
        cell: ({ cell, row }) => {
          const raw = cell.getValue();
          return column.render ? column.render(raw, row.original) : raw;
        },
      }),
    ),
  );

  /** tanstack table api */
  const table = useTable({
    features,
    data: rows,
    columns: columnDefinitions,
    autoResetPageIndex: true,
    initialState: {
      sorting: sort,
      pagination: {
        pageIndex: 0,
        pageSize: Number(perPage),
      },
    },
  });

  /** download data, in json form */
  const json = table.getPrePaginatedRowModel().rows.map((row) => row.original);

  /** download data, in tabular form */
  const tabular: Tabular = [
    columns.map((column) =>
      typeof column.name === "string" ? column.name : String(column.key),
    ),
    ...table.getPrePaginatedRowModel().rows.map((row) =>
      row.getAllCells().map((cell) => {
        const value = cell.getValue();
        if (["string", "number", "boolean"].includes(typeof value))
          return String(value);
        if (value === null || value === undefined) return "";
        if (Array.isArray(value)) return value.join(", ");
        if (typeof value === "object") return JSON.stringify(value);
        return "";
      }),
    ),
  ];

  return (
    <div className="flex max-w-max flex-col gap-4 self-center">
      {/* controls */}
      {(pageControls || actionControls) && (
        <div
          className={clsx(
            "flex flex-wrap items-center justify-between gap-x-8 gap-y-4 max-md:flex-col",
          )}
        >
          {pageControls && (
            <>
              {/* per page */}
              <Select
                label="Per page"
                options={perPageOptions}
                defaultValue={defaultPerPage.value}
                onChange={(value) => table.setPageSize(Number(value))}
              />

              {/* pagination */}
              <div className="flex flex-wrap items-center justify-center *:p-1 *:hover:text-primary">
                <Button
                  onClick={() => table.setPageIndex(0)}
                  aria-disabled={!table.getCanPreviousPage()}
                  aria-label="First page"
                >
                  <IconChevronsLeft />
                </Button>
                <Button
                  onClick={() => table.previousPage()}
                  aria-disabled={!table.getCanPreviousPage()}
                  aria-label="Previous page"
                >
                  <IconChevronLeft />
                </Button>
                <Button
                  onClick={() => {
                    const page = parseInt(window.prompt("Jump to page") || "");
                    if (Number.isNaN(page)) return;
                    table.setPageIndex(
                      clamp(page, 0, Math.ceil(table.getPageCount())) - 1,
                    );
                  }}
                >
                  {(() => {
                    const rows = table.getPrePaginatedRowModel().rows.length;
                    const { pageIndex, pageSize } = table.state.pagination;
                    return [
                      rows ? pageIndex * pageSize + 1 : 0,
                      "–",
                      Math.min(
                        (pageIndex + 1) * Math.min(pageSize, rows),
                        rows,
                      ),
                      "of",
                      rows,
                    ].join(" ");
                  })()}
                </Button>
                <Button
                  onClick={() => table.nextPage()}
                  aria-disabled={!table.getCanNextPage()}
                  aria-label="Next page"
                >
                  <IconChevronRight />
                </Button>
                <Button
                  onClick={() => table.setPageIndex(table.getPageCount() - 1)}
                  aria-disabled={!table.getCanNextPage()}
                  aria-label="Last page"
                >
                  <IconChevronsRight />
                </Button>
              </div>
            </>
          )}

          {/* actions */}
          {actionControls && (
            <div
              className={clsx(
                "flex flex-wrap items-center gap-4",
                !pageControls && "ml-auto",
              )}
            >
              <Popover
                content={
                  <div className="flex flex-col gap-2">
                    <Button
                      design="plain"
                      onClick={() => downloadJson(json, itemNames)}
                    >
                      JSON
                      <IconDownload />
                    </Button>
                    <Button
                      design="plain"
                      onClick={() => downloadCsv(tabular, itemNames)}
                    >
                      CSV
                      <IconDownload />
                    </Button>
                    <Button
                      design="plain"
                      onClick={() => downloadTsv(tabular, itemNames)}
                    >
                      TSV
                      <IconDownload />
                    </Button>
                  </div>
                }
              >
                <Button design="plain">Download</Button>
              </Popover>
            </div>
          )}
        </div>
      )}

      <div className="w-full overflow-x-auto rounded-md shadow-md">
        {/* table */}
        <table
          className={className}
          aria-rowcount={table.getPrePaginatedRowModel().rows.length}
          aria-colcount={columns.length}
        >
          {/* head */}
          <thead>
            {table.getHeaderGroups().map((headerGroup) => (
              <tr key={headerGroup.id}>
                {headerGroup.headers.map((header) => (
                  <th
                    key={header.id}
                    aria-colindex={Number(header.id) + 1}
                    className="p-0"
                  >
                    {/* wrapper */}
                    <div
                      className={clsx(
                        "flex items-center justify-center gap-2 p-2 text-center",
                        header.column.columnDef.meta?.className,
                      )}
                    >
                      {/* header label */}
                      <span>
                        <table.FlexRender header={header} />
                      </span>

                      {/* sort button */}
                      {header.column.getCanSort() && (
                        <Button
                          className={clsx(
                            "text-dark-gray hover:bg-transparent hover:text-primary",
                            header.column.getIsSorted() === false &&
                              "not-hover:opacity-25",
                          )}
                          data-active={
                            header.column.getIsSorted() ? "" : undefined
                          }
                          onClick={header.column.getToggleSortingHandler()}
                          aria-label="Sort this column"
                        >
                          {header.column.getIsSorted() === "asc" && (
                            <IconSortAscending />
                          )}
                          {header.column.getIsSorted() === "desc" && (
                            <IconSortDescending />
                          )}
                          {header.column.getIsSorted() === false && (
                            <IconArrowsSort />
                          )}
                        </Button>
                      )}
                    </div>
                  </th>
                ))}
              </tr>
            ))}
          </thead>

          {/* body */}
          <tbody>
            {table.getRowModel().rows.length ? (
              table.getRowModel().rows.map((row, index) => (
                <tr
                  key={row.id}
                  aria-rowindex={
                    table.state.pagination.pageIndex *
                      table.state.pagination.pageSize +
                    index +
                    1
                  }
                >
                  {row.getAllCells().map((cell) => (
                    <td key={cell.id} className="p-0">
                      {/* wrapper */}
                      <div
                        className={clsx(
                          "flex items-center justify-center p-2 text-center",
                          cell.column.columnDef.meta?.className,
                        )}
                      >
                        <table.FlexRender cell={cell} />
                      </div>
                    </td>
                  ))}
                </tr>
              ))
            ) : (
              <tr>
                <td
                  className="p-4 text-center text-dark-gray"
                  colSpan={columns.length}
                >
                  No Rows
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
}
