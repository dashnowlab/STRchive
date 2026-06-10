import type { CSSProperties, HTMLAttributes, ReactNode } from "react";
import type { SortingState } from "@tanstack/react-table";
import { useState } from "react";
import {
  FaAngleLeft,
  FaAngleRight,
  FaAnglesLeft,
  FaAnglesRight,
  FaSort,
  FaSortDown,
  FaSortUp,
} from "react-icons/fa6";
import Button from "@/components/Button";
import Select from "@/components/Select";
import {
  createColumnHelper,
  flexRender,
  getCoreRowModel,
  getFacetedMinMaxValues,
  getFacetedRowModel,
  getFacetedUniqueValues,
  getFilteredRowModel,
  getPaginationRowModel,
  getSortedRowModel,
  useReactTable,
} from "@tanstack/react-table";
import clsx from "clsx";

type Props<Datum extends object> = {
  cols: _Col<Datum>[];
  rows: Datum[];
  sort?: SortingState;
  showControls?: boolean;
};

type Col<
  Datum extends object = object,
  Key extends keyof Datum = keyof Datum,
> = {
  /** key of row object to access as cell value */
  key: Key;
  /** label for header */
  name: string;
  /** is sortable (default true) */
  sortable?: boolean;
  /** cell attributes */
  attrs?: HTMLAttributes<HTMLTableCellElement>;
  /** cell style */
  style?: CSSProperties;
  /** custom render function for cell */
  render?: (cell: NoInfer<Datum[Key]>, row: Datum) => ReactNode;
};

/**
 * https://stackoverflow.com/questions/68274805/typescript-reference-type-of-property-by-other-property-of-same-object
 * https://github.com/vuejs/core/discussions/8851
 */
type _Col<Datum extends object> = {
  [Key in keyof Datum]: Col<Datum, Key>;
}[keyof Datum];

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
  cols,
  rows,
  sort,
  showControls = true,
}: Props<Datum>) {
  /** current per-page selection */
  const [perPage] = useState(defaultPerPage.value);

  /** column definitions */
  const columnHelper = createColumnHelper<Datum>();
  const columns = cols.map((col, index) =>
    columnHelper.accessor((row) => row[col.key], {
      id: String(index),
      header: col.name,
      enableSorting: col.sortable ?? true,
      enableColumnFilter: true,
      enableGlobalFilter: true,
      meta: {
        attrs: col.attrs,
        style: col.style,
      },
      /** render func for cell */
      cell: ({ cell, row }) =>
        col.render
          ? col.render?.(cell.getValue() as Datum[keyof Datum], row.original)
          : cell.getValue(),
    }),
  );

  /** get original col def */
  const getCol = (id: string) => cols[Number(id)];

  /** tanstack table api */
  const table = useReactTable({
    data: rows,
    columns,
    getCoreRowModel: getCoreRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
    getPaginationRowModel: getPaginationRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFacetedRowModel: getFacetedRowModel(),
    getFacetedUniqueValues: getFacetedUniqueValues(),
    getFacetedMinMaxValues: getFacetedMinMaxValues(),
    getColumnCanGlobalFilter: () => true,
    autoResetPageIndex: true,
    columnResizeMode: "onChange",
    initialState: {
      sorting: sort,
      pagination: {
        pageIndex: 0,
        pageSize: Number(perPage),
      },
    },
  });

  return (
    <div className="flex w-full flex-col items-center gap-2">
      <div className="w-full overflow-x-auto rounded-md shadow-md">
        {/* table */}
        <table
          className="w-full"
          aria-rowcount={table.getPrePaginationRowModel().rows.length}
          aria-colcount={cols.length}
        >
          {/* head */}
          <thead>
            {table.getHeaderGroups().map((headerGroup) => (
              <tr key={headerGroup.id}>
                {headerGroup.headers.map((header) => (
                  <th
                    key={header.id}
                    aria-colindex={Number(header.id) + 1}
                    style={getCol(header.column.id)?.style}
                    align="left"
                    {...getCol(header.column.id)?.attrs}
                  >
                    {header.isPlaceholder ? null : (
                      <Button
                        disabled={!header.column.getCanSort()}
                        className="size-full justify-start gap-1 border-none bg-transparent px-2 py-1 text-dark-gray hover:bg-transparent hover:text-primary"
                        data-active={
                          header.column.getIsSorted() ? "" : undefined
                        }
                        onClick={header.column.getToggleSortingHandler()}
                        aria-label="Sort this column"
                      >
                        {/* header label */}
                        <span>
                          {flexRender(
                            header.column.columnDef.header,
                            header.getContext(),
                          )}
                        </span>

                        {header.column.getIsSorted() === "asc" && <FaSortUp />}
                        {header.column.getIsSorted() === "desc" && (
                          <FaSortDown />
                        )}
                        {header.column.getIsSorted() === false && (
                          <FaSort style={{ opacity: 0.1 }} />
                        )}
                      </Button>
                    )}
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
                    table.getState().pagination.pageIndex *
                      table.getState().pagination.pageSize +
                    index +
                    1
                  }
                >
                  {row.getVisibleCells().map((cell) => (
                    <td
                      key={cell.id}
                      style={getCol(cell.column.id)?.style}
                      align="left"
                      {...getCol(cell.column.id)?.attrs}
                    >
                      {flexRender(
                        cell.column.columnDef.cell,
                        cell.getContext(),
                      )}
                    </td>
                  ))}
                </tr>
              ))
            ) : (
              <tr>
                <td
                  className="p-4 text-center text-dark-gray"
                  colSpan={cols.length}
                >
                  No Rows
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {/* controls */}
      {showControls && (
        <div
          className={clsx(
            "flex w-full flex-wrap items-center justify-between gap-x-4 gap-y-2 max-md:flex-col",
          )}
        >
          {/* pagination */}
          <div className="flex flex-wrap items-center justify-center">
            <Button
              onClick={() => table.setPageIndex(0)}
              disabled={!table.getCanPreviousPage()}
              aria-label="First page"
            >
              <FaAnglesLeft />
            </Button>
            <Button
              onClick={() => table.previousPage()}
              disabled={!table.getCanPreviousPage()}
              aria-label="Previous page"
            >
              <FaAngleLeft />
            </Button>
            <Button
              onClick={() => {
                const page = parseInt(window.prompt("Jump to page") || "");
                if (Number.isNaN(page)) return;
                table.setPageIndex(page);
              }}
            >
              Page{" "}
              {(table.getState().pagination.pageIndex + 1).toLocaleString()} of{" "}
              {(table.getPageCount() || 1).toLocaleString()}
            </Button>
            <Button
              onClick={() => table.nextPage()}
              disabled={!table.getCanNextPage()}
              aria-label="Next page"
            >
              <FaAngleRight />
            </Button>
            <Button
              onClick={() => table.setPageIndex(table.getPageCount() - 1)}
              disabled={!table.getCanNextPage()}
              aria-label="Last page"
            >
              <FaAnglesRight />
            </Button>
          </div>

          {/* row count */}
          <>{table.getRowCount().toLocaleString()} rows</>

          {/* per page */}
          <Select
            label="Per page"
            options={perPageOptions}
            defaultValue={defaultPerPage.value}
            onChange={(value) => table.setPageSize(Number(value))}
          />
        </div>
      )}
    </div>
  );
}
