import type { ReactNode } from "react";
import type { ColumnDef, RowData, SortingState } from "@tanstack/react-table";
import type { ValueOf } from "type-fest";
import { useState } from "react";
import {
  LuArrowDownWideNarrow,
  LuArrowUpDown,
  LuArrowUpNarrowWide,
  LuChevronLeft,
  LuChevronRight,
  LuChevronsLeft,
  LuChevronsRight,
} from "react-icons/lu";
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

declare module "@tanstack/react-table" {
  // eslint-disable-next-line
  interface ColumnMeta<TData extends RowData, TValue> {
    className?: string;
  }
}

/**
 * https://stackoverflow.com/questions/68274805/typescript-reference-type-of-property-by-other-property-of-same-object
 * https://github.com/vuejs/core/discussions/8851
 */
type _Column<Datum extends object> = {
  [Key in keyof Datum]: Column<Datum, Key>;
}[keyof Datum];

type Column<
  Datum extends object = object,
  Key extends keyof Datum = keyof Datum,
> = {
  /** key of row object to access as cell value */
  key: Key;
  /** label for header */
  name?: string;
  /** is sortable (default true) */
  sortable?: boolean;
  /** class on cells */
  className?: string;
  /** custom render function for cell */
  render?: (cell: Datum[Key], row: Datum) => ReactNode;
};

/** helper to define both rows and columns on table in type-safe way */
export const defineData = <Datum extends object>(
  rows: Datum[],
  columnFunc: (
    column: <Key extends keyof Datum>(
      column: Column<Datum, Key>,
    ) => Column<Datum, Key>,
  ) => _Column<Datum>[],
) => {
  const helper = createColumnHelper<Datum>();

  const columns = columnFunc(
    <Key extends keyof Datum>(column: Column<Datum, Key>) => column,
  ).map((column, index) =>
    helper.accessor((row) => row[column.key], {
      id: String(index),
      header: column.name,
      enableSorting: column.sortable ?? true,
      enableColumnFilter: true,
      enableGlobalFilter: true,
      meta: {
        className: column.className,
      },
      /** render func for cell */
      cell: ({ cell, row }) =>
        column.render
          ? column.render(cell.getValue() as ValueOf<Datum>, row.original)
          : cell.getValue(),
    }),
  );

  return { rows, columns };
};

type Props<Datum extends object> = {
  columns: ColumnDef<Datum, Datum[keyof Datum]>[];
  rows: Datum[];
  sort?: SortingState;
  showControls?: boolean;
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
  showControls = true,
}: Props<Datum>) {
  /** current per-page selection */
  const [perPage] = useState(defaultPerPage.value);

  /** tanstack table api */
  // eslint-disable-next-line
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
                        "flex items-center justify-center p-2 text-center",
                        header.column.columnDef.meta?.className,
                      )}
                    >
                      {/* header label */}
                      <span>
                        {flexRender(
                          header.column.columnDef.header,
                          header.getContext(),
                        )}
                      </span>

                      {/* sort button */}
                      {header.column.getCanSort() && (
                        <Button
                          className="p-0! text-dark-gray hover:bg-transparent hover:text-primary"
                          data-active={
                            header.column.getIsSorted() ? "" : undefined
                          }
                          onClick={header.column.getToggleSortingHandler()}
                          aria-label="Sort this column"
                        >
                          {header.column.getIsSorted() === "asc" && (
                            <LuArrowUpNarrowWide />
                          )}
                          {header.column.getIsSorted() === "desc" && (
                            <LuArrowDownWideNarrow />
                          )}
                          {header.column.getIsSorted() === false && (
                            <LuArrowUpDown className="opacity-25" />
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
                    table.getState().pagination.pageIndex *
                      table.getState().pagination.pageSize +
                    index +
                    1
                  }
                >
                  {row.getVisibleCells().map((cell) => (
                    <td key={cell.id} className="p-0">
                      {/* wrapper */}
                      <div
                        className={clsx(
                          "flex items-center justify-center p-2 text-center",
                          cell.column.columnDef.meta?.className,
                        )}
                      >
                        {flexRender(
                          cell.column.columnDef.cell,
                          cell.getContext(),
                        )}
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
              aria-disabled={!table.getCanPreviousPage()}
              aria-label="First page"
            >
              <LuChevronsLeft />
            </Button>
            <Button
              onClick={() => table.previousPage()}
              aria-disabled={!table.getCanPreviousPage()}
              aria-label="Previous page"
            >
              <LuChevronLeft />
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
              aria-disabled={!table.getCanNextPage()}
              aria-label="Next page"
            >
              <LuChevronRight />
            </Button>
            <Button
              onClick={() => table.setPageIndex(table.getPageCount() - 1)}
              aria-disabled={!table.getCanNextPage()}
              aria-label="Last page"
            >
              <LuChevronsRight />
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
