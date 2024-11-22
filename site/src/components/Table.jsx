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
import { preserveScroll } from "@/util/dom";
import Select from "./Select";
import classes from "./Table.module.css";

/** options for per-page select */
const perPageOptions = [
  { value: 5, label: "5" },
  { value: 10, label: "10" },
  { value: 25, label: "25" },
  { value: 50, label: "50" },
  { value: 100, label: "100" },
  { value: 9999, label: "All" },
];

/** per-page option selected at start */
const defaultPerPage = perPageOptions.at(-1);

/** table component with sorting, filtering, and more */
const Table = ({ cols, rows, sort = undefined, showControls = true }) => {
  /** current per-page selection */
  const [perPage] = useState(defaultPerPage.value);

  /** column definitions */
  const columnHelper = createColumnHelper();
  const columns = cols.map((col, index) =>
    columnHelper.accessor((row) => row[col.key], {
      id: String(index),
      header: col.name,
      enableSorting: col.sortable ?? true,
      enableColumnFilter: true,
      enableGlobalFilter: true,
      meta: {
        filterType: col.filterType,
        attrs: col.attrs,
        style: col.style,
      },
      /** render func for cell */
      cell: ({ cell, row }) =>
        col.render
          ? col.render?.(cell.getValue(), row.original)
          : cell.getValue(),
    }),
  );

  /** get original col def */
  const getCol = (id) => cols[Number(id)];

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
        pageSize: perPage,
      },
    },
  });

  return (
    <div className={classes.root}>
      <div className={classes.scroll}>
        {/* table */}
        <table
          className={classes.table}
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
                      <button
                        type="button"
                        disabled={!header.column.getCanSort()}
                        className={classes.th}
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
                      </button>
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
                <td className={classes.empty} colSpan={cols.length}>
                  No Rows
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {/* controls */}
      {showControls && (
        <div className={classes.controls}>
          {/* pagination */}
          <div className={classes["control-row"]}>
            <button
              type="button"
              className={classes["page-button"]}
              onClick={(event) => {
                table.setPageIndex(0);
                preserveScroll(event.currentTarget);
              }}
              disabled={!table.getCanPreviousPage()}
              aria-label="First page"
            >
              <FaAnglesLeft />
            </button>
            <button
              type="button"
              className={classes["page-button"]}
              onClick={(event) => {
                table.previousPage();
                preserveScroll(event.currentTarget);
              }}
              disabled={!table.getCanPreviousPage()}
              aria-label="Previous page"
            >
              <FaAngleLeft />
            </button>
            <button
              type="button"
              className={classes["page-button"]}
              onClick={(event) => {
                const page = parseInt(window.prompt("Jump to page") || "");
                if (Number.isNaN(page)) return;
                table.setPageIndex(page);
                preserveScroll(event.currentTarget);
              }}
            >
              Page{" "}
              {(table.getState().pagination.pageIndex + 1).toLocaleString()} of{" "}
              {(table.getPageCount() || 1).toLocaleString()}
            </button>
            <button
              type="button"
              className={classes["page-button"]}
              onClick={(event) => {
                table.nextPage();
                preserveScroll(event.currentTarget);
              }}
              disabled={!table.getCanNextPage()}
              aria-label="Next page"
            >
              <FaAngleRight />
            </button>
            <button
              type="button"
              className={classes["page-button"]}
              onClick={(event) => {
                table.setPageIndex(table.getPageCount() - 1);
                preserveScroll(event.currentTarget);
              }}
              disabled={!table.getCanNextPage()}
              aria-label="Last page"
            >
              <FaAnglesRight />
            </button>
          </div>

          {/* row count */}
          <>{table.getRowCount().toLocaleString()} rows</>

          {/* per page */}
          <Select
            label="Per page"
            options={perPageOptions}
            defaultValue={defaultPerPage.value}
            onChange={(value, event) => {
              table.setPageSize(Number(value));
              preserveScroll(event.currentTarget);
            }}
          />
        </div>
      )}
    </div>
  );
};

export default Table;
