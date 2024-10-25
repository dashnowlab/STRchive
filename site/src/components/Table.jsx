import SortUp from "@/assets/sort-up.svg?react";
import SortDown from "@/assets/sort-down.svg?react";
import Sort from "@/assets/sort.svg?react";
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
import AngleRight from "@/assets/angle-right.svg?react";
import AngleLeft from "@/assets/angle-left.svg?react";
import AnglesRight from "@/assets/angles-right.svg?react";
import AnglesLeft from "@/assets/angles-left.svg?react";
import classes from "./Table.module.css";
import { useState } from "react";

const perPageOptions = [5, 10, 20, 50, 100];
const defaultPerPage = perPageOptions[0];

const Table = ({ cols, rows, sort = undefined }) => {
  const [perPage, setPerPage] = useState(defaultPerPage);
  const [search, setSearch] = useState("");

  const columnHelper = createColumnHelper();
  /** column definitions */
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
      sorting: sort ?? [{ id: "2", desc: false }],
      pagination: {
        pageIndex: 0,
        pageSize: perPage,
      },
    },
    state: {
      globalFilter: search,
    },
  });

  return (
    <>
      {/* controls */}
      <div className={classes.controls}>
        {/* per page */}
        <label className={classes["control-row"]}>
          Per page
          <select
            defaultValue={defaultPerPage}
            onChange={(event) => table.setPageSize(Number(event.target.value))}
          >
            {perPageOptions.map((value, index) => (
              <option key={index} value={value}>
                {value}
              </option>
            ))}
          </select>
        </label>

        {/* search */}
        <div className={classes["control-row"]}>
          <input
            placeholder="Search"
            value={search}
            onChange={(event) => setSearch(event.target.value)}
          />
          {table.getRowCount().toLocaleString()} results
        </div>

        {/* pagination */}
        <div className={classes["control-row"]}>
          <button
            type="button"
            className={classes["page-button"]}
            onClick={() => table.setPageIndex(0)}
            disabled={!table.getCanPreviousPage()}
            aria-label="First page"
          >
            <AnglesLeft />
          </button>
          <button
            type="button"
            className={classes["page-button"]}
            onClick={table.previousPage}
            disabled={!table.getCanPreviousPage()}
            aria-label="Previous page"
          >
            <AngleLeft />
          </button>
          <button
            type="button"
            className={classes["page-button"]}
            onClick={() => {
              const page = parseInt(window.prompt("Jump to page") || "");
              if (Number.isNaN(page)) return;
              table.setPageIndex(page);
            }}
          >
            Page {(table.getState().pagination.pageIndex + 1).toLocaleString()}{" "}
            of {table.getPageCount().toLocaleString()}
          </button>
          <button
            type="button"
            className={classes["page-button"]}
            onClick={table.nextPage}
            disabled={!table.getCanNextPage()}
            aria-label="Next page"
          >
            <AngleRight />
          </button>
          <button
            type="button"
            className={classes["page-button"]}
            onClick={() => table.setPageIndex(table.getPageCount() - 1)}
            disabled={!table.getCanNextPage()}
            aria-label="Last page"
          >
            <AnglesRight />
          </button>
        </div>
      </div>

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

                        {header.column.getIsSorted() === "asc" && <SortUp />}
                        {header.column.getIsSorted() === "desc" && <SortDown />}
                        {header.column.getIsSorted() === false && (
                          <Sort style={{ opacity: 0 }} />
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
    </>
  );
};

export default Table;
