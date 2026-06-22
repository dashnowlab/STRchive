import { useEffect, useRef } from "react";
import { fitViewBox } from "@/util/dom";
import { lerp } from "@/util/math";
import { useElementSize } from "@reactuses/core";
import { maxBy } from "lodash-es";

type Props = {
  rows: { name: string; className: string; values: number[] }[];
  xAxis: string;
  fontSize?: number;
  min?: number;
  max?: number;
  className?: string;
};

/** basic 2d line chart */
export default function LineChart({
  rows,
  xAxis,
  fontSize = 16,
  min = undefined,
  max = undefined,
  className = "",
}: Props) {
  const ref = useRef<SVGSVGElement>(null);

  /** available size */
  let [availableWidth] = useElementSize(ref);
  availableWidth ||= 200;

  /** main chart area size */
  const width =
    availableWidth -
    /** estimate space used by y-axis labels */
    maxBy(rows, (row) => row.name.length)!.name.length * 0.6 * fontSize;
  const height = (rows.length + 0.5) * 2 * fontSize;

  /** rows with missing values filtered out */
  const filteredRows = rows.filter(({ values }) =>
    values.every((value) => value !== undefined && value !== null),
  );

  /** map value to svg units */
  const values = filteredRows.map((row) => row.values).flat();
  min ??= Math.min(...values);
  max ??= Math.max(...values);
  const scaleX = (value: number) =>
    lerp(value, min, max, 1.5 * fontSize, width - 1.5 * fontSize);
  const scaleY = (index: number) => (index + 0.75) * 2 * fontSize;

  /** rows with derived props */
  const mappedRows = rows
    /** don't show row if we're missing values */
    .filter(({ values }) =>
      values.every((value) => value !== undefined && value !== null),
    )
    .map((row, index) => {
      /** value x to svg coords */
      const x = row.values.map(scaleX);
      /** if same value, nudge apart to give bar nominal width */
      if (x[0] === x[1]) {
        x[0] -= fontSize * 0.1;
        x[1] += fontSize * 0.1;
      }
      /** value y to svg coords */
      const y = scaleY(index);
      return { ...row, x, y };
    });

  /** fit to contents */
  useEffect(() => {
    if (!ref.current) return;
    const { x, y, width, height } = fitViewBox(ref.current);
    ref.current.setAttribute("viewBox", [x, y, width, height].join(" "));
  });

  /** if no rows, don't render anything */
  if (!mappedRows.length) return <></>;

  return (
    <svg
      ref={ref}
      className={className}
      style={{
        fontSize: fontSize + "px",
      }}
    >
      {/* axis lines */}
      <path
        className="fill-none stroke-current"
        d={`M 0 0 L 0 ${height} L ${width} ${height}`}
      />

      {/* y-axis label */}
      <g textAnchor="end" dominantBaseline="central">
        {mappedRows.map(({ name, y }, index) => (
          <text x={-fontSize * 0.5} y={y} key={index}>
            {name}
          </text>
        ))}
      </g>

      {/* x-axis label */}
      <text
        x={width * 0.5}
        y={height + fontSize * 0.5}
        textAnchor="middle"
        dominantBaseline="hanging"
      >
        {xAxis}
      </text>

      {/* colored bars */}
      {mappedRows.map(({ x, y, className }, index) => (
        <rect
          key={index}
          className={className}
          x={x[0]}
          y={y - fontSize * 0.5}
          width={x[1] - x[0]}
          height={fontSize * 1}
        />
      ))}

      {/* bar labels */}
      {mappedRows.map(({ x, y, values: [min, max] }, index) => (
        <g key={index} dominantBaseline="central">
          <text x={x[0]} y={y} textAnchor="end">
            {min.toLocaleString()}&nbsp;
          </text>
          <text x={x[1]} y={y}>
            &nbsp;{max.toLocaleString()}
          </text>
        </g>
      ))}
    </svg>
  );
}
