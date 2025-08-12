import { useEffect, useRef } from "react";
import clsx from "clsx";
import { fitViewBox } from "@/util/dom";
import { lerp } from "@/util/math";
import classes from "./LineChart.module.css";

/** basic 2d line chart */
const LineChart = ({
  rows,
  xAxis,
  width = 200,
  height = 50,
  fontSize = 10,
  min = undefined,
  max = undefined,
  className = "",
}) => {
  const ref = useRef(null);

  /** rows with missing values filtered out */
  const filteredRows = rows.filter(({ values }) =>
    values.every((value) => value !== undefined && value !== null),
  );

  /** map value to svg units */
  const values = filteredRows.map((row) => row.values).flat();
  min ??= Math.min(...values);
  max ??= Math.max(...values);
  const scale = (value) =>
    lerp(value, min, max, fontSize * 2, width - fontSize * 2);

  /** rows with derived props */
  const mappedRows = rows
    /** don't show row if we're missing values */
    .filter(({ values }) =>
      values.every((value) => value !== undefined && value !== null),
    )
    .map((row, index, rows) => {
      /** value x svg coords */
      const x = row.values.map(scale);
      /** if same value, nudge apart to give bar nominal width */
      if (x[0] === x[1]) {
        x[0] -= fontSize * 0.1;
        x[1] += fontSize * 0.1;
      }
      /** y svg coord */
      const y = height * ((index + 1) / (rows.length + 1));
      return { ...row, x, y };
    });

  /** fit to contents */
  useEffect(() => {
    if (!ref.current) return;
    const { height } = fitViewBox(ref.current);
    ref.current.style.height = 1.9 * height + "px";
  });

  /** if no rows, don't render anything */
  if (!mappedRows.length) return <></>;

  return (
    <svg
      ref={ref}
      className={clsx(className, classes.chart)}
      style={{
        fontSize: fontSize + "px",
      }}
    >
      {/* axis lines */}
      <path
        fill="none"
        stroke="var(--black)"
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
      {mappedRows.map(({ x, y, color, values: [min, max] }, index) => (
        <rect
          key={index}
          fill={color}
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
};

export default LineChart;
