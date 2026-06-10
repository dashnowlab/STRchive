import type { ComponentProps } from "react";
import { useEffect, useRef } from "react";
import Plotly from "plotly.js-dist";

type Props = {
  data: Plotly.Data[];
  layout: Plotly.Layout;
} & ComponentProps<"div">;

/** plotly chart that supports any plotly options */
export default function PlotlyChart({ data, layout, ...props }: Props) {
  const ref = useRef(null);

  useEffect(() => {
    if (!ref.current) return;
    Plotly.newPlot(ref.current, data, layout, { responsive: true });
  }, [data, layout]);

  return (
    <div
      ref={ref}
      className="w-full min-w-full rounded-md shadow-md"
      {...props}
    />
  );
}
