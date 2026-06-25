import type { ComponentProps } from "react";
import { useEffect, useRef } from "react";
import Plotly from "plotly.js-dist";

type Props = {
  // eslint-disable-next-line -- trust input structure
  data: any;
  // eslint-disable-next-line -- trust input structure
  layout: any;
  // data: Plotly.Data[];
  // layout: Partial<Plotly.Layout>;
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
