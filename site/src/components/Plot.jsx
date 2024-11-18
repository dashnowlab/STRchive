import { useEffect, useRef } from "react";
import Plotly from "plotly.js-dist";
import classes from "./Plot.module.css";

const Plot = ({ data, layout, ...props }) => {
  const ref = useRef();

  useEffect(() => {
    if (!ref.current) return;
    Plotly.newPlot(ref.current, {
      data,
      layout: {
        ...layout,
        margin: { l: 50, r: 50, b: 100, t: 50, pad: 4 },
      },
      config: { responsive: true },
    });
  }, []);

  return <div ref={ref} className={classes.plot} {...props}></div>;
};

export default Plot;
