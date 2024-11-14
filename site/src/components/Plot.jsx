import { useEffect, useRef } from "react";
import Plotly from "plotly.js-dist";

const Plot = ({ data, layout }) => {
  const ref = useRef();

  useEffect(() => {
    if (!ref.current) return;
    Plotly.newPlot(ref.current, { data, layout });
  }, []);

  return <div ref={ref}></div>;
};

export default Plot;
