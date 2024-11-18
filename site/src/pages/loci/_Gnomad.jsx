import { useState } from "react";
import { startCase } from "lodash-es";
import RangeChart from "@/components/RangeChart";
import Select from "@/components/Select";

const Gnomad = ({ title, data = {} }) => {
  const options = Object.keys(data)
    .map((key) => ({
      value: key,
      label: startCase(key),
    }))
    .sort()
    .reverse();

  const [sex, setSex] = useState(options[0].value);

  const d = data[sex];

  if (!d) return <></>;

  return (
    <>
      {options.length > 1 && (
        <Select label="Sex" options={options} value={sex} onChange={setSex} />
      )}

      <div className="charts">
        <RangeChart
          title={`${title} (${sex.replace("_", " ")})`}
          xLabel="Pathogenic Genotype (%)"
          yLabels={d.labels}
          counts={d.counts}
          values={d.values}
          lowerBounds={d.confidence_lower_bounds}
          upperBounds={d.confidence_upper_bounds}
          tooltip={(index) => {
            const label = d.labels[index];
            const count = d.counts[index];
            const value = d.values[index];
            const lower = d.confidence_lower_bounds[index];
            const upper = d.confidence_upper_bounds[index];
            return `
              <p>
                <b>${label}</b>
              </p>
              <p>
                <span>Value:</span><br/>
                <b>${value.toFixed(2)}%</b>
              </p>
              <p>
                <span>95% Confidence Interval:</span><br/>
                <b>${lower.toFixed(2)}% &ndash; ${upper.toFixed(2)}%</b>
              </p>
              <p>
                <span>Count:</span><br/>
                <b>${count.toLocaleString()}</b>
              </p>
            `;
          }}
        />
      </div>
    </>
  );
};

export default Gnomad;
