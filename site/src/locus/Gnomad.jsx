import { useState } from "react";
import { startCase } from "lodash-es";
import RangeChart from "@/components/RangeChart";
import Select from "@/components/Select";

/** charts for gnomad data */
const Gnomad = ({ title, data = {} }) => {
  /** sex options */
  const sexes = Object.keys(data)
    .map((key) => ({
      value: key,
      label: startCase(key),
    }))
    .sort()
    .reverse();

  /** selected sex */
  const [sex, setSex] = useState(sexes[0].value);

  /** datum object */
  const d = data[sex];

  if (!d) return <></>;

  return (
    <>
      {/* sex dropdown */}
      {sexes.length > 1 && (
        <Select label="Sex" options={sexes} value={sex} onChange={setSex} />
      )}

      <div className="charts">
        <RangeChart
          title={
            sexes.length > 1 ? `${title} (${sex.replace("_", " ")})` : title
          }
          xLabel="Pathogenic Genotype (%)"
          yLabels={d.labels}
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

      <p className="center">
        <b>Pathogenic Genotype (%):</b> % of individuals predicted to be
        affected based on their genotype
      </p>
    </>
  );
};

export default Gnomad;
