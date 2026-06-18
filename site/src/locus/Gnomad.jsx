import { useState } from "react";
import RangeChart from "@/components/RangeChart";
import Select from "@/components/Select";
import { startCase } from "lodash-es";

/** charts for gnomad data */
const Gnomad = ({ title, data = {} }) => {
  /** sex options */
  const sexes = [
    { value: "both", label: "Both" },
    { value: "XX", label: "XX" },
    { value: "XY", label: "XY" },
  ].filter(({ value }) => value in data);

  /** selected sex */
  const [sex, setSex] = useState(sexes[0].value);

  if (!sexes.length) return <></>;

  /** datum object */
  const d = data[sex];

  return (
    <>
      {/* sex dropdown */}
      {!(sexes.length === 1 && sexes[0].value === "both") && (
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
                <span>${label}</span>
              </p>
              <p>
                <span>Value:</span>
                <br/>
                <span>${value.toFixed(2)}%</span>
              </p>
              <p>
                <span>95% Confidence Interval:</span>
                <br/>
                <span>${lower.toFixed(2)}% &ndash; ${upper.toFixed(2)}%</span>
              </p>
              <p>
                <span>Count:</span>
                <br/>
                <span>${count.toLocaleString()}</span>
              </p>
            `;
          }}
        />
      </div>

      <p className="text-center text-balance">
        Pathogenic Genotype (%): % of individuals predicted to be affected based
        on their genotype
      </p>
    </>
  );
};

export default Gnomad;
