import type { ValueOf } from "type-fest";
import type gnomad from "~/plots/gnomad.json";
import { useState } from "react";
import RangeChart from "@/components/RangeChart";
import Select from "@/components/Select";

type Props = {
  title: string;
  data: ValueOf<typeof gnomad>;
};

/** charts for gnomad data */
export default function Gnomad({ title, data }: Props) {
  /** sex options */
  const sexes = [
    { value: "both", label: "Both" },
    { value: "XX", label: "XX" },
    { value: "XY", label: "XY" },
  ].filter(({ value }) => value in data);

  /** selected sex */
  const [sex, setSex] = useState("both");

  if (!sexes.length) return <></>;

  /** datum object */
  const datum = data[sex as keyof typeof data];

  if (typeof datum !== "object" || datum === null) return <></>;

  const {
    labels,
    values,
    confidence_lower_bounds,
    confidence_upper_bounds,
    counts,
  } = datum;

  return (
    <>
      {/* sex dropdown */}
      {!(sexes.length === 1 && sexes[0].value === "both") && (
        <Select label="Sex" options={sexes} value={sex} onChange={setSex} />
      )}

      <div className="grid grid-cols-2 gap-8 max-md:grid-cols-1">
        <RangeChart
          title={
            sexes.length > 1 ? `${title} (${sex.replace("_", " ")})` : title
          }
          xLabel="Pathogenic Genotype (%)"
          yLabels={labels}
          values={values}
          lowerBounds={confidence_lower_bounds}
          upperBounds={confidence_upper_bounds}
          tooltip={(index) => {
            const label = labels[index];
            const count = counts[index];
            const value = values[index];
            const lower = confidence_lower_bounds[index];
            const upper = confidence_upper_bounds[index];
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
}
