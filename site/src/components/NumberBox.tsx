import type { ReactNode } from "react";
import clsx from "clsx";
import { sortBy, uniq } from "lodash-es";

type Props = {
  label?: ReactNode;
  value?: number;
  min?: number;
  max?: number;
  step?: number;
  onChange?: (value: number) => void;
  snapValues?: number[];
  className?: string;
  tooltip?: string;
};

/** number input with label */
export default function NumberBox({
  label,
  value,
  onChange,
  snapValues,
  className,
  tooltip,
  ...props
}: Props) {
  return (
    <label data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <input
        type="number"
        className={clsx(
          className,
          "field-sizing-content min-h-10 min-w-10 rounded-md bg-white px-4 py-2 leading-normal ring-2 ring-inset",
        )}
        value={typeof value === "number" && !Number.isNaN(value) ? value : ""}
        onChange={(event) => {
          const string = event.target.value;
          if (!string.trim() || Number.isNaN(Number(string))) return;
          else {
            let number = Number(string);
            if (value && snapValues?.length) {
              const sortedSnaps = sortBy(uniq(snapValues));
              if (number < value)
                number =
                  sortedSnaps.findLast((value) => value <= number) ?? number;
              if (number > value)
                number = sortedSnaps.find((value) => value >= number) ?? number;
            }
            onChange?.(number);
          }
        }}
        {...props}
      />
    </label>
  );
}
