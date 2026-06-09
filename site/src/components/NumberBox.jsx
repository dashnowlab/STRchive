import clsx from "clsx";
import { sortBy, uniq } from "lodash-es";

/** number input with label */
const NumberBox = ({
  label,
  value,
  onChange,
  snapValues,
  className,
  tooltip,
  ...props
}) => (
  <label data-tooltip={tooltip}>
    {label && <span>{label}</span>}
    <input
      type="number"
      className={clsx(
        className,
        "min-h-10 rounded-md bg-white px-[15px] py-[7.5px] [box-shadow:inset_0_0_0_2px_var(--color-dark-gray)] [field-sizing:content]",
      )}
      value={typeof value === "number" && !Number.isNaN(value) ? value : ""}
      onChange={(event) => {
        let newValue = event.target.value;
        if (
          !newValue.trim() ||
          newValue === null ||
          Number.isNaN(Number(newValue))
        )
          onChange?.(null);
        else {
          newValue = Number(newValue);
          if (value !== undefined && snapValues?.length) {
            snapValues = sortBy(uniq(snapValues));
            if (newValue < value)
              newValue = snapValues.findLast((v) => v <= newValue);
            if (newValue > value)
              newValue = snapValues.find((v) => v >= newValue);
          }
          onChange?.(newValue);
        }
      }}
      {...props}
    />
  </label>
);

export default NumberBox;
