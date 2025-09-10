import clsx from "clsx";
import { sortBy, uniq } from "lodash-es";
import classes from "./NumberBox.module.css";

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
      className={clsx(className, classes.input)}
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
