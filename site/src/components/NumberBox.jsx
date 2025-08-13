import clsx from "clsx";
import { sortBy, uniq } from "lodash-es";
import classes from "./NumberBox.module.css";

/** number input with label */
const NumberBox = ({
  label,
  value,
  onChange,
  className,
  snapValues,
  ...props
}) => (
  <label>
    {label && <span>{label}</span>}
    <input
      type="number"
      className={clsx(className, classes.box)}
      value={typeof value === "number" && !Number.isNaN(value) ? value : ""}
      onChange={(event) => {
        let newValue = Number(event.target.value) || 0;
        if (snapValues?.length) {
          snapValues = sortBy(uniq(snapValues));
          if (newValue < value)
            newValue = snapValues.findLast((v) => v <= newValue);
          if (newValue > value)
            newValue = snapValues.find((v) => v >= newValue);
        }
        onChange?.(newValue);
      }}
      {...props}
    />
  </label>
);

export default NumberBox;
