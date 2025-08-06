import clsx from "clsx";
import classes from "./NumberBox.module.css";

/** number input with label */
const NumberBox = ({ label, value, onChange, className, ...props }) => (
  <label>
    <span>{label}</span>
    <input
      type="number"
      className={clsx(className, classes.box)}
      value={typeof value === "number" && !Number.isNaN(value) ? value : ""}
      onChange={(event) => onChange?.(Number(event.target.value) || 0)}
      {...props}
    />
  </label>
);

export default NumberBox;
