import clsx from "clsx";
import classes from "./CheckBox.module.css";

/** checkbox with label */
const CheckBox = ({ label, onChange, tooltip, className, ...props }) => (
  <label data-tooltip={tooltip}>
    <input
      type="checkbox"
      className={clsx(className, classes.box)}
      onChange={(event) => onChange?.(event.target.checked)}
      aria-label={tooltip}
      {...props}
    />
    <span>{label}</span>
  </label>
);

export default CheckBox;
