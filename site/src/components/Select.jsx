import { FaAngleDown } from "react-icons/fa";
import clsx from "clsx";
import { preserveScroll } from "@/util/dom";
import classes from "./Select.module.css";

/** dropdown select with label */
const Select = ({ label, options, onChange, className, tooltip, ...props }) => (
  <label data-tooltip={tooltip}>
    {label && <span>{label}</span>}
    <div className={classes.container}>
      <select
        className={clsx(className, classes.select)}
        onChange={(event) => {
          onChange?.(event.target.value);
          preserveScroll(event.target);
        }}
        {...props}
      >
        {options.map(({ value, label }, index) => (
          <option key={index} value={value}>
            {label}
          </option>
        ))}
      </select>
      <FaAngleDown className={classes.icon} />
    </div>
  </label>
);

export default Select;
