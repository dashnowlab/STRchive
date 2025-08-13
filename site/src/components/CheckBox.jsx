import { useEffect, useRef } from "react";
import clsx from "clsx";
import classes from "./CheckBox.module.css";

/** checkbox with label */
const CheckBox = ({
  label,
  checked,
  onChange,
  tooltip,
  className,
  ...props
}) => {
  const ref = useRef(null);

  useEffect(() => {
    ref.current.indeterminate = checked === "mixed";
  }, [checked]);

  return (
    <label data-tooltip={tooltip}>
      <input
        ref={ref}
        type="checkbox"
        className={clsx(className, classes.box)}
        checked={!!checked}
        onChange={() => {
          if (checked === "mixed") onChange?.(true);
          else if (checked === true) onChange?.(false);
          else onChange?.("mixed");
        }}
        aria-checked={checked}
        aria-label={tooltip}
        {...props}
      />
      {label && <span>{label}</span>}
    </label>
  );
};

export default CheckBox;
