import { useRef } from "react";
import { FaXmark } from "react-icons/fa6";
import clsx from "clsx";
import Button from "@/components/Button";
import classes from "./TextBox.module.css";

/** text input with label */
const TextBox = ({ label, multi, onChange, className, tooltip, ...props }) => {
  const ref = useRef(null);

  const Component = multi ? "textarea" : "input";

  return (
    <label data-tooltip={tooltip}>
      {label && <span>{label}</span>}
      <div className={classes.container}>
        <Component
          ref={ref}
          type="text"
          className={clsx(className, classes.input)}
          onChange={(event) => onChange?.(event.target.value)}
          {...props}
        />

        <div className={classes.buttons}>
          <Button
            onClick={() => {
              if (ref.current) ref.current.value = "";
              onChange?.("");
            }}
            aria-label="Clear textbox"
          >
            <FaXmark />
          </Button>
        </div>
      </div>
    </label>
  );
};

export default TextBox;
