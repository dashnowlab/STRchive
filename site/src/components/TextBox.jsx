import { useRef } from "react";
import { FaXmark } from "react-icons/fa6";
import classes from "./TextBox.module.css";

/** text input with label */
const TextBox = ({ label, onChange, ...props }) => {
  const ref = useRef();

  return (
    <label>
      <span>{label}</span>
      <div className={classes.container}>
        <input
          ref={ref}
          {...props}
          type="text"
          onChange={(event) => onChange?.(event.target.value)}
        />

        <div className={classes.buttons}>
          <button
            type="button"
            onClick={() => {
              if (ref.current) ref.current.value = "";
              onChange?.("");
            }}
            aria-label="Clear textbox"
          >
            <FaXmark />
          </button>
        </div>
      </div>
    </label>
  );
};

export default TextBox;
