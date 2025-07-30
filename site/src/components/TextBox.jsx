import { FaXmark } from "react-icons/fa6";
import classes from "./TextBox.module.css";

/** text input with label */
const TextBox = ({ label, onChange, ...props }) => (
  <label>
    {label}
    <div className={classes.container}>
      <input
        {...props}
        type="text"
        onChange={(event) => onChange?.(event.target.value)}
      />

      <div className={classes.buttons}>
        <button
          onClick={(event) => {
            event.currentTarget.previousElementSibling.value == "";
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

export default TextBox;
