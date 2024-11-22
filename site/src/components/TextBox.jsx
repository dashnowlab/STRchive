import { FaXmark } from "react-icons/fa6";

/** text input with label */
const TextBox = ({ label, onChange, ...props }) => (
  <label>
    {label}
    <input
      {...props}
      type="text"
      onChange={(event) => onChange?.(event.target.value)}
    />
    <button
      onClick={(event) => {
        event.currentTarget.previousElementSibling.value == "";
        onChange("");
      }}
    >
      <FaXmark />
    </button>
  </label>
);

export default TextBox;
