/** number input with label */
const NumberBox = ({ label, value, onChange, ...props }) => (
  <label>
    {label}
    <input
      {...props}
      type="number"
      value={typeof value === "number" && !Number.isNaN(value) ? value : ""}
      onChange={(event) => onChange?.(Number(event.target.value) || 0)}
    />
  </label>
);

export default NumberBox;
