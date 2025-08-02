/** checkbox with label */
const CheckBox = ({ label, onChange, tooltip, ...props }) => (
  <label data-tooltip={tooltip}>
    <input
      {...props}
      type="checkbox"
      onChange={(event) => onChange?.(event.target.checked)}
      aria-label={tooltip}
    />
    <span>{label}</span>
  </label>
);

export default CheckBox;
