const CheckBox = ({ label, value, onChange, ...props }) => (
  <label>
    <input
      {...props}
      type="checkbox"
      onChange={(event) => onChange?.(event.target.checked)}
    />
    {label}
  </label>
);

export default CheckBox;
