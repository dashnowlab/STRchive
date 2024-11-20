/** dropdown select with label */
const Select = ({ label, options, onChange, ...props }) => (
  <label>
    {label}
    <select
      {...props}
      onChange={(event) => onChange?.(event.target.value, event)}
    >
      {options.map(({ value, label }, index) => (
        <option key={index} value={value}>
          {label}
        </option>
      ))}
    </select>
  </label>
);

export default Select;
