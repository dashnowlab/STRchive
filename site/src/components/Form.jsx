/** form wrapper around set of fields */
const Form = ({ onSubmit, ...props }) => (
  <form
    style={{ display: "contents" }}
    onSubmit={(event) => {
      /** prevent page navigation */
      event.preventDefault();

      /** get data from form */
      const form = event.currentTarget;
      const formData = new FormData(form);

      const data = {};

      /** transform form data into nicer format */
      for (const [key, value] of formData.entries()) {
        /** if we can parse as number, do it */
        if (
          typeof value === "string" &&
          value.trim() &&
          !Number.isNaN(Number(value))
        )
          data[key] = Number(value);
        else
          /** raw (string) value */
          data[key] = String(value);
      }

      /** call callback with data */
      onSubmit(data);
    }}
    {...props}
  />
);

export default Form;
